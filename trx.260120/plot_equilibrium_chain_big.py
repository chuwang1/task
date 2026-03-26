#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_sumavir2_from_gfile import read_gfile_2d, read_input_profiles_mapping, bilinear_interp
from plot_tts_f_relation import read_gfile_fpol

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_OUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
GFILE = os.path.join(OMFIT_OUT, "g260206.20000_teq_0114")
INPUT_PROFILES = os.path.join(OMFIT_OUT, "input.profiles")


def read_tr_eqcalq(path):
    lines = open(path).read().splitlines()
    out = {"rho": [], "psip": [], "psit": [], "dvdpsit": [], "dvdpsip": [], "qps": [], "tts": [], "sumavir2": []}
    for line in lines:
        if not line.strip() or line.startswith('#') or line.startswith('nr,'):
            continue
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["rho"].append(vals[0])
        out["psip"].append(vals[1])
        out["psit"].append(vals[2])
        out["dvdpsit"].append(vals[3])
        out["dvdpsip"].append(vals[4])
        out["qps"].append(vals[5])
        out["tts"].append(vals[6])
        out["sumavir2"].append(vals[7])
    return {k: np.array(v) for k, v in out.items()}


def read_tr_dvrho(path):
    lines = open(path).read().splitlines()
    rho = []
    dvrho = []
    for line in lines[2:]:
        a, b = [float(x) for x in line.split(',')]
        rho.append(a)
        dvrho.append(b)
    return np.array(rho), np.array(dvrho)


def integrate_psit(psip, q):
    psit = np.zeros_like(psip)
    for i in range(1, len(psip)):
        qeff = 2.0 * q[i] * q[i - 1] / (q[i] + q[i - 1])
        psit[i] = psit[i - 1] + qeff * (psip[i] - psip[i - 1])
    return psit


def read_gfile_q(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    hdr = lines[0].split()
    nw = int(hdr[-2]); nh = int(hdr[-1])
    def parse(line):
        res=[]
        for i in range(0,len(line),16):
            c=line[i:i+16].strip()
            if not c:
                continue
            if '-' in c[1:] and 'E' not in c and 'e' not in c:
                idx=c.rfind('-'); c=c[:idx]+'E'+c[idx:]
            elif '+' in c[1:] and 'E' not in c and 'e' not in c:
                idx=c.rfind('+'); c=c[:idx]+'E'+c[idx:]
            try:
                res.append(float(c))
            except ValueError:
                pass
        return res
    data=[]
    for line in lines[5:]:
        data.extend(parse(line.rstrip('\n')))
    qpsi = np.abs(np.array(data[4*nw + nw*nh : 4*nw + nw*nh + nw]))
    psi_n = np.linspace(0.0, 1.0, nw)
    return psi_n, qpsi


def read_gfile_psipa(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    def parse(line):
        res=[]
        for i in range(0,len(line),16):
            c=line[i:i+16].strip()
            if not c:
                continue
            if '-' in c[1:] and 'E' not in c and 'e' not in c:
                idx=c.rfind('-'); c=c[:idx]+'E'+c[idx:]
            elif '+' in c[1:] and 'E' not in c and 'e' not in c:
                idx=c.rfind('+'); c=c[:idx]+'E'+c[idx:]
            try:
                res.append(float(c))
            except ValueError:
                pass
        return res
    l3 = parse(lines[2])
    simag, sibry = l3[2], l3[3]
    return 2.0 * np.pi * (sibry - simag)


def compute_sumavir2_and_dvdpsip_from_gfile():
    R, Z, psi_n_2d, psi_p_2d, rmaxis, zmaxis = read_gfile_2d(GFILE)
    psi_n_map, rho_map = read_input_profiles_mapping(INPUT_PROFILES)
    dpsi_dZ, dpsi_dR = np.gradient(psi_p_2d, Z, R)
    fig_dummy, ax = plt.subplots()
    sumavir2 = np.zeros_like(rho_map)
    dvdpsip = np.zeros_like(rho_map)
    from matplotlib.path import Path
    for i, p_target in enumerate(psi_n_map):
        if p_target <= 1e-6:
            sumavir2[i] = 0.0
            dvdpsip[i] = 0.0
            continue
        cs = ax.contour(R, Z, psi_n_2d, levels=[p_target])
        paths=[]
        if len(cs.allsegs)>0 and len(cs.allsegs[0])>0:
            for seg in cs.allsegs[0]:
                paths.append(Path(seg))
        if not paths:
            sumavir2[i] = np.nan; dvdpsip[i] = np.nan; continue
        valid=[]
        for p in paths:
            v=p.vertices
            if np.min(v[:,0]) < rmaxis < np.max(v[:,0]) and np.min(v[:,1]) < zmaxis < np.max(v[:,1]):
                valid.append(p)
        if not valid:
            valid=paths
        path=max(valid,key=lambda x: len(x.vertices))
        v=path.vertices
        if not np.allclose(v[0],v[-1]):
            v=np.vstack([v,v[0]])
        s1=0.0; s2=0.0
        for k in range(1,len(v)):
            x1,y1=v[k-1]; x2,y2=v[k]
            dl=np.hypot(x2-x1,y2-y1)
            rm=0.5*(x1+x2); zm=0.5*(y1+y2)
            dpr=bilinear_interp(R,Z,dpsi_dR,rm,zm)
            dpz=bilinear_interp(R,Z,dpsi_dZ,rm,zm)
            bp=np.sqrt(dpr*dpr+dpz*dpz)/(2*np.pi*rm)
            if bp <= 1e-30 or rm <= 1e-30:
                continue
            s1 += dl/(bp*rm*rm)
            s2 += dl/bp
        sumavir2[i]=s1
        dvdpsip[i]=s2
    plt.close(fig_dummy)
    return rho_map, sumavir2, dvdpsip


def main():
    tr = read_tr_eqcalq(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))
    rho_dvrho_tr, dvrho_tr = read_tr_dvrho(os.path.join(BASE, 'tr_data_133.csv'))

    # gfile/OMFIT references
    rho_prof=[]; polflux_prof=[]
    with open(INPUT_PROFILES) as f:
        for line in f:
            if line.startswith('#') or 'N_ION' in line or '=' in line or not line.strip():
                continue
            parts=line.split()
            if len(parts)==5:
                try:
                    rho_prof.append(float(parts[0]))
                    polflux_prof.append(float(parts[2]))
                except ValueError:
                    pass
            if len(rho_prof)==201:
                break
    rho_prof=np.array(rho_prof)
    polflux_prof=np.array(polflux_prof)
    psi_n_q, q_g = read_gfile_q(GFILE)
    psipa_gf = read_gfile_psipa(GFILE)
    # Use input.profiles psi-grid shape, but gfile absolute PSIPA scale
    psi_n_prof = np.abs(polflux_prof) / np.abs(polflux_prof[-1])
    x_om = psi_n_prof
    psip_om = x_om * psipa_gf

    q_om = np.interp(x_om, psi_n_q, q_g)
    psit_om = integrate_psit(psip_om, q_om)

    extra = pd.read_csv(os.path.join(OMFIT_OUT, 'outputs_new', 'input_profiles_extra.csv'))
    vol = extra['EXPRO_vol'].values
    dvdpsip_om = np.abs(np.gradient(vol, psip_om))
    dvdpsit_om = dvdpsip_om / q_om

    psi_n_f, fpol = read_gfile_fpol(GFILE)
    tts_om = np.interp(x_om, psi_n_f, -2*np.pi*fpol)
    sumavir2_om = 4*np.pi**2 * q_om / np.maximum(tts_om, 1e-30)

    # gfile contour direct geometry
    rho_map, sumavir2_g, dvdpsip_g = compute_sumavir2_and_dvdpsip_from_gfile()
    x_g = np.interp(rho_map, rho_prof, x_om)

    x_tr = tr['psip'] / tr['psip'][-1]
    refs = {
        'psip': np.interp(x_tr, x_om, psip_om),
        'psit': np.interp(x_tr, x_om, psit_om),
        'qps': np.interp(x_tr, x_om, q_om),
        'tts': np.interp(x_tr, x_om, tts_om),
        'sumavir2': np.interp(x_tr, x_g, sumavir2_g),
        'dvdpsip': np.interp(x_tr, x_g, dvdpsip_g),
        'dvdpsit': np.interp(x_tr, x_om, dvdpsit_om),
    }

    # DVRHO references
    psita_tr = tr['psit'][-1]
    dvrho_from_eq_tr = tr['dvdpsit'] * 2.0 * psita_tr * tr['rho']
    refs['dvrho_eqchain_tr'] = dvrho_from_eq_tr
    if os.path.exists(os.path.join(BASE, 'gfile_dvrho.csv')):
        gfd = pd.read_csv(os.path.join(BASE, 'gfile_dvrho.csv'))
        refs['dvrho_gfile'] = np.interp(x_tr, np.interp(gfd['rho'].values, rho_prof, x_om), gfd['DVRHO_gfile'].values)
    else:
        refs['dvrho_gfile'] = np.interp(x_tr, x_om, np.gradient(extra['EXPRO_vol'].values, rho_prof))

    items = [
        ('psip', 'PSIP'), ('psit', 'PSIT'), ('qps', 'QPS / q'), ('tts', 'TTS'),
        ('sumavir2', 'SUMAVIR2'), ('dvdpsip', 'dV/dpsi_p'), ('dvdpsit', 'dV/dpsi_t'), ('dvrho', 'DVRHO')
    ]
    fig, axes = plt.subplots(len(items), 2, figsize=(13, 22), sharex='col')
    for i,(key,title) in enumerate(items):
        if key == 'dvrho':
            x_dv = np.interp(rho_dvrho_tr, tr['rho'], x_tr)
            axes[i,0].plot(x_dv, dvrho_tr, 'o-', ms=2.8, lw=1.6, label='TR DVRHO from tr_data_133')
            axes[i,0].plot(x_tr, refs['dvrho_eqchain_tr'], '--', lw=1.8, color='tab:green', label='TR eq-chain DVRHO')
            axes[i,0].plot(x_tr, refs['dvrho_gfile'], ':', lw=2.0, color='tab:red', label='gfile DVRHO')
            ratio = np.interp(x_tr, x_dv, dvrho_tr) / np.where(np.abs(refs['dvrho_gfile']) > 1e-30, refs['dvrho_gfile'], np.nan)
        else:
            axes[i,0].plot(x_tr, tr[key], 'o-', ms=2.8, lw=1.6, label=f'TR {title}')
            axes[i,0].plot(x_tr, refs[key], '--', lw=2.0, color='tab:red', label=f'gfile/OMFIT {title}')
            ratio = tr[key] / np.where(np.abs(refs[key]) > 1e-30, refs[key], np.nan)
        axes[i,0].set_ylabel(title)
        axes[i,0].set_title(f'{title} on normalized psi_p grid')
        axes[i,0].grid(True, alpha=0.3)
        axes[i,0].legend(fontsize=8)
        axes[i,1].plot(x_tr, ratio, 'd-', color='tab:blue', lw=1.4, ms=2.6)
        axes[i,1].axhline(1.0, color='k', lw=1, alpha=0.5)
        axes[i,1].set_ylabel('TR / ref')
        axes[i,1].set_title(f'{title} ratio')
        axes[i,1].grid(True, alpha=0.3)

    axes[-1,0].set_xlabel('psi_p / psi_pa')
    axes[-1,1].set_xlabel('psi_p / psi_pa')
    out = os.path.join(BASE, 'equilibrium_chain_big_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=170)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()

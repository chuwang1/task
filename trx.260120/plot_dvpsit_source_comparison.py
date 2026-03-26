#!/usr/bin/env python3
"""
Compare the source quantities of dvpsit between TR EQ module and gfile/OMFIT.

dvpsit = dV/dpsi_t = (dV/dpsi_p) / q

This script compares:
  1. dV/dpsi_p : TR vs OMFIT
  2. q         : TR vs OMFIT/gfile
  3. dV/dpsi_t : TR vs OMFIT (derived)
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new"
GFILE = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114"


def read_eqcalq_raw(path):
    with open(path) as f:
        lines = f.readlines()
    data = {"nr": [], "rhot": [], "psip": [], "dvdpsit": [], "dvdpsip": [], "qps": []}
    for line in lines[1:]:
        parts = line.strip().split(",")
        nr = int(parts[0])
        vals = [float(x) for x in parts[1].split()]
        data["nr"].append(nr)
        data["rhot"].append(vals[0])
        data["psip"].append(vals[1])
        data["dvdpsit"].append(vals[2])
        data["dvdpsip"].append(vals[3])
        data["qps"].append(vals[4])
    return {k: np.array(v) for k, v in data.items()}


def read_gfile_qpsi(path):
    """Read q profile from gfile."""
    with open(path) as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    nw = int(header[-2])

    def parse_floats(line):
        res = []
        for i in range(0, len(line), 16):
            chunk = line[i:i + 16].strip()
            if not chunk:
                continue
            if '-' in chunk[1:] and 'E' not in chunk and 'e' not in chunk:
                idx = chunk.rfind('-')
                chunk = chunk[:idx] + 'E' + chunk[idx:]
            elif '+' in chunk[1:] and 'E' not in chunk and 'e' not in chunk:
                idx = chunk.rfind('+')
                chunk = chunk[:idx] + 'E' + chunk[idx:]
            try:
                res.append(float(chunk))
            except ValueError:
                pass
        return res

    l2 = parse_floats(lines[1])
    l3 = parse_floats(lines[2])
    simag, sibry = l3[2], l3[3]

    data = []
    for line in lines[5:]:
        data.extend(parse_floats(line.rstrip('\n')))

    # Skip fpol, pres, ffprime, pprime, psirz, qpsi
    # fpol: nw, pres: nw, ffprime: nw, pprime: nw
    # psirz: nw*nh (already skipped by offset 4*nw + nw*nh)
    nh = int(header[-1])
    cur = 4 * nw + nw * nh
    qpsi = np.array(data[cur: cur + nw])

    psi_n = np.linspace(0, 1, nw)
    psi_p = simag + psi_n * (sibry - simag)
    return psi_n, psi_p, qpsi, simag, sibry


def main():
    import pandas as pd

    # --- TR data ---
    tr = read_eqcalq_raw(os.path.join(BASE, "eqcalq_dvdpsit_raw.csv"))

    # --- OMFIT data ---
    allp = pd.read_csv(os.path.join(OMFIT, "all_profiles.csv"))
    extra = pd.read_csv(os.path.join(OMFIT, "input_profiles_extra.csv"))

    rho_om = allp["rho"].values
    vol_om = extra["EXPRO_vol"].values

    # polflux from input.profiles
    inp_path = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/input.profiles"
    rho_inp = []
    polflux_inp = []
    with open(inp_path) as f:
        for line in f:
            if line.startswith('#') or 'N_ION' in line or '=' in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) == 5:
                try:
                    rho_inp.append(float(parts[0]))
                    polflux_inp.append(float(parts[2]))
                except ValueError:
                    pass
            if len(rho_inp) == 201:
                break
    rho_inp = np.array(rho_inp)
    polflux_inp = np.array(polflux_inp)  # psi_p in physical units

    # dV/dpsi_p from OMFIT
    # TR uses PSIP = 2*pi * |polflux|, so convert OMFIT polflux to TR convention
    psip_tr_convention = np.abs(polflux_inp) * 2.0 * np.pi
    dvdpsip_om = np.abs(np.gradient(vol_om, psip_tr_convention))

    # q from gfile
    psi_n_gf, psi_p_gf, q_gf, simag, sibry = read_gfile_qpsi(GFILE)
    # interpolate gfile q onto OMFIT rho grid
    # need psi_n for OMFIT points
    psi_n_om = (polflux_inp - polflux_inp[0]) / (polflux_inp[-1] - polflux_inp[0])
    q_om = np.interp(psi_n_om, psi_n_gf, q_gf)

    # dV/dpsi_t from OMFIT
    dvdpsit_om = dvdpsip_om / q_om

    # --- Interpolate OMFIT onto TR rho grid ---
    rho_tr = tr["rhot"]
    dvdpsip_om_tr = np.interp(rho_tr, rho_om, dvdpsip_om)
    q_om_tr = np.interp(rho_tr, rho_om, q_om)
    dvdpsit_om_tr = np.interp(rho_tr, rho_om, dvdpsit_om)

    # --- Plot ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    # Row 1: absolute values
    axes[0, 0].plot(rho_tr, tr["dvdpsip"], "o-", ms=3, lw=2, label="TR dV/dpsi_p")
    axes[0, 0].plot(rho_tr, dvdpsip_om_tr, "--", lw=2, color="tab:red", label="OMFIT dV/dpsi_p")
    axes[0, 0].set_ylabel("dV/dpsi_p")
    axes[0, 0].set_title("dV/dpsi_p")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(rho_tr, tr["qps"], "o-", ms=3, lw=2, label="TR q")
    axes[0, 1].plot(rho_tr, q_om_tr, "--", lw=2, color="tab:red", label="gfile q")
    axes[0, 1].set_ylabel("q")
    axes[0, 1].set_title("Safety Factor q")
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    axes[0, 2].plot(rho_tr, tr["dvdpsit"], "o-", ms=3, lw=2, label="TR dV/dpsi_t")
    axes[0, 2].plot(rho_tr, dvdpsit_om_tr, "--", lw=2, color="tab:red", label="OMFIT dV/dpsi_t")
    axes[0, 2].set_ylabel("dV/dpsi_t")
    axes[0, 2].set_title("dV/dpsi_t = (dV/dpsi_p)/q")
    axes[0, 2].grid(True, alpha=0.3)
    axes[0, 2].legend(fontsize=8)

    # Row 2: ratios TR/OMFIT
    eps = 1e-30
    r1 = tr["dvdpsip"] / (dvdpsip_om_tr + eps)
    r2 = tr["qps"] / (q_om_tr + eps)
    r3 = tr["dvdpsit"] / (dvdpsit_om_tr + eps)

    axes[1, 0].plot(rho_tr, r1, "d-", ms=3, lw=1.5, color="tab:blue")
    axes[1, 0].axhline(1.0, color="k", lw=1, alpha=0.5)
    axes[1, 0].set_xlabel("rho")
    axes[1, 0].set_ylabel("TR / OMFIT")
    axes[1, 0].set_title("dV/dpsi_p ratio")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(rho_tr, r2, "d-", ms=3, lw=1.5, color="tab:blue")
    axes[1, 1].axhline(1.0, color="k", lw=1, alpha=0.5)
    axes[1, 1].set_xlabel("rho")
    axes[1, 1].set_ylabel("TR / gfile")
    axes[1, 1].set_title("q ratio")
    axes[1, 1].grid(True, alpha=0.3)

    axes[1, 2].plot(rho_tr, r3, "d-", ms=3, lw=1.5, color="tab:blue")
    axes[1, 2].axhline(1.0, color="k", lw=1, alpha=0.5)
    axes[1, 2].set_xlabel("rho")
    axes[1, 2].set_ylabel("TR / OMFIT")
    axes[1, 2].set_title("dV/dpsi_t ratio")
    axes[1, 2].grid(True, alpha=0.3)

    out = os.path.join(BASE, "dvpsit_source_comparison.png")
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f"Wrote {out}")

    # Print edge comparison
    print("\n=== Edge comparison ===")
    print(f"{'rho':>6s} {'TR_dvdpsip':>12s} {'OM_dvdpsip':>12s} {'TR_q':>8s} {'OM_q':>8s} {'TR_dvdpsit':>12s} {'OM_dvdpsit':>12s}")
    for i in range(len(rho_tr) - 8, len(rho_tr)):
        print(f"{rho_tr[i]:6.3f} {tr['dvdpsip'][i]:12.4f} {dvdpsip_om_tr[i]:12.4f} "
              f"{tr['qps'][i]:8.4f} {q_om_tr[i]:8.4f} "
              f"{tr['dvdpsit'][i]:12.4f} {dvdpsit_om_tr[i]:12.4f}")


if __name__ == "__main__":
    main()

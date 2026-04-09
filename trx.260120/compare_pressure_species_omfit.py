#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
from collections import OrderedDict

import f90nml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from solve_gs_from_gfile import read_gfile, rho_from_qpsi

RKEV = 1.602176634e-16
DEFAULT_GFILE = Path(__file__).resolve().parent / "g260206.20000_teq_0114"
DEFAULT_NAMELIST = Path("/Users/dengxiaoya/CFEDRSW/OMFIT_out/profiles_CFEDR.namelist")
DEFAULT_TR_DENSITY_CSV = Path(__file__).resolve().parent / "tr_data_074.csv"
DEFAULT_TR_TEMPERATURE_CSV = Path(__file__).resolve().parent / "tr_data_076.csv"
DEFAULT_TR_FAST_DENSITY_CSV = Path(__file__).resolve().parent / "tr_data_075.csv"
DEFAULT_TR_FAST_TEMPERATURE_CSV = Path(__file__).resolve().parent / "tr_data_030.csv"


def species_pressure_pa(block: dict) -> np.ndarray:
    density = np.asarray(block["density"], dtype=float) / 1e20
    temperature = np.asarray(block["temperature"], dtype=float)
    return density * temperature * 1e20 * RKEV


def species_density(block: dict) -> np.ndarray:
    return np.asarray(block["density"], dtype=float) / 1e20


def species_temperature(block: dict) -> np.ndarray:
    return np.asarray(block["temperature"], dtype=float)


def load_tr_pressure_csv(base_dir: Path) -> tuple[np.ndarray, np.ndarray, Path] | None:
    csv_paths = sorted(base_dir.glob("tr_data_*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
    for path in csv_paths:
        with path.open("r", encoding="utf-8") as fh:
            title = fh.readline().strip()
        if "@P(ns) [MPa]  vs r@" not in title:
            continue

        df = pd.read_csv(path, skiprows=1)
        df.columns = [str(c).strip() for c in df.columns]
        species_cols = [c for c in df.columns if c != "X"]
        if not species_cols:
            continue

        rho = df["X"].to_numpy(dtype=float)
        p_total_pa = df[species_cols].sum(axis=1).to_numpy(dtype=float) * 1e6
        return rho, p_total_pa, path
    return None


def load_tr_pressure_from_density_temperature(
    density_csv: Path, temperature_csv: Path
) -> tuple[np.ndarray, np.ndarray, str]:
    dn = pd.read_csv(density_csv, skiprows=1)
    tt = pd.read_csv(temperature_csv, skiprows=1)
    dn.columns = [str(c).strip() for c in dn.columns]
    tt.columns = [str(c).strip() for c in tt.columns]

    rho_n = dn["X"].to_numpy(dtype=float)
    rho_t = tt["X"].to_numpy(dtype=float)
    if rho_n.shape != rho_t.shape or not np.allclose(rho_n, rho_t, rtol=0.0, atol=1e-10):
        raise ValueError("TR density/temperature CSV grids do not match")

    pressure_cols = [("nE", "TE"), ("nD", "TD"), ("nT", "TT"), ("nA", "TA")]
    p_total_pa = np.zeros_like(rho_n, dtype=float)
    for ncol, tcol in pressure_cols:
        p_total_pa += dn[ncol].to_numpy(dtype=float) * tt[tcol].to_numpy(dtype=float) * 1e20 * RKEV
    label = f"TR nT ({density_csv.name}+{temperature_csv.name})"
    return rho_n, p_total_pa, label


def load_tr_species_profiles(density_csv: Path, temperature_csv: Path) -> tuple[np.ndarray, dict[str, dict[str, np.ndarray]]]:
    dn = pd.read_csv(density_csv, skiprows=1)
    tt = pd.read_csv(temperature_csv, skiprows=1)
    dn.columns = [str(c).strip() for c in dn.columns]
    tt.columns = [str(c).strip() for c in tt.columns]

    rho_n = dn["X"].to_numpy(dtype=float)
    rho_t = tt["X"].to_numpy(dtype=float)
    if rho_n.shape != rho_t.shape or not np.allclose(rho_n, rho_t, rtol=0.0, atol=1e-10):
        raise ValueError("TR density/temperature CSV grids do not match")

    tr_species = OrderedDict()
    for name, ncol, tcol in [
        ("electron", "nE", "TE"),
        ("d", "nD", "TD"),
        ("t", "nT", "TT"),
        ("he4", "nA", "TA"),
    ]:
        density = dn[ncol].to_numpy(dtype=float)
        temperature = tt[tcol].to_numpy(dtype=float)
        pressure = density * temperature * 1e20 * RKEV
        tr_species[name] = {
            "density": density,
            "temperature": temperature,
            "pressure": pressure,
        }
    return rho_n, tr_species


def load_tr_fast_alpha_profiles(
    density_csv: Path, temperature_csv: Path
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    dn = pd.read_csv(density_csv, skiprows=1)
    tt = pd.read_csv(temperature_csv, skiprows=1)
    dn.columns = [str(c).strip() for c in dn.columns]
    tt.columns = [str(c).strip() for c in tt.columns]

    rho_n = dn["X"].to_numpy(dtype=float)
    rho_t = tt["X"].to_numpy(dtype=float)
    if rho_n.shape != rho_t.shape or not np.allclose(rho_n, rho_t, rtol=0.0, atol=1e-10):
        raise ValueError("TR fast-alpha density/temperature CSV grids do not match")

    density = dn["NF"].to_numpy(dtype=float)
    temperature = tt["TF"].to_numpy(dtype=float)
    pressure = density * temperature * 1e20 * RKEV
    return rho_n, {
        "density": density,
        "temperature": temperature,
        "pressure": pressure,
    }


def aggregate_omfit_species(nml: dict) -> OrderedDict[str, dict[str, np.ndarray]]:
    species = OrderedDict()
    he4_count = 0

    def add_species(name: str, density: np.ndarray, temperature: np.ndarray, pressure: np.ndarray) -> None:
        if name in species:
            prev = species[name]
            density_sum = prev["density"] + density
            pressure_sum = prev["pressure"] + pressure
            temperature_eff = np.divide(
                pressure_sum,
                density_sum * 1e20 * RKEV,
                out=np.zeros_like(density_sum),
                where=density_sum > 1e-30,
            )
            species[name] = {
                "density": density_sum,
                "temperature": temperature_eff,
                "pressure": pressure_sum,
            }
            return

        species[name] = {
            "density": density.copy(),
            "temperature": temperature.copy(),
            "pressure": pressure.copy(),
        }

    e_block = nml["electron"]
    add_species("electron", species_density(e_block), species_temperature(e_block), species_pressure_pa(e_block))

    for key, block in nml.items():
        if not (isinstance(block, dict) and str(key).lower().startswith("ions_")):
            continue
        name = str(block.get("name", key)).strip().lower()
        if name == "he4":
            he4_count += 1
            if he4_count >= 2:
                name = "he4_fast"
        add_species(name, species_density(block), species_temperature(block), species_pressure_pa(block))

    return species


def main():
    base_dir = Path(__file__).resolve().parent
    g = read_gfile(str(DEFAULT_GFILE))
    nml = f90nml.read(str(DEFAULT_NAMELIST))
    tr_pressure = load_tr_pressure_csv(base_dir)
    tr_pressure_nt = load_tr_pressure_from_density_temperature(
        DEFAULT_TR_DENSITY_CSV,
        DEFAULT_TR_TEMPERATURE_CSV,
    )
    rho_tr_species, tr_species = load_tr_species_profiles(
        DEFAULT_TR_DENSITY_CSV,
        DEFAULT_TR_TEMPERATURE_CSV,
    )
    rho_tr_fast, tr_fast_alpha = load_tr_fast_alpha_profiles(
        DEFAULT_TR_FAST_DENSITY_CSV,
        DEFAULT_TR_FAST_TEMPERATURE_CSV,
    )
    if rho_tr_species.shape != rho_tr_fast.shape or not np.allclose(rho_tr_species, rho_tr_fast, rtol=0.0, atol=1e-10):
        raise ValueError("TR thermal and fast-alpha grids do not match")
    tr_species["he4_fast"] = tr_fast_alpha

    rho = np.asarray(nml["powers_particle_flux_onetwo"]["rho"], dtype=float)
    rho_g = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)
    psi_n_g = np.asarray(g.psi_n, dtype=float)
    p_g = np.asarray(g.pres, dtype=float)
    p_g_on_rho = np.interp(rho, rho_g, p_g)
    psi_of_rho = interp1d(
        rho_g[np.argsort(rho_g)],
        psi_n_g[np.argsort(rho_g)],
        bounds_error=False,
        fill_value="extrapolate",
    )
    psi_n = psi_of_rho(rho)

    species = aggregate_omfit_species(nml)
    p_total = np.zeros_like(rho)
    for vals in species.values():
        p_total += vals["pressure"]

    fig, axes = plt.subplots(3, 2, figsize=(14, 15), constrained_layout=True)
    species_colors: dict[str, str] = {}

    for name, vals in species.items():
        line = axes[0, 0].plot(rho, vals["pressure"] / 1e6, lw=1.8, label=name)[0]
        species_colors[name] = line.get_color()
        if name in tr_species:
            axes[0, 0].plot(
                rho_tr_species,
                tr_species[name]["pressure"] / 1e6,
                lw=1.8,
                ls="--",
                color=species_colors[name],
                label=f"TR {name}",
            )
    axes[0, 0].plot(rho, p_total / 1e6, "k-", lw=2.6, label="total nT")
    axes[0, 0].plot(rho, p_g_on_rho / 1e6, "k--", lw=2.2, label="gfile PRES")
    if tr_pressure is not None:
        rho_tr, p_tr, tr_path = tr_pressure
        axes[0, 0].plot(rho_tr, p_tr / 1e6, color="tab:purple", lw=2.2, ls="-.", label=f"TR total P ({tr_path.name})")
    rho_tr_nt, p_tr_nt, tr_nt_label = tr_pressure_nt
    p_tr_total = np.zeros_like(rho_tr_species)
    for vals in tr_species.values():
        p_tr_total += vals["pressure"]
    axes[0, 0].plot(rho_tr_nt, p_tr_nt / 1e6, color="tab:green", lw=2.0, ls=":", label=tr_nt_label)
    axes[0, 0].plot(rho_tr_species, p_tr_total / 1e6, color="tab:brown", lw=2.0, ls=":", label="TR total incl. fast alpha")
    axes[0, 0].set_title("Pressure Components vs rho_tor")
    axes[0, 0].set_xlabel("rho_tor")
    axes[0, 0].set_ylabel("Pressure (MPa)")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=3)

    for name, vals in species.items():
        axes[0, 1].plot(psi_n, vals["pressure"] / 1e6, lw=1.8, color=species_colors[name], label=name)
        if name in tr_species:
            psi_tr_species = psi_of_rho(rho_tr_species)
            axes[0, 1].plot(
                psi_tr_species,
                tr_species[name]["pressure"] / 1e6,
                lw=1.8,
                ls="--",
                color=species_colors[name],
                label=f"TR {name}",
            )
    axes[0, 1].plot(psi_n, p_total / 1e6, "k-", lw=2.6, label="total nT")
    axes[0, 1].plot(psi_n, np.interp(psi_n, psi_n_g, p_g) / 1e6, "k--", lw=2.2, label="gfile PRES")
    if tr_pressure is not None:
        rho_tr, p_tr, tr_path = tr_pressure
        psi_tr = psi_of_rho(rho_tr)
        axes[0, 1].plot(psi_tr, p_tr / 1e6, color="tab:purple", lw=2.2, ls="-.", label=f"TR total P ({tr_path.name})")
    psi_tr_nt = psi_of_rho(rho_tr_nt)
    axes[0, 1].plot(psi_tr_nt, p_tr_nt / 1e6, color="tab:green", lw=2.0, ls=":", label=tr_nt_label)
    psi_tr_total = psi_of_rho(rho_tr_species)
    axes[0, 1].plot(psi_tr_total, p_tr_total / 1e6, color="tab:brown", lw=2.0, ls=":", label="TR total incl. fast alpha")
    axes[0, 1].set_title("Pressure Components vs psi_n")
    axes[0, 1].set_xlabel("psi_n")
    axes[0, 1].set_ylabel("Pressure (MPa)")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend(fontsize=8, ncol=3)

    frac_axis = []
    labels = []
    p_axis_total = float(p_total[0])
    for name, vals in species.items():
        labels.append(name)
        frac_axis.append(100.0 * float(vals["pressure"][0]) / max(p_axis_total, 1e-30))
    axes[1, 0].bar(labels, frac_axis, color="tab:blue")
    axes[1, 0].set_title("Species Fraction at Axis")
    axes[1, 0].set_ylabel("Percent of total pressure (%)")
    axes[1, 0].tick_params(axis="x", rotation=30)
    axes[1, 0].grid(alpha=0.25, axis="y")

    diff = p_total - p_g_on_rho
    axes[1, 1].plot(rho, diff / 1e3, color="tab:red", lw=2.0, label="OMFIT total nT - gfile")
    if tr_pressure is not None:
        rho_tr, p_tr, _tr_path = tr_pressure
        p_g_on_tr = np.interp(rho_tr, rho_g, p_g)
        axes[1, 1].plot(rho_tr, (p_tr - p_g_on_tr) / 1e3, color="tab:purple", lw=2.0, ls="-.", label="TR total P - gfile")
    p_g_on_tr_nt = np.interp(rho_tr_nt, rho_g, p_g)
    axes[1, 1].plot(rho_tr_nt, (p_tr_nt - p_g_on_tr_nt) / 1e3, color="tab:green", lw=2.0, ls=":", label="TR nT - gfile")
    p_g_on_tr_total = np.interp(rho_tr_species, rho_g, p_g)
    axes[1, 1].plot(rho_tr_species, (p_tr_total - p_g_on_tr_total) / 1e3, color="tab:brown", lw=2.0, ls=":", label="TR total incl. fast alpha - gfile")
    axes[1, 1].axhline(0.0, color="k", ls="--", lw=1.0)
    axes[1, 1].set_title("Pressure Difference to gfile PRES")
    axes[1, 1].set_xlabel("rho_tor")
    axes[1, 1].set_ylabel("Difference (kPa)")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=8)

    for name, vals in species.items():
        axes[2, 0].plot(rho, vals["density"], lw=1.8, color=species_colors[name], label=name)
        if name in tr_species:
            axes[2, 0].plot(
                rho_tr_species,
                tr_species[name]["density"],
                lw=1.8,
                ls="--",
                color=species_colors[name],
                label=f"TR {name}",
            )
    axes[2, 0].set_title("Species Density vs rho_tor")
    axes[2, 0].set_xlabel("rho_tor")
    axes[2, 0].set_ylabel(r"Density ($10^{20}$ m$^{-3}$)")
    axes[2, 0].grid(alpha=0.25)
    axes[2, 0].legend(fontsize=8, ncol=3)

    for name, vals in species.items():
        axes[2, 1].plot(rho, vals["temperature"], lw=1.8, color=species_colors[name], label=name)
        if name in tr_species:
            axes[2, 1].plot(
                rho_tr_species,
                tr_species[name]["temperature"],
                lw=1.8,
                ls="--",
                color=species_colors[name],
                label=f"TR {name}",
            )
    axes[2, 1].set_title("Species Temperature vs rho_tor")
    axes[2, 1].set_xlabel("rho_tor")
    axes[2, 1].set_ylabel("Temperature (keV)")
    axes[2, 1].grid(alpha=0.25)
    axes[2, 1].legend(fontsize=8, ncol=3)

    fig.suptitle("OMFIT Species Pressure vs gfile Pressure")
    out = Path(__file__).resolve().parent.parent / "omfit_species_pressure_breakdown.png"
    fig.savefig(out, dpi=180)
    print(out)


if __name__ == "__main__":
    main()

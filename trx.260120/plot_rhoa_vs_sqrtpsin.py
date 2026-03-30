#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_sumavir2_from_gfile import read_input_profiles_mapping


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_OUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
INPUT_PROFILES = os.path.join(OMFIT_OUT, "input.profiles")


def main():
    psi_n, rho_a = read_input_profiles_mapping(INPUT_PROFILES)
    sqrt_psin = np.sqrt(np.clip(psi_n, 0.0, None))

    diff = rho_a - sqrt_psin
    abs_diff = np.abs(diff)

    drho_dra = np.gradient(sqrt_psin, rho_a)
    dra_drho = np.gradient(rho_a, sqrt_psin, edge_order=1)

    out_df = pd.DataFrame({
        "r_over_a": rho_a,
        "psi_n": psi_n,
        "sqrt_psi_n": sqrt_psin,
        "r_over_a_minus_sqrt_psi_n": diff,
        "abs_diff": abs_diff,
        "d_sqrtpsin_d_rhoa": drho_dra,
        "d_rhoa_d_sqrtpsin": dra_drho,
    })
    csv_out = os.path.join(BASE, "rhoa_vs_sqrtpsin.csv")
    out_df.to_csv(csv_out, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    ax.plot(rho_a, rho_a, "k--", lw=1.2, label=r"$r/a$")
    ax.plot(rho_a, sqrt_psin, "r-", lw=2.0, label=r"$\sqrt{\psi_n}$")
    ax.set_xlabel(r"$r/a$")
    ax.set_ylabel("value")
    ax.set_title(r"$r/a$ and $\sqrt{\psi_n}$")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(rho_a, diff, "b-", lw=2.0)
    ax.axhline(0.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel(r"$r/a$")
    ax.set_ylabel(r"$(r/a) - \sqrt{\psi_n}$")
    ax.set_title("Difference")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(rho_a, drho_dra, "m-", lw=2.0)
    ax.axhline(1.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel(r"$r/a$")
    ax.set_ylabel(r"$d\sqrt{\psi_n}/d(r/a)$")
    ax.set_title("Slope")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(sqrt_psin, rho_a, "g-", lw=2.0)
    ax.plot([0, 1], [0, 1], "k--", lw=1.0)
    ax.set_xlabel(r"$\sqrt{\psi_n}$")
    ax.set_ylabel(r"$r/a$")
    ax.set_title(r"Mapping: $\sqrt{\psi_n} \rightarrow r/a$")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_out = os.path.join(BASE, "rhoa_vs_sqrtpsin.png")
    plt.savefig(fig_out, dpi=170)
    plt.close(fig)

    idx_max = int(np.nanargmax(abs_diff))
    print(f"Wrote {fig_out}")
    print(f"Wrote {csv_out}")
    print(f"max |r/a - sqrt(psi_n)| = {abs_diff[idx_max]:.6f} at r/a = {rho_a[idx_max]:.3f}")
    print(f"mean |r/a - sqrt(psi_n)| = {np.nanmean(abs_diff):.6f}")
    print(f"max |d(sqrt(psi_n))/d(r/a) - 1| = {np.nanmax(np.abs(drho_dra - 1.0)):.6f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Compare psi_N definitions/mappings:
  - OMFIT namelist psi_normalized
  - rho^2 from namelist rho (rho = sqrt(normalized toroidal flux))
  - rho (for reference)
"""

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_top_level_array(filepath: str, varname: str):
    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    pattern = re.compile(rf"^\s*{re.escape(varname)}\s*=\s*", re.MULTILINE | re.IGNORECASE)
    m = pattern.search(text)
    if m is None:
        return None

    vals = []
    for line in text[m.end() :].splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("&") or re.match(r"^[A-Za-z_]\w*\s*=", s):
            break
        for token in s.replace(",", " ").split():
            if "*" in token:
                c, v = token.split("*", 1)
                try:
                    vals.extend([float(v)] * int(c))
                except Exception:
                    pass
            else:
                try:
                    vals.append(float(token))
                except Exception:
                    pass
    if not vals:
        return None
    return np.asarray(vals, dtype=float)


def err(a, b):
    d = a - b
    rmse = float(np.sqrt(np.mean(d**2)))
    mae = float(np.mean(np.abs(d)))
    mx = float(np.max(np.abs(d)))
    return rmse, mae, mx


def main():
    parser = argparse.ArgumentParser(description="Compare psi_N mappings from namelist.")
    parser.add_argument(
        "--namelist",
        default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/profiles_CFEDR.namelist",
    )
    parser.add_argument("--out-prefix", default="in/psiN_mapping_compare")
    args = parser.parse_args()

    if not os.path.exists(args.namelist):
        raise FileNotFoundError(args.namelist)

    rho = parse_top_level_array(args.namelist, "rho")
    psi_n = parse_top_level_array(args.namelist, "psi_normalized")
    if rho is None or psi_n is None:
        raise RuntimeError("rho or psi_normalized not found in namelist top-level.")

    n = min(len(rho), len(psi_n))
    rho = rho[:n]
    psi_n = psi_n[:n]

    psi_from_rho2 = rho**2
    psi_from_rho = rho

    e2 = err(psi_from_rho2, psi_n)
    e1 = err(psi_from_rho, psi_n)

    # Best shift for (rho+s)^2
    best_rmse = 1e99
    best_s = None
    for s in np.linspace(-0.2, 0.2, 801):
        ps = np.clip(rho + s, 0.0, 1.0) ** 2
        r, _, _ = err(ps, psi_n)
        if r < best_rmse:
            best_rmse = r
            best_s = float(s)

    print(f"n={n}")
    print(
        f"psi_n vs rho^2: rmse={e2[0]:.4e}, mae={e2[1]:.4e}, max={e2[2]:.4e}"
    )
    print(
        f"psi_n vs rho  : rmse={e1[0]:.4e}, mae={e1[1]:.4e}, max={e1[2]:.4e}"
    )
    print(f"best shift s for (rho+s)^2: s={best_s:.4f}, rmse={best_rmse:.4e}")

    out_prefix = os.path.abspath(args.out_prefix)
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    out_csv = f"{out_prefix}.csv"
    out_png = f"{out_prefix}.png"

    df = pd.DataFrame(
        {
            "rho": rho,
            "psi_n_namelist": psi_n,
            "rho_squared": psi_from_rho2,
            "rho_linear": psi_from_rho,
            "diff_rho2_minus_psi_n": psi_from_rho2 - psi_n,
            "diff_rho_minus_psi_n": psi_from_rho - psi_n,
        }
    )
    df.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    axes[0].plot(rho, psi_n, "k-", lw=2, label="psi_normalized (namelist)")
    axes[0].plot(rho, psi_from_rho2, "b--", lw=1.8, label="rho^2")
    axes[0].plot(rho, psi_from_rho, "r-.", lw=1.8, label="rho")
    axes[0].set_xlabel("rho")
    axes[0].set_ylabel("psi_N")
    axes[0].set_title("psi_N Mapping")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(rho, psi_from_rho2 - psi_n, "b--", lw=1.8, label="rho^2 - psi_n")
    axes[1].plot(rho, psi_from_rho - psi_n, "r-.", lw=1.8, label="rho - psi_n")
    axes[1].axhline(0.0, color="k", lw=1)
    axes[1].set_xlabel("rho")
    axes[1].set_ylabel("difference")
    axes[1].set_title("Difference")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    fig.suptitle("Namelist psi_N vs rho mappings")
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    print(f"saved: {out_csv}")
    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = "/Users/dengxiaoya/TASK/latest/task/trx"


def main():
    path = os.path.join(BASE, "eqbpsd_metric1d_debug.csv")
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    header = [h.strip() for h in lines[0].split(",")]
    cols = {h: [] for h in header}
    for line in lines[1:]:
        vals = [float(x) for x in line.split()]
        for h, v in zip(header, vals):
            cols[h].append(v)

    rho = cols["rho"]
    dvpsit = cols["dvpsit"]

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})

    axes[0].plot(rho, dvpsit, "o-", lw=2, ms=4, label="dvpsit = dV/dpsit")
    axes[0].set_ylabel("dV/dpsit")
    axes[0].set_title("eqbpsd metric1D dvpsit")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    rel = [v / dvpsit[0] for v in dvpsit]
    axes[1].plot(rho, rel, "d-", color="tab:red", lw=1.5, ms=3)
    axes[1].axhline(1.0, color="k", lw=1, alpha=0.5)
    axes[1].set_xlabel("rho")
    axes[1].set_ylabel("/axis")
    axes[1].grid(True, alpha=0.3)

    out = os.path.join(BASE, "eqbpsd_dvpsit.png")
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

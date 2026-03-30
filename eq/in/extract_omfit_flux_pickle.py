#!/usr/bin/env python3
import argparse
import csv
import sys

import numpy as np


def _to_array(x):
    try:
        return np.asarray(x, dtype=float)
    except Exception:
        try:
            return np.asarray(list(x), dtype=float)
        except Exception:
            return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickle", default="in/OMFITpickled_80f062a3fe")
    ap.add_argument("--omfit-src", default="/Users/dengxiaoya/CFEDRSW/OMFIT-source_unstable/omfit")
    ap.add_argument("--out", default="in/omfit_pickled_80f062a3fe_flux_export.csv")
    args = ap.parse_args()

    import scipy.interpolate.polyint as poly
    import scipy.interpolate._polyint as _p

    if not hasattr(poly, "_Interpolator1D") and hasattr(_p, "_Interpolator1D"):
        poly._Interpolator1D = _p._Interpolator1D

    if args.omfit_src not in sys.path:
        sys.path.insert(0, args.omfit_src)

    import pickle

    with open(args.pickle, "rb") as f:
        obj = pickle.load(f)

    fields = {
        "geo_psi": _to_array(obj["geo"]["psi"]),
        "geo_psin": _to_array(obj["geo"]["psin"]),
        "geo_rho": _to_array(obj["geo"]["rho"]),
        "geo_rhon": _to_array(obj["geo"]["rhon"]),
        "geo_R": _to_array(obj["geo"]["R"]),
        "geo_a": _to_array(obj["geo"]["a"]),
        "geo_Z": _to_array(obj["geo"]["Z"]),
        "geo_vol": _to_array(obj["geo"]["vol"]),
        "geo_bunit": _to_array(obj["geo"]["bunit"]),
        "avg_q": _to_array(obj["avg"]["q"]),
        "avg_PPRIME": _to_array(obj["avg"]["PPRIME"]),
        "avg_F": _to_array(obj["avg"]["F"]),
        "avg_1_over_R2": _to_array(obj["avg"]["1/R**2"]),
        "avg_vp": _to_array(obj["avg"]["vp"]),
        "avg_Bp2R2": _to_array(obj["avg"]["Bp**2*R**2"]),
    }

    n = min(len(v) for v in fields.values() if v is not None)
    for k in fields:
        fields[k] = fields[k][:n]

    r0 = float(obj["R0"])
    b0 = float(obj["BCENTR"])
    fields["derived_AVIR2_R0sq_times_avg1oR2"] = (r0 * r0) * fields["avg_1_over_R2"]

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["R0", "BCENTR"] + list(fields.keys())
        w.writerow(header)
        for i in range(n):
            row = [r0, b0] + [fields[k][i] for k in fields]
            w.writerow(row)

    print(f"loaded: {type(obj)}")
    print(f"R0={r0}, BCENTR={b0}, n={n}")
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()

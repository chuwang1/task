#!/usr/bin/env python3
import argparse
import csv
import math
import os


def load_grid(csv_path):
    rs = []
    zs = []
    values = {}
    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header or header[:3] != ["r", "z", "psi"]:
            raise ValueError("unexpected header in psi_rz.csv")
        for row in reader:
            r = float(row[0])
            z = float(row[1])
            psi = float(row[2])
            if r not in values:
                rs.append(r)
                values[r] = {}
            if z not in values[r]:
                values[r][z] = psi
            if z not in zs:
                zs.append(z)
    rs.sort()
    zs.sort()
    nr = len(rs)
    nz = len(zs)
    grid = [[0.0] * nr for _ in range(nz)]
    vmin = float("inf")
    vmax = float("-inf")
    for j, z in enumerate(zs):
        for i, r in enumerate(rs):
            psi = values[r][z]
            grid[j][i] = psi
            vmin = min(vmin, psi)
            vmax = max(vmax, psi)
    return rs, zs, grid, vmin, vmax


def lerp(a, b, t):
    return a + (b - a) * t


def colormap(t):
    t = max(0.0, min(1.0, t))
    if t < 0.25:
        # navy -> blue
        tt = t / 0.25
        return (
            int(lerp(0, 0, tt)),
            int(lerp(0, 64, tt)),
            int(lerp(64, 200, tt)),
        )
    if t < 0.5:
        # blue -> cyan
        tt = (t - 0.25) / 0.25
        return (
            int(lerp(0, 0, tt)),
            int(lerp(64, 200, tt)),
            int(lerp(200, 255, tt)),
        )
    if t < 0.75:
        # cyan -> yellow
        tt = (t - 0.5) / 0.25
        return (
            int(lerp(0, 255, tt)),
            int(lerp(200, 255, tt)),
            int(lerp(255, 0, tt)),
        )
    # yellow -> red
    tt = (t - 0.75) / 0.25
    return (
        int(lerp(255, 200, tt)),
        int(lerp(255, 0, tt)),
        int(lerp(0, 0, tt)),
    )


def write_ppm(path, grid, vmin, vmax, contour_levels):
    nz = len(grid)
    nr = len(grid[0])
    span = vmax - vmin if vmax > vmin else 1.0
    levels = [
        vmin + (span * i) / (contour_levels + 1)
        for i in range(1, contour_levels + 1)
    ]
    eps = span * 0.002
    with open(path, "wb") as fh:
        fh.write(f"P6\n{nr} {nz}\n255\n".encode("ascii"))
        for j in reversed(range(nz)):
            row = grid[j]
            for i in range(nr):
                val = row[i]
                t = (val - vmin) / span
                r, g, b = colormap(t)
                for level in levels:
                    if abs(val - level) <= eps:
                        r, g, b = 0, 0, 0
                        break
                fh.write(bytes((r, g, b)))


def main():
    parser = argparse.ArgumentParser(
        description="Render a psi(R,Z) contour heatmap from psi_rz.csv"
    )
    parser.add_argument("csv", help="Path to psi_rz.csv")
    parser.add_argument(
        "-o", "--output", default="psi_rz.ppm", help="Output PPM image path"
    )
    parser.add_argument(
        "--levels", type=int, default=15, help="Number of contour levels"
    )
    args = parser.parse_args()

    rs, zs, grid, vmin, vmax = load_grid(args.csv)
    _ = rs
    _ = zs
    write_ppm(args.output, grid, vmin, vmax, args.levels)
    print(
        f"Wrote {args.output} ({len(rs)} x {len(zs)}), "
        f"psi range [{vmin:.6g}, {vmax:.6g}]"
    )


if __name__ == "__main__":
    main()

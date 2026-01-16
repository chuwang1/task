#!/usr/bin/env python3
import argparse
import csv
import math


def read_ppm(path):
    with open(path, "rb") as fh:
        magic = fh.readline().strip()
        if magic != b"P6":
            raise ValueError(f"{path}: unsupported magic {magic!r}")
        def _next_token():
            while True:
                line = fh.readline()
                if not line:
                    raise ValueError("unexpected EOF in header")
                line = line.strip()
                if not line or line.startswith(b"#"):
                    continue
                return line
        dims = _next_token().split()
        if len(dims) == 2:
            width, height = map(int, dims)
        else:
            width = int(dims[0])
            height = int(_next_token())
        maxval = int(_next_token())
        if maxval != 255:
            raise ValueError(f"{path}: unsupported maxval {maxval}")
        data = fh.read(width * height * 3)
        if len(data) != width * height * 3:
            raise ValueError(f"{path}: truncated data")
        return width, height, bytearray(data)


def write_ppm(path, width, height, data):
    with open(path, "wb") as fh:
        fh.write(f"P6\n{width} {height}\n255\n".encode("ascii"))
        fh.write(data)


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


def find_plot_bbox(width, height, data, threshold=250):
    xmin, ymin = width, height
    xmax, ymax = -1, -1
    for y in range(height):
        row_offset = y * width * 3
        for x in range(width):
            idx = row_offset + x * 3
            r, g, b = data[idx], data[idx + 1], data[idx + 2]
            if r < threshold or g < threshold or b < threshold:
                if x < xmin:
                    xmin = x
                if x > xmax:
                    xmax = x
                if y < ymin:
                    ymin = y
                if y > ymax:
                    ymax = y
    if xmax < xmin or ymax < ymin:
        return 0, 0, width - 1, height - 1
    return xmin, ymin, xmax, ymax


def draw_line(data, width, height, x0, y0, x1, y1, color):
    x0 = int(round(x0))
    y0 = int(round(y0))
    x1 = int(round(x1))
    y1 = int(round(y1))
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    x, y = x0, y0
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    if dx >= dy:
        err = dx / 2
        while x != x1:
            if 0 <= x < width and 0 <= y < height:
                idx = (y * width + x) * 3
                data[idx:idx + 3] = color
            err -= dy
            if err < 0:
                y += sy
                err += dx
            x += sx
    else:
        err = dy / 2
        while y != y1:
            if 0 <= x < width and 0 <= y < height:
                idx = (y * width + x) * 3
                data[idx:idx + 3] = color
            err -= dx
            if err < 0:
                x += sx
                err += dy
            y += sy
    if 0 <= x1 < width and 0 <= y1 < height:
        idx = (y1 * width + x1) * 3
        data[idx:idx + 3] = color


def interpolate(p1, p2, v1, v2, level):
    if v2 == v1:
        return p1
    t = (level - v1) / (v2 - v1)
    return p1 + (p2 - p1) * t


def marching_squares(grid, levels):
    nz = len(grid)
    nr = len(grid[0])
    segments = {lvl: [] for lvl in levels}
    for j in range(nz - 1):
        for i in range(nr - 1):
            v00 = grid[j][i]
            v10 = grid[j][i + 1]
            v01 = grid[j + 1][i]
            v11 = grid[j + 1][i + 1]
            for lvl in levels:
                idx = 0
                if v00 >= lvl:
                    idx |= 1
                if v10 >= lvl:
                    idx |= 2
                if v11 >= lvl:
                    idx |= 4
                if v01 >= lvl:
                    idx |= 8
                if idx == 0 or idx == 15:
                    continue
                # Edges: 0=bottom,1=right,2=top,3=left
                edges = []
                if (idx & 1) != (idx & 2):
                    x = interpolate(i, i + 1, v00, v10, lvl)
                    edges.append((x, j))
                if (idx & 2) != (idx & 4):
                    y = interpolate(j, j + 1, v10, v11, lvl)
                    edges.append((i + 1, y))
                if (idx & 4) != (idx & 8):
                    x = interpolate(i + 1, i, v11, v01, lvl)
                    edges.append((x, j + 1))
                if (idx & 8) != (idx & 1):
                    y = interpolate(j + 1, j, v01, v00, lvl)
                    edges.append((i, y))
                if len(edges) == 2:
                    segments[lvl].append((edges[0], edges[1]))
                elif len(edges) == 4:
                    segments[lvl].append((edges[0], edges[1]))
                    segments[lvl].append((edges[2], edges[3]))
    return segments


def main():
    parser = argparse.ArgumentParser(
        description="Overlay contour lines from psi_rz.csv onto a PPM image."
    )
    parser.add_argument("csv", help="Path to psi_rz.csv")
    parser.add_argument("base_ppm", help="Base PPM image (equ plot)")
    parser.add_argument(
        "-o", "--output", default="overlay.ppm", help="Output PPM path"
    )
    parser.add_argument(
        "--levels", type=int, default=15, help="Number of contour levels"
    )
    parser.add_argument(
        "--swap-axes",
        action="store_true",
        help="Swap X/Y axes when mapping CSV grid to image",
    )
    parser.add_argument(
        "--no-flip-y",
        action="store_true",
        help="Do not flip Y axis (image origin at top-left by default)",
    )
    args = parser.parse_args()

    _, _, grid, vmin, vmax = load_grid(args.csv)
    width, height, data = read_ppm(args.base_ppm)
    xmin, ymin, xmax, ymax = find_plot_bbox(width, height, data)
    if xmax <= xmin or ymax <= ymin:
        xmin, ymin, xmax, ymax = 0, 0, width - 1, height - 1

    span = vmax - vmin if vmax > vmin else 1.0
    levels = [
        vmin + (span * i) / (args.levels + 1)
        for i in range(1, args.levels + 1)
    ]
    segments = marching_squares(grid, levels)

    nr = len(grid[0])
    nz = len(grid)
    plot_w = xmax - xmin
    plot_h = ymax - ymin

    color = bytes((0, 0, 0))
    def map_point(p):
        if args.swap_axes:
            x_norm = p[1] / (nz - 1)
            y_norm = p[0] / (nr - 1)
        else:
            x_norm = p[0] / (nr - 1)
            y_norm = p[1] / (nz - 1)
        if not args.no_flip_y:
            y_norm = 1.0 - y_norm
        x = xmin + x_norm * plot_w
        y = ymin + y_norm * plot_h
        return x, y

    for lvl in levels:
        for (p0, p1) in segments[lvl]:
            x0, y0 = map_point(p0)
            x1, y1 = map_point(p1)
            draw_line(data, width, height, x0, y0, x1, y1, color)

    write_ppm(args.output, width, height, data)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

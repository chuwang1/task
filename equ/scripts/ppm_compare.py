#!/usr/bin/env python3
import argparse


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
        return width, height, data


def write_ppm(path, width, height, data):
    with open(path, "wb") as fh:
        fh.write(f"P6\n{width} {height}\n255\n".encode("ascii"))
        fh.write(data)


def resize_nearest(width, height, data, new_width, new_height):
    if new_width == width and new_height == height:
        return data
    out = bytearray(new_width * new_height * 3)
    for y in range(new_height):
        src_y = int(y * height / new_height)
        for x in range(new_width):
            src_x = int(x * width / new_width)
            src_idx = (src_y * width + src_x) * 3
            dst_idx = (y * new_width + x) * 3
            out[dst_idx:dst_idx + 3] = data[src_idx:src_idx + 3]
    return bytes(out)


def main():
    parser = argparse.ArgumentParser(
        description="Compare two PPM images side-by-side."
    )
    parser.add_argument("left", help="Left PPM path")
    parser.add_argument("right", help="Right PPM path")
    parser.add_argument(
        "-o", "--output", default="ppm_compare.ppm", help="Output PPM path"
    )
    parser.add_argument(
        "--label-left", default="", help="Left label (ignored, for future use)"
    )
    parser.add_argument(
        "--label-right", default="", help="Right label (ignored, for future use)"
    )
    args = parser.parse_args()

    lw, lh, ldata = read_ppm(args.left)
    rw, rh, rdata = read_ppm(args.right)

    target_h = max(lh, rh)
    lnew_w = int(lw * target_h / lh)
    rnew_w = int(rw * target_h / rh)
    ldata = resize_nearest(lw, lh, ldata, lnew_w, target_h)
    rdata = resize_nearest(rw, rh, rdata, rnew_w, target_h)

    out_w = lnew_w + rnew_w
    out_h = target_h
    out = bytearray(out_w * out_h * 3)
    for y in range(out_h):
        for x in range(lnew_w):
            src_idx = (y * lnew_w + x) * 3
            dst_idx = (y * out_w + x) * 3
            out[dst_idx:dst_idx + 3] = ldata[src_idx:src_idx + 3]
        for x in range(rnew_w):
            src_idx = (y * rnew_w + x) * 3
            dst_idx = (y * out_w + lnew_w + x) * 3
            out[dst_idx:dst_idx + 3] = rdata[src_idx:src_idx + 3]

    write_ppm(args.output, out_w, out_h, bytes(out))
    print(f"Wrote {args.output} ({out_w} x {out_h})")


if __name__ == "__main__":
    main()

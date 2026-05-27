#!/usr/bin/env python3
"""Validate the Fortran ONETWO PRC kernel against PyMak's ONETWO routine."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import subprocess
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_BT_T = 6.0
DEFAULT_MINOR_RADIUS_M = 2.72
DEFAULT_MAJOR_RADIUS_M = 8.03
DEFAULT_WALL_REFLECTION = 0.8


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_pymak_onetwo(root: Path):
    """Load PyMak/cyclotron_radiation.py without importing PyMak.__init__."""

    pymak_pkg = types.ModuleType("PyMak")
    pymak_pkg.__path__ = [str(root / "PyMak")]

    gfile_reader_stub = types.ModuleType("PyMak.gfile_reader")
    gfile_reader_stub.build_cytran_input_arrays = _unused_pymak_stub
    gfile_reader_stub.map_psi_n_to_rho_tor = _unused_pymak_stub

    previous_pymak = sys.modules.get("PyMak")
    previous_gfile_reader = sys.modules.get("PyMak.gfile_reader")
    sys.modules["PyMak"] = pymak_pkg
    sys.modules["PyMak.gfile_reader"] = gfile_reader_stub

    try:
        module_path = root / "PyMak" / "cyclotron_radiation.py"
        spec = importlib.util.spec_from_file_location(
            "PyMak.cyclotron_radiation", module_path
        )
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not load module spec for {module_path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules["PyMak.cyclotron_radiation"] = module
        spec.loader.exec_module(module)
        return module.calculate_cyclotron_radiation_onetwo
    finally:
        if previous_pymak is None:
            sys.modules.pop("PyMak", None)
        else:
            sys.modules["PyMak"] = previous_pymak
        if previous_gfile_reader is None:
            sys.modules.pop("PyMak.gfile_reader", None)
        else:
            sys.modules["PyMak.gfile_reader"] = previous_gfile_reader


def _unused_pymak_stub(*args, **kwargs):
    raise RuntimeError("This validation only uses PyMak's ONETWO analytic routine.")


def run_command(command: list[str], cwd: Path) -> None:
    printable = " ".join(command)
    print(f"+ {printable}")
    subprocess.run(command, cwd=cwd, check=True)


def compile_fortran(root: Path, validation_dir: Path, build_dir: Path) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    bpsd_kinds = root.parent.parent / "bpsd" / "bpsd_kinds.f90"
    kernel = root / "tr_onetwo_prc_kernel.f90"
    driver = validation_dir / "onetwo_prc_driver.f90"

    objects = [
        build_dir / "bpsd_kinds.o",
        build_dir / "tr_onetwo_prc_kernel.o",
        build_dir / "onetwo_prc_driver.o",
    ]
    exe = build_dir / "onetwo_prc_driver"

    run_command(
        [
            "gfortran",
            "-c",
            str(bpsd_kinds),
            "-J",
            str(build_dir),
            "-o",
            str(objects[0]),
        ],
        cwd=root,
    )
    run_command(
        [
            "gfortran",
            "-c",
            str(kernel),
            "-I",
            str(build_dir),
            "-J",
            str(build_dir),
            "-o",
            str(objects[1]),
        ],
        cwd=root,
    )
    run_command(
        [
            "gfortran",
            "-c",
            str(driver),
            "-I",
            str(build_dir),
            "-J",
            str(build_dir),
            "-o",
            str(objects[2]),
        ],
        cwd=root,
    )
    run_command(["gfortran", *(str(obj) for obj in objects), "-o", str(exe)], cwd=root)
    return exe


def read_input_profile(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        "rho": np.atleast_1d(data["rho"]).astype(float),
        "ne": np.atleast_1d(data["ne_1e20_m3"]).astype(float),
        "te": np.atleast_1d(data["te_kev"]).astype(float),
    }


def read_fortran_output(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {
        "rho": np.atleast_1d(data["rho"]).astype(float),
        "ne": np.atleast_1d(data["ne_1e20_m3"]).astype(float),
        "te": np.atleast_1d(data["te_kev"]).astype(float),
        "prc_w_m3": np.atleast_1d(data["prc_w_m3"]).astype(float),
        "prc_mw_m3": np.atleast_1d(data["prc_mw_m3"]).astype(float),
        "phi_bar": np.atleast_1d(data["phi_bar"]).astype(float),
    }


def write_comparison_csv(
    path: Path,
    rho: np.ndarray,
    ne: np.ndarray,
    te: np.ndarray,
    fortran_prc_w_m3: np.ndarray,
    pymak_prc_w_m3: np.ndarray,
    fortran_phi: np.ndarray,
    pymak_phi: np.ndarray,
) -> None:
    abs_err_prc = fortran_prc_w_m3 - pymak_prc_w_m3
    rel_err_prc = np.abs(abs_err_prc) / np.maximum(np.abs(pymak_prc_w_m3), 1.0e-300)
    abs_err_phi = fortran_phi - pymak_phi
    rel_err_phi = np.abs(abs_err_phi) / np.maximum(np.abs(pymak_phi), 1.0e-300)

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "rho",
                "ne_1e20_m3",
                "te_kev",
                "fortran_prc_w_m3",
                "pymak_prc_w_m3",
                "abs_err_prc_w_m3",
                "rel_err_prc",
                "fortran_phi_bar",
                "pymak_phi_bar",
                "abs_err_phi_bar",
                "rel_err_phi_bar",
            ]
        )
        for row in zip(
            rho,
            ne,
            te,
            fortran_prc_w_m3,
            pymak_prc_w_m3,
            abs_err_prc,
            rel_err_prc,
            fortran_phi,
            pymak_phi,
            abs_err_phi,
            rel_err_phi,
        ):
            writer.writerow([f"{value:.16e}" for value in row])


def plot_comparison(
    path: Path,
    rho: np.ndarray,
    fortran_prc_w_m3: np.ndarray,
    pymak_prc_w_m3: np.ndarray,
    fortran_phi: np.ndarray,
    pymak_phi: np.ndarray,
) -> None:
    prc_rel_err = np.abs(fortran_prc_w_m3 - pymak_prc_w_m3) / np.maximum(
        np.abs(pymak_prc_w_m3), 1.0e-300
    )
    phi_rel_err = np.abs(fortran_phi - pymak_phi) / np.maximum(
        np.abs(pymak_phi), 1.0e-300
    )

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)

    axes[0, 0].plot(rho, pymak_prc_w_m3 * 1.0e-6, "o-", label="PyMak")
    axes[0, 0].plot(rho, fortran_prc_w_m3 * 1.0e-6, "s--", label="Fortran")
    axes[0, 0].set_xlabel("rho")
    axes[0, 0].set_ylabel("PRC [MW/m^3]")
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].semilogy(rho, np.maximum(prc_rel_err, 1.0e-18), "o-")
    axes[0, 1].set_xlabel("rho")
    axes[0, 1].set_ylabel("relative error PRC")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(rho, pymak_phi, "o-", label="PyMak")
    axes[1, 0].plot(rho, fortran_phi, "s--", label="Fortran")
    axes[1, 0].set_xlabel("rho")
    axes[1, 0].set_ylabel("phi_bar")
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].semilogy(rho, np.maximum(phi_rel_err, 1.0e-18), "o-")
    axes[1, 1].set_xlabel("rho")
    axes[1, 1].set_ylabel("relative error phi_bar")
    axes[1, 1].grid(True, alpha=0.3)

    fig.suptitle("ONETWO PRC kernel validation")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).with_name("input_profile.csv"),
        help="Input CSV with rho, ne_1e20_m3, te_kev columns.",
    )
    parser.add_argument("--bt-T", type=float, default=DEFAULT_BT_T)
    parser.add_argument("--minor-radius-m", type=float, default=DEFAULT_MINOR_RADIUS_M)
    parser.add_argument("--major-radius-m", type=float, default=DEFAULT_MAJOR_RADIUS_M)
    parser.add_argument("--wall-reflection", type=float, default=DEFAULT_WALL_REFLECTION)
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=Path(__file__).with_name("build"),
        help="Directory for compiled Fortran validation objects.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).with_name("results"),
        help="Directory for validation CSV and PNG outputs.",
    )
    args = parser.parse_args()

    root = repo_root()
    validation_dir = Path(__file__).resolve().parent
    args.output_dir.mkdir(parents=True, exist_ok=True)

    exe = compile_fortran(root, validation_dir, args.build_dir)
    fortran_csv = args.output_dir / "fortran_onetwo_prc.csv"
    comparison_csv = args.output_dir / "onetwo_prc_comparison.csv"
    figure_path = args.output_dir / "onetwo_prc_comparison.png"

    run_command(
        [
            str(exe),
            str(args.input),
            str(fortran_csv),
            f"{args.bt_T:.17g}",
            f"{args.minor_radius_m:.17g}",
            f"{args.major_radius_m:.17g}",
            f"{args.wall_reflection:.17g}",
        ],
        cwd=root,
    )

    profile = read_input_profile(args.input)
    fortran = read_fortran_output(fortran_csv)

    calculate_onetwo = load_pymak_onetwo(root)
    geometry = SimpleNamespace(bcentr=args.bt_T, rmaxis=args.major_radius_m)
    pymak = calculate_onetwo(
        geometry=geometry,
        ne=profile["ne"],
        Te=profile["te"],
        bt_T=args.bt_T,
        wall_reflection=args.wall_reflection,
        minor_radius_m=args.minor_radius_m,
        major_radius_m=args.major_radius_m,
    )
    pymak_prc_w_m3 = np.asarray(pymak["power_density_total_MW_m3"], dtype=float) * 1.0e6
    pymak_phi = np.asarray(pymak["phi_bar"], dtype=float)

    write_comparison_csv(
        comparison_csv,
        profile["rho"],
        profile["ne"],
        profile["te"],
        fortran["prc_w_m3"],
        pymak_prc_w_m3,
        fortran["phi_bar"],
        pymak_phi,
    )
    plot_comparison(
        figure_path,
        profile["rho"],
        fortran["prc_w_m3"],
        pymak_prc_w_m3,
        fortran["phi_bar"],
        pymak_phi,
    )

    prc_rel_err = np.abs(fortran["prc_w_m3"] - pymak_prc_w_m3) / np.maximum(
        np.abs(pymak_prc_w_m3), 1.0e-300
    )
    phi_rel_err = np.abs(fortran["phi_bar"] - pymak_phi) / np.maximum(
        np.abs(pymak_phi), 1.0e-300
    )

    print(f"Wrote {fortran_csv}")
    print(f"Wrote {comparison_csv}")
    print(f"Wrote {figure_path}")
    print(f"max_rel_err_prc={np.max(prc_rel_err):.6e}")
    print(f"max_rel_err_phi_bar={np.max(phi_rel_err):.6e}")


if __name__ == "__main__":
    main()

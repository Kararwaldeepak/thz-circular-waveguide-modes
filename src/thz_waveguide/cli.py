"""Command-line interface for the broadband waveguide simulation."""

from __future__ import annotations

import argparse
from pathlib import Path

from .model import Lens, Waveguide
from .plotting import make_all_plots
from .simulation import export_results, run_sweep


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Semi-analytical Gaussian coupling to a circular copper THz waveguide."
    )
    parser.add_argument("--diameter-mm", type=float, default=2.0)
    parser.add_argument("--length-cm", type=float, default=40.0)
    parser.add_argument("--conductivity", type=float, default=5.8e7, help="Copper S/m")
    parser.add_argument("--focal-length-mm", type=float, default=100.0)
    parser.add_argument("--lens-diameter-mm", type=float, default=50.8)
    parser.add_argument(
        "--fill-factor",
        type=float,
        default=1.0,
        help="1/e beam radius on lens divided by lens radius",
    )
    parser.add_argument("--minimum-thz", type=float, default=0.10)
    parser.add_argument("--maximum-thz", type=float, default=2.00)
    parser.add_argument("--points", type=int, default=96)
    parser.add_argument("--radial-orders", type=int, default=12)
    parser.add_argument("--grid-size", type=int, default=161)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    waveguide = Waveguide(
        radius_m=args.diameter_mm * 1e-3 / 2.0,
        length_m=args.length_cm * 1e-2,
        conductivity_s_m=args.conductivity,
    )
    lens = Lens(
        focal_length_m=args.focal_length_mm * 1e-3,
        diameter_m=args.lens_diameter_mm * 1e-3,
        fill_factor=args.fill_factor,
    )
    result = run_sweep(
        waveguide=waveguide,
        lens=lens,
        minimum_thz=args.minimum_thz,
        maximum_thz=args.maximum_thz,
        points=args.points,
        radial_orders=args.radial_orders,
        grid_size=args.grid_size,
    )
    export_results(result, args.output)
    if not args.skip_plots:
        make_all_plots(result, args.output)

    te11_cutoff = result.modes[0].cutoff_hz(waveguide) / 1e12
    print(f"Completed {len(result.frequencies_hz)} frequencies.")
    print(
        f"Waveguide: {waveguide.diameter_m * 1e3:.3f} mm ID "
        f"x {waveguide.length_m * 100:.1f} cm"
    )
    print(f"Lens: f/{lens.f_number:.2f}, NA={lens.numerical_aperture:.3f}")
    print(f"TE11 cutoff: {te11_cutoff:.5f} THz")
    print(f"Results: {args.output.resolve()}")


if __name__ == "__main__":
    main()

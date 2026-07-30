"""Broadband simulation orchestration and CSV export."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .model import (
    Z0,
    Lens,
    ModeFields,
    ModeSpec,
    Waveguide,
    aperture_throughput,
    conductor_attenuation,
    coupling_amplitude,
    focused_waist,
    gaussian_field,
    mode_fields,
    power_loss_db,
    relative_group_delay_s,
)


@dataclass
class ModalPoint:
    frequency_thz: float
    mode: str
    cutoff_thz: float
    beta_rad_m: float
    coupling_fraction_full_beam: float
    coupling_fraction_aperture: float
    attenuation_np_m: float
    copper_loss_db: float
    output_power_fraction: float
    relative_group_delay_ps: float
    amplitude_real: float
    amplitude_imag: float


@dataclass
class SpectrumPoint:
    frequency_thz: float
    waist_mm: float
    aperture_throughput: float
    total_coupled_fraction: float
    total_output_power_fraction: float
    dominant_mode: str
    dominant_output_fraction: float


@dataclass
class SimulationResult:
    waveguide: Waveguide
    lens: Lens
    frequencies_hz: np.ndarray
    modes: list[ModeSpec]
    modal_rows: list[ModalPoint]
    spectrum_rows: list[SpectrumPoint]


def coupling_modes(radial_orders: int = 12) -> list[ModeSpec]:
    """Modes allowed by a centered, x-polarized, cylindrically symmetric beam."""

    modes: list[ModeSpec] = []
    for n in range(1, radial_orders + 1):
        # These orientations have predominantly x-directed transverse E.
        modes.append(ModeSpec("TE", 1, n, "sin"))
        modes.append(ModeSpec("TM", 1, n, "cos"))
    return modes


def _full_gaussian_power(waist_m: float) -> float:
    """Power for peak electric field 1 V/m, integrated over the infinite plane."""

    return float(np.pi * waist_m**2 / (4.0 * Z0))


def simulate_frequency(
    frequency_hz: float,
    waveguide: Waveguide,
    lens: Lens,
    modes: list[ModeSpec],
    grid_size: int,
) -> tuple[list[ModalPoint], SpectrumPoint]:
    waist = focused_waist(frequency_hz, lens)
    throughput = aperture_throughput(waist, waveguide)
    full_power = _full_gaussian_power(waist)
    aperture_power = full_power * throughput

    working: list[dict[str, object]] = []
    for mode in modes:
        if not mode.propagates(frequency_hz, waveguide):
            continue
        fields = mode_fields(mode, frequency_hz, waveguide, grid_size=grid_size)
        inc_ex, inc_ey = gaussian_field(fields, waist)
        amplitude, coupled_power = coupling_amplitude(fields, inc_ex, inc_ey)
        working.append(
            {
                "mode": mode,
                "fields": fields,
                "amplitude": amplitude,
                "coupled_power": coupled_power,
            }
        )

    # Remove sub-percent quadrature excess when a finite grid approximates a
    # complete orthogonal expansion.  The modal input power cannot exceed the
    # Gaussian power geometrically admitted by the bore.
    raw_total = sum(float(item["coupled_power"]) for item in working)
    energy_scale = min(1.0, aperture_power / raw_total) if raw_total > 0.0 else 1.0

    rows: list[ModalPoint] = []
    dominant_mode = "none"
    dominant_output = 0.0
    total_coupled = 0.0
    total_output = 0.0
    for item in working:
        mode = item["mode"]
        fields = item["fields"]
        assert isinstance(mode, ModeSpec)
        assert isinstance(fields, ModeFields)
        amplitude = complex(item["amplitude"]) * np.sqrt(energy_scale)
        coupled_power = float(item["coupled_power"]) * energy_scale
        alpha = conductor_attenuation(mode, frequency_hz, waveguide, fields=fields)
        output_power = coupled_power * np.exp(-2.0 * alpha * waveguide.length_m)
        input_fraction = coupled_power / full_power
        aperture_fraction = coupled_power / aperture_power if aperture_power > 0.0 else 0.0
        output_fraction = output_power / full_power
        total_coupled += input_fraction
        total_output += output_fraction
        if output_fraction > dominant_output:
            dominant_output = output_fraction
            dominant_mode = mode.label
        rows.append(
            ModalPoint(
                frequency_thz=frequency_hz / 1e12,
                mode=mode.label,
                cutoff_thz=mode.cutoff_hz(waveguide) / 1e12,
                beta_rad_m=fields.beta_rad_m,
                coupling_fraction_full_beam=input_fraction,
                coupling_fraction_aperture=aperture_fraction,
                attenuation_np_m=alpha,
                copper_loss_db=power_loss_db(alpha, waveguide.length_m),
                output_power_fraction=output_fraction,
                relative_group_delay_ps=relative_group_delay_s(
                    mode, frequency_hz, waveguide
                )
                * 1e12,
                amplitude_real=amplitude.real,
                amplitude_imag=amplitude.imag,
            )
        )

    spectrum = SpectrumPoint(
        frequency_thz=frequency_hz / 1e12,
        waist_mm=waist * 1e3,
        aperture_throughput=throughput,
        total_coupled_fraction=total_coupled,
        total_output_power_fraction=total_output,
        dominant_mode=dominant_mode,
        dominant_output_fraction=dominant_output,
    )
    return rows, spectrum


def run_sweep(
    waveguide: Waveguide,
    lens: Lens,
    minimum_thz: float = 0.10,
    maximum_thz: float = 2.00,
    points: int = 96,
    radial_orders: int = 12,
    grid_size: int = 161,
) -> SimulationResult:
    frequencies_hz = np.linspace(minimum_thz, maximum_thz, points) * 1e12
    modes = coupling_modes(radial_orders)
    modal_rows: list[ModalPoint] = []
    spectrum_rows: list[SpectrumPoint] = []
    for frequency_hz in frequencies_hz:
        modal, spectrum = simulate_frequency(
            frequency_hz, waveguide, lens, modes, grid_size
        )
        modal_rows.extend(modal)
        spectrum_rows.append(spectrum)
    return SimulationResult(
        waveguide=waveguide,
        lens=lens,
        frequencies_hz=frequencies_hz,
        modes=modes,
        modal_rows=modal_rows,
        spectrum_rows=spectrum_rows,
    )


def cutoff_table(
    waveguide: Waveguide, max_m: int = 4, radial_orders: int = 4
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for family in ("TE", "TM"):
        for m in range(max_m + 1):
            for n in range(1, radial_orders + 1):
                mode = ModeSpec(family, m, n, "cos")
                rows.append(
                    {
                        "mode": mode.label,
                        "family": family,
                        "m": m,
                        "n": n,
                        "bessel_root": mode.root,
                        "cutoff_thz": mode.cutoff_hz(waveguide) / 1e12,
                    }
                )
    return sorted(rows, key=lambda row: float(row["cutoff_thz"]))


def _write_dataclass_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def export_results(result: SimulationResult, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_dataclass_csv(output_dir / "modal_spectrum.csv", result.modal_rows)
    _write_dataclass_csv(output_dir / "total_spectrum.csv", result.spectrum_rows)

    cutoffs = cutoff_table(result.waveguide)
    with (output_dir / "cutoff_frequencies.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cutoffs[0]))
        writer.writeheader()
        writer.writerows(cutoffs)

    metadata = {
        "model": "semi-analytical circular metallic waveguide",
        "phasor_convention": "exp(+j*omega*t - j*beta*z)",
        "waveguide": asdict(result.waveguide),
        "lens": asdict(result.lens),
        "frequency_min_thz": float(result.frequencies_hz.min() / 1e12),
        "frequency_max_thz": float(result.frequencies_hz.max() / 1e12),
        "frequency_points": len(result.frequencies_hz),
        "coupling_modes": [mode.oriented_label for mode in result.modes],
        "assumptions": [
            "air-filled circular bore",
            "perfect-conductor eigenfields",
            "copper loss from first-order surface-resistance perturbation",
            "centered x-polarized diffraction-limited Gaussian",
            "lens fill_factor defines the 1/e field radius on the lens",
            "no flange reflection, roughness, taper, misalignment, or atmospheric absorption",
        ],
    }
    (output_dir / "simulation_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

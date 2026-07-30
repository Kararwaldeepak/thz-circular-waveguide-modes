"""Publication-style plots for the broadband waveguide simulation."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from .model import (
    Lens,
    ModeSpec,
    Waveguide,
    conductor_attenuation,
    coupling_amplitude,
    focused_waist,
    gaussian_field,
    mode_fields,
)
from .simulation import SimulationResult, coupling_modes


COLORS = {
    "TE11": "#0068B4",
    "TM11": "#D1495B",
    "TE12": "#00A676",
    "TM12": "#F28E2B",
    "TE13": "#8E5EA2",
    "TM13": "#8C564B",
}


def _style_axis(axis: plt.Axes) -> None:
    axis.tick_params(which="both", direction="in", top=True, right=True)
    axis.minorticks_on()
    axis.grid(alpha=0.22, linewidth=0.6)


def plot_cutoffs(result: SimulationResult, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.7), constrained_layout=True)
    displayed = sorted(
        {
            (mode.label, mode.cutoff_hz(result.waveguide) / 1e12)
            for mode in result.modes
            if mode.cutoff_hz(result.waveguide) <= result.frequencies_hz.max()
        },
        key=lambda item: item[1],
    )
    for index, (label, cutoff) in enumerate(displayed):
        color = COLORS.get(label, "0.55")
        ax.vlines(cutoff, index - 0.36, index + 0.36, color=color, linewidth=2.0)
    ax.set_yticks(range(len(displayed)), [label for label, _ in displayed])
    ax.set_xlim(0.0, result.frequencies_hz.max() / 1e12)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Gaussian-accessible modes")
    ax.set_title("Cutoff frequencies for centered x-polarized excitation")
    _style_axis(ax)
    fig.savefig(output_dir / "01_mode_cutoffs.png", dpi=220)
    plt.close(fig)


def plot_coupling_and_output(result: SimulationResult, output_dir: Path) -> None:
    frequency = np.array([row.frequency_thz for row in result.spectrum_rows])
    throughput = np.array([row.aperture_throughput for row in result.spectrum_rows])
    coupled = np.array([row.total_coupled_fraction for row in result.spectrum_rows])
    output = np.array([row.total_output_power_fraction for row in result.spectrum_rows])

    fig, ax = plt.subplots(figsize=(8.0, 4.7), constrained_layout=True)
    ax.plot(frequency, throughput, color="0.25", linestyle="--", label="Aperture throughput")
    ax.plot(frequency, coupled, color="#0068B4", label="Power coupled into modes")
    ax.plot(frequency, output, color="#D1495B", label="Power after 40 cm")
    ax.set(xlabel="Frequency (THz)", ylabel="Fraction of incident Gaussian power", ylim=(0, 1.05))
    ax.set_title("Gaussian coupling and propagated modal power")
    ax.legend(frameon=False, loc="lower right")
    _style_axis(ax)
    fig.savefig(output_dir / "02_coupling_and_output.png", dpi=220)
    plt.close(fig)


def plot_modal_power(result: SimulationResult, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.7), constrained_layout=True)
    for label in ("TE11", "TM11", "TE12", "TM12", "TE13", "TM13"):
        rows = [row for row in result.modal_rows if row.mode == label]
        if not rows:
            continue
        ax.plot(
            [row.frequency_thz for row in rows],
            [row.output_power_fraction for row in rows],
            label=label,
            color=COLORS[label],
        )
    ax.set(
        xlabel="Frequency (THz)",
        ylabel="Output power / incident Gaussian power",
        ylim=(0, None),
    )
    ax.set_title("Modal content at the waveguide output")
    ax.legend(frameon=False, ncol=3)
    _style_axis(ax)
    fig.savefig(output_dir / "03_output_modal_power.png", dpi=220)
    plt.close(fig)


def plot_copper_loss(result: SimulationResult, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.7), constrained_layout=True)
    for label in ("TE11", "TM11", "TE12", "TM12"):
        rows = [row for row in result.modal_rows if row.mode == label]
        if not rows:
            continue
        ax.plot(
            [row.frequency_thz for row in rows],
            [row.copper_loss_db for row in rows],
            label=label,
            color=COLORS[label],
        )
    ax.set(xlabel="Frequency (THz)", ylabel="Copper-wall loss over 40 cm (dB)")
    ax.set_title("Perturbative conductor loss")
    ax.legend(frameon=False, ncol=2)
    _style_axis(ax)
    fig.savefig(output_dir / "04_copper_loss.png", dpi=220)
    plt.close(fig)


def _output_field(
    frequency_hz: float,
    waveguide: Waveguide,
    lens: Lens,
    grid_size: int,
    radial_orders: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ex_total: np.ndarray | None = None
    ey_total: np.ndarray | None = None
    axis: np.ndarray | None = None
    mask: np.ndarray | None = None
    for mode in coupling_modes(radial_orders):
        if not mode.propagates(frequency_hz, waveguide):
            continue
        fields = mode_fields(mode, frequency_hz, waveguide, grid_size)
        waist = focused_waist(frequency_hz, lens)
        inc_ex, inc_ey = gaussian_field(fields, waist)
        amplitude, _ = coupling_amplitude(fields, inc_ex, inc_ey)
        alpha = conductor_attenuation(mode, frequency_hz, waveguide, fields)
        propagated = amplitude * np.exp(
            -alpha * waveguide.length_m - 1j * fields.beta_rad_m * waveguide.length_m
        )
        if ex_total is None:
            ex_total = np.zeros_like(fields.ex)
            ey_total = np.zeros_like(fields.ey)
            axis = fields.x_m
            mask = fields.mask
        ex_total += propagated * fields.ex
        ey_total += propagated * fields.ey
    if ex_total is None or ey_total is None or axis is None or mask is None:
        raise ValueError("No propagating mode at the requested frequency")

    # Choose a global phase that makes the strongest pixel instantaneously real.
    magnitude = np.abs(ex_total) ** 2 + np.abs(ey_total) ** 2
    peak = np.unravel_index(np.argmax(magnitude), magnitude.shape)
    reference = ex_total[peak] if abs(ex_total[peak]) >= abs(ey_total[peak]) else ey_total[peak]
    phase = np.exp(-1j * np.angle(reference))
    return axis, mask, ex_total * phase, ey_total * phase, magnitude


def plot_output_fields(
    result: SimulationResult,
    output_dir: Path,
    selected_thz: tuple[float, ...] = (0.20, 0.50, 1.00, 1.50),
    grid_size: int = 181,
    radial_orders: int = 12,
) -> None:
    fig, axes = plt.subplots(
        1, len(selected_thz), figsize=(3.6 * len(selected_thz), 3.7), constrained_layout=True
    )
    for ax, frequency_thz in zip(np.atleast_1d(axes), selected_thz):
        axis, mask, ex, ey, intensity = _output_field(
            frequency_thz * 1e12,
            result.waveguide,
            result.lens,
            grid_size,
            radial_orders,
        )
        normalized = intensity / np.nanmax(intensity)
        normalized[~mask] = np.nan
        extent = [axis[0] * 1e3, axis[-1] * 1e3, axis[0] * 1e3, axis[-1] * 1e3]
        image = ax.imshow(
            normalized,
            origin="lower",
            extent=extent,
            cmap="magma",
            vmin=0.0,
            vmax=1.0,
        )
        step = max(1, grid_size // 17)
        sampled = np.s_[::step, ::step]
        x, y = np.meshgrid(axis * 1e3, axis * 1e3)
        u, v = np.real(ex), np.real(ey)
        norm = np.hypot(u, v)
        valid = mask & (norm > 0.08 * np.nanmax(norm))
        u_plot = np.where(valid, u / np.maximum(norm, 1e-30), np.nan)
        v_plot = np.where(valid, v / np.maximum(norm, 1e-30), np.nan)
        ax.quiver(
            x[sampled],
            y[sampled],
            u_plot[sampled],
            v_plot[sampled],
            color="white",
            pivot="mid",
            scale=20,
            width=0.006,
        )
        ax.add_patch(Circle((0, 0), result.waveguide.radius_m * 1e3, fill=False, color="white"))
        ax.set(
            title=f"{frequency_thz:.2f} THz",
            xlabel="x (mm)",
            ylabel="y (mm)",
            aspect="equal",
        )
        ax.tick_params(direction="in", top=True, right=True)
    fig.colorbar(image, ax=np.atleast_1d(axes), label="Normalized output intensity", shrink=0.82)
    fig.suptitle("Coherent output fields after 40 cm (arrows: instantaneous transverse E)")
    fig.savefig(output_dir / "05_output_field_maps.png", dpi=220)
    plt.close(fig)


def make_all_plots(result: SimulationResult, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_cutoffs(result, output_dir)
    plot_coupling_and_output(result, output_dir)
    plot_modal_power(result, output_dir)
    plot_copper_loss(result, output_dir)
    plot_output_fields(result, output_dir)

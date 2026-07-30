"""Circular-waveguide fields, copper loss, and Gaussian coupling.

The phasor convention is exp(+j*omega*t - j*beta*z).  Perfect-conductor
eigenfields are used for the modal profiles.  Finite copper conductivity is
included perturbatively through the surface-resistance wall-loss integral.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.constants import c as C0
from scipy.constants import epsilon_0 as EPS0
from scipy.constants import mu_0 as MU0
from scipy.special import jn_zeros, jnp_zeros, jv, jvp

Z0 = float(np.sqrt(MU0 / EPS0))
Family = Literal["TE", "TM"]
Orientation = Literal["cos", "sin"]
ComplexArray = NDArray[np.complex128]
RealArray = NDArray[np.float64]


@dataclass(frozen=True)
class Waveguide:
    """Air-filled circular copper waveguide."""

    radius_m: float = 1.0e-3
    length_m: float = 0.40
    conductivity_s_m: float = 5.8e7

    @property
    def diameter_m(self) -> float:
        return 2.0 * self.radius_m


@dataclass(frozen=True)
class Lens:
    """Gaussian focusing lens.

    ``fill_factor`` is w_lens / (D/2), where w_lens is the 1/e electric-field
    radius of the collimated Gaussian on the lens.
    """

    focal_length_m: float = 0.100
    diameter_m: float = 0.0508
    fill_factor: float = 1.0
    beam_quality_m2: float = 1.0

    @property
    def beam_radius_on_lens_m(self) -> float:
        return self.fill_factor * self.diameter_m / 2.0

    @property
    def f_number(self) -> float:
        return self.focal_length_m / self.diameter_m

    @property
    def numerical_aperture(self) -> float:
        half_angle = np.arctan2(self.diameter_m / 2.0, self.focal_length_m)
        return float(np.sin(half_angle))


@dataclass(frozen=True)
class ModeSpec:
    """One TE or TM circular-waveguide mode and one angular orientation."""

    family: Family
    m: int
    n: int
    orientation: Orientation = "cos"

    def __post_init__(self) -> None:
        family = self.family.upper()
        if family not in {"TE", "TM"}:
            raise ValueError("family must be 'TE' or 'TM'")
        object.__setattr__(self, "family", family)
        if self.m < 0 or self.n < 1:
            raise ValueError("Require m >= 0 and n >= 1")
        if self.orientation not in {"cos", "sin"}:
            raise ValueError("orientation must be 'cos' or 'sin'")
        if self.m == 0 and self.orientation == "sin":
            raise ValueError("The sin orientation vanishes for m = 0")

    @cached_property
    def root(self) -> float:
        zeros = jnp_zeros(self.m, self.n) if self.family == "TE" else jn_zeros(self.m, self.n)
        return float(zeros[-1])

    @property
    def label(self) -> str:
        radial = str(self.n) if self.n < 10 else f",{self.n}"
        return f"{self.family}{self.m}{radial}"

    @property
    def oriented_label(self) -> str:
        return f"{self.label}-{self.orientation}"

    def cutoff_hz(self, waveguide: Waveguide) -> float:
        return C0 * self.root / (2.0 * np.pi * waveguide.radius_m)

    def propagates(self, frequency_hz: float, waveguide: Waveguide) -> bool:
        return frequency_hz > self.cutoff_hz(waveguide)


@dataclass
class ModeFields:
    """Complex phasor fields on a Cartesian transverse grid."""

    x_m: RealArray
    y_m: RealArray
    mask: NDArray[np.bool_]
    ex: ComplexArray
    ey: ComplexArray
    ez: ComplexArray
    hx: ComplexArray
    hy: ComplexArray
    hz: ComplexArray
    beta_rad_m: float
    power_w: float

    @property
    def dx_m(self) -> float:
        return float(self.x_m[1] - self.x_m[0])


def focused_waist(frequency_hz: float, lens: Lens) -> float:
    """Return the paraxial 1/e electric-field waist at the lens focus."""

    if frequency_hz <= 0.0:
        raise ValueError("frequency_hz must be positive")
    if lens.beam_radius_on_lens_m <= 0.0:
        raise ValueError("The illuminated lens radius must be positive")
    wavelength = C0 / frequency_hz
    return (
        lens.beam_quality_m2
        * wavelength
        * lens.focal_length_m
        / (np.pi * lens.beam_radius_on_lens_m)
    )


def aperture_throughput(waist_m: float, waveguide: Waveguide) -> float:
    """Gaussian power fraction geometrically incident on the circular aperture."""

    if waist_m <= 0.0:
        raise ValueError("waist_m must be positive")
    return float(1.0 - np.exp(-2.0 * (waveguide.radius_m / waist_m) ** 2))


def propagation_constant(
    mode: ModeSpec, frequency_hz: float, waveguide: Waveguide
) -> float:
    """Return real beta for a propagating mode."""

    k0 = 2.0 * np.pi * frequency_hz / C0
    kc = mode.root / waveguide.radius_m
    if k0 <= kc:
        raise ValueError(f"{mode.label} is below cutoff at {frequency_hz / 1e12:.4g} THz")
    return float(np.sqrt(k0**2 - kc**2))


def group_velocity(mode: ModeSpec, frequency_hz: float, waveguide: Waveguide) -> float:
    """Lossless-guide group velocity."""

    fc = mode.cutoff_hz(waveguide)
    if frequency_hz <= fc:
        raise ValueError(f"{mode.label} is below cutoff")
    return float(C0 * np.sqrt(1.0 - (fc / frequency_hz) ** 2))


def _angular_terms(
    m: int, phi: RealArray, orientation: Orientation
) -> tuple[RealArray, RealArray]:
    if orientation == "cos":
        return np.cos(m * phi), -m * np.sin(m * phi)
    return np.sin(m * phi), m * np.cos(m * phi)


def _spatial_terms(
    mode: ModeSpec,
    radius: RealArray,
    phi: RealArray,
    waveguide: Waveguide,
) -> tuple[RealArray, RealArray, RealArray]:
    """Return psi, d(psi)/dr, and (1/r)d(psi)/dphi."""

    kc = mode.root / waveguide.radius_m
    angular, d_angular = _angular_terms(mode.m, phi, mode.orientation)
    kr = kc * radius
    radial = jv(mode.m, kr)
    psi = radial * angular
    dpsi_dr = kc * jvp(mode.m, kr) * angular
    dpsi_dphi_over_r = np.zeros_like(radius)
    nonzero = radius > 0.0
    dpsi_dphi_over_r[nonzero] = (
        radial[nonzero] * d_angular[nonzero] / radius[nonzero]
    )
    return psi, dpsi_dr, dpsi_dphi_over_r


def mode_fields(
    mode: ModeSpec,
    frequency_hz: float,
    waveguide: Waveguide,
    grid_size: int = 181,
) -> ModeFields:
    """Evaluate a propagating mode on a square grid spanning the circular bore."""

    if grid_size < 41:
        raise ValueError("grid_size must be at least 41")
    a = waveguide.radius_m
    axis = np.linspace(-a, a, grid_size, dtype=float)
    x, y = np.meshgrid(axis, axis)
    radius = np.hypot(x, y)
    phi = np.arctan2(y, x)
    mask = radius <= a

    psi, dpsi_dr, dpsi_dphi_over_r = _spatial_terms(mode, radius, phi, waveguide)
    kc = mode.root / a
    beta = propagation_constant(mode, frequency_hz, waveguide)
    omega = 2.0 * np.pi * frequency_hz

    # Convert grad_t(psi) from cylindrical to Cartesian components.
    grad_x = dpsi_dr * np.cos(phi) - dpsi_dphi_over_r * np.sin(phi)
    grad_y = dpsi_dr * np.sin(phi) + dpsi_dphi_over_r * np.cos(phi)

    # The polar formula has a removable singularity at r = 0 for m = 1.
    at_center = radius == 0.0
    if mode.m == 1:
        if mode.orientation == "cos":
            grad_x[at_center], grad_y[at_center] = kc / 2.0, 0.0
        else:
            grad_x[at_center], grad_y[at_center] = 0.0, kc / 2.0
    else:
        grad_x[at_center], grad_y[at_center] = 0.0, 0.0

    zeros = np.zeros_like(radius, dtype=np.complex128)
    if mode.family == "TE":
        # E_t = j omega mu / kc^2 (z-hat x grad_t H_z)
        # H_t = -j beta / kc^2 grad_t H_z
        ex = -1j * omega * MU0 * grad_y / kc**2
        ey = +1j * omega * MU0 * grad_x / kc**2
        ez = zeros.copy()
        hx = -1j * beta * grad_x / kc**2
        hy = -1j * beta * grad_y / kc**2
        hz = psi.astype(np.complex128)
    else:
        # E_t = -j beta / kc^2 grad_t E_z
        # H_t = -j omega eps / kc^2 (z-hat x grad_t E_z)
        ex = -1j * beta * grad_x / kc**2
        ey = -1j * beta * grad_y / kc**2
        ez = psi.astype(np.complex128)
        hx = +1j * omega * EPS0 * grad_y / kc**2
        hy = -1j * omega * EPS0 * grad_x / kc**2
        hz = zeros.copy()

    arrays = (ex, ey, ez, hx, hy, hz)
    for array in arrays:
        array[~mask] = 0.0

    dx = float(axis[1] - axis[0])
    poynting_z = 0.5 * np.real(ex * np.conj(hy) - ey * np.conj(hx))
    power = float(np.sum(poynting_z) * dx**2)
    if power <= 0.0:
        raise RuntimeError(f"Computed non-positive power for {mode.oriented_label}")

    return ModeFields(
        x_m=axis,
        y_m=axis,
        mask=mask,
        ex=ex,
        ey=ey,
        ez=ez,
        hx=hx,
        hy=hy,
        hz=hz,
        beta_rad_m=beta,
        power_w=power,
    )


def gaussian_field(
    fields: ModeFields, waist_m: float
) -> tuple[ComplexArray, ComplexArray]:
    """Centered, x-polarized Gaussian field clipped by the guide aperture."""

    x, y = np.meshgrid(fields.x_m, fields.y_m)
    radius2 = x**2 + y**2
    ex = np.exp(-radius2 / waist_m**2).astype(np.complex128)
    ex[~fields.mask] = 0.0
    ey = np.zeros_like(ex)
    return ex, ey


def coupling_amplitude(
    fields: ModeFields, incident_ex: ComplexArray, incident_ey: ComplexArray
) -> tuple[complex, float]:
    """Project an aperture field onto one forward waveguide mode.

    Returns the complex modal amplitude and coupled power.  The projection is
    a = integral[(E_inc x H_mode*) . z] / (2 P_mode).
    """

    if incident_ex.shape != fields.ex.shape or incident_ey.shape != fields.ey.shape:
        raise ValueError("Incident and modal fields must share the same grid")
    overlap = np.sum(
        incident_ex * np.conj(fields.hy) - incident_ey * np.conj(fields.hx)
    ) * fields.dx_m**2
    amplitude = complex(overlap / (2.0 * fields.power_w))
    coupled_power = float(abs(amplitude) ** 2 * fields.power_w)
    return amplitude, coupled_power


def surface_resistance(frequency_hz: float, conductivity_s_m: float) -> float:
    """Good-conductor surface resistance in ohms."""

    return float(np.sqrt(np.pi * frequency_hz * MU0 / conductivity_s_m))


def conductor_attenuation(
    mode: ModeSpec,
    frequency_hz: float,
    waveguide: Waveguide,
    fields: ModeFields | None = None,
    boundary_points: int = 1440,
) -> float:
    """Return the field attenuation constant alpha in Np/m.

    The wall power dissipated per unit length is
    (Rs/2) integral |H_t|^2 ds, and alpha = P_loss/(2 P_mode).
    """

    if fields is None:
        fields = mode_fields(mode, frequency_hz, waveguide)
    a = waveguide.radius_m
    kc = mode.root / a
    beta = fields.beta_rad_m
    omega = 2.0 * np.pi * frequency_hz
    phi = np.linspace(0.0, 2.0 * np.pi, boundary_points, endpoint=False)
    radius = np.full_like(phi, a)
    psi, dpsi_dr, dpsi_dphi_over_r = _spatial_terms(mode, radius, phi, waveguide)

    if mode.family == "TE":
        h_phi = -1j * beta * dpsi_dphi_over_r / kc**2
        h_z = psi
    else:
        h_phi = -1j * omega * EPS0 * dpsi_dr / kc**2
        h_z = np.zeros_like(h_phi)

    rs = surface_resistance(frequency_hz, waveguide.conductivity_s_m)
    dphi = 2.0 * np.pi / boundary_points
    wall_integral = a * dphi * np.sum(np.abs(h_phi) ** 2 + np.abs(h_z) ** 2)
    loss_per_length = 0.5 * rs * wall_integral
    return float(loss_per_length / (2.0 * fields.power_w))


def power_loss_db(alpha_np_m: float, length_m: float) -> float:
    """Power loss in dB for a field attenuation constant alpha."""

    return float(20.0 / np.log(10.0) * alpha_np_m * length_m)


def relative_group_delay_s(
    mode: ModeSpec, frequency_hz: float, waveguide: Waveguide
) -> float:
    """Mode delay relative to the same length in free space."""

    vg = group_velocity(mode, frequency_hz, waveguide)
    return float(waveguide.length_m * (1.0 / vg - 1.0 / C0))

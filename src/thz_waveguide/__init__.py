"""Semi-analytical circular metallic waveguide model for broadband THz beams."""

from .model import (
    C0,
    EPS0,
    MU0,
    Z0,
    Lens,
    ModeFields,
    ModeSpec,
    Waveguide,
    aperture_throughput,
    conductor_attenuation,
    coupling_amplitude,
    focused_waist,
    mode_fields,
)

__all__ = [
    "C0",
    "EPS0",
    "MU0",
    "Z0",
    "Lens",
    "ModeFields",
    "ModeSpec",
    "Waveguide",
    "aperture_throughput",
    "conductor_attenuation",
    "coupling_amplitude",
    "focused_waist",
    "mode_fields",
]

__version__ = "1.0.0"

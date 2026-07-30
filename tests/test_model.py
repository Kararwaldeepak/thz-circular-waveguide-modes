"""Numerical sanity checks for the waveguide model."""

import unittest

import numpy as np

from thz_waveguide.model import (
    Lens,
    ModeSpec,
    Waveguide,
    aperture_throughput,
    conductor_attenuation,
    focused_waist,
    mode_fields,
)
from thz_waveguide.simulation import coupling_modes, simulate_frequency


class CircularWaveguideTests(unittest.TestCase):
    def setUp(self) -> None:
        self.waveguide = Waveguide()
        self.lens = Lens()

    def test_te11_cutoff(self) -> None:
        cutoff_thz = ModeSpec("TE", 1, 1, "sin").cutoff_hz(self.waveguide) / 1e12
        self.assertAlmostEqual(cutoff_thz, 0.0878492, places=6)

    def test_tm01_cutoff(self) -> None:
        cutoff_thz = ModeSpec("TM", 0, 1).cutoff_hz(self.waveguide) / 1e12
        self.assertAlmostEqual(cutoff_thz, 0.1147425, places=6)

    def test_focused_waist_at_one_thz(self) -> None:
        waist_mm = focused_waist(1e12, self.lens) * 1e3
        self.assertAlmostEqual(waist_mm, 0.375696, places=5)

    def test_gaussian_throughput_is_bounded(self) -> None:
        throughput = aperture_throughput(0.5e-3, self.waveguide)
        self.assertGreater(throughput, 0.0)
        self.assertLessEqual(throughput, 1.0)

    def test_mode_power_and_loss_are_positive(self) -> None:
        mode = ModeSpec("TE", 1, 1, "sin")
        fields = mode_fields(mode, 0.5e12, self.waveguide, grid_size=121)
        alpha = conductor_attenuation(mode, 0.5e12, self.waveguide, fields)
        self.assertGreater(fields.power_w, 0.0)
        self.assertGreater(alpha, 0.0)
        self.assertTrue(np.isfinite(alpha))

    def test_energy_guard(self) -> None:
        rows, spectrum = simulate_frequency(
            1.0e12,
            self.waveguide,
            self.lens,
            coupling_modes(8),
            grid_size=101,
        )
        self.assertTrue(rows)
        self.assertLessEqual(
            spectrum.total_coupled_fraction,
            spectrum.aperture_throughput + 1e-12,
        )
        self.assertLessEqual(
            spectrum.total_output_power_fraction,
            spectrum.total_coupled_fraction + 1e-12,
        )


if __name__ == "__main__":
    unittest.main()

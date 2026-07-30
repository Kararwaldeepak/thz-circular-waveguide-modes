# Default simulation summary

These values were generated with the default 2.0 mm-ID, 40 cm copper guide and
the 100 mm focal-length, 50.8 mm-diameter lens. The Gaussian fills the lens
(`fill_factor = 1.0`).

## Selected cutoff frequencies

| Mode | Cutoff (THz) |
|---|---:|
| TE11 | 0.087849 |
| TM01 | 0.114743 |
| TM11 | 0.182824 |
| TE12 | 0.254382 |

## Selected broadband results

| Frequency (THz) | Focused 1/e field radius (mm) | Aperture throughput | Coupled modal power | Output power after 40 cm | Dominant output mode |
|---:|---:|---:|---:|---:|---|
| 0.20 | 1.878 | 0.433 | 0.343 | 0.289 | TE11 |
| 0.50 | 0.751 | 0.971 | 0.959 | 0.782 | TE11 |
| 1.00 | 0.376 | 1.000 | 1.000 | 0.733 | TE11 |
| 1.50 | 0.250 | 1.000 | 1.000 | 0.702 | TE12 |
| 2.00 | 0.188 | 1.000 | 1.000 | 0.681 | TE12 |

All power values are normalized to the full Gaussian beam immediately before
the waveguide aperture.

The high total coupling at higher frequencies does **not** imply single-mode
operation. The focused waist becomes much smaller than the bore, so the field
is represented by an increasing number of radial TE1n and TM1n modes. The
individual modal spectrum in `results/modal_spectrum.csv` and
`results/03_output_modal_power.png` should therefore be examined alongside the
total power.

These are idealized semi-analytical predictions. Experimental insertion loss
also includes alignment, reflection, flange, roughness, lens, and detection
effects that are outside this model.

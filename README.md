# THz Circular Copper Waveguide Modes

A reproducible semi-analytical simulation of a broadband, linearly polarized
Gaussian THz beam coupled through a lens into an air-filled circular copper
waveguide.

## Default physical system

| Parameter | Value |
|---|---:|
| Waveguide material | Copper |
| Inner diameter | 2.0 mm |
| Length | 40 cm |
| Copper conductivity | 5.8 × 10⁷ S m⁻¹ |
| Lens focal length | 100 mm |
| Lens clear diameter | 50.8 mm (2 inch) |
| Lens f-number | 1.97 |
| Lens numerical aperture | 0.246 |
| Frequency sweep | 0.10–2.00 THz |
| Input polarization | Linear, along x |

The default Gaussian 1/e **electric-field radius on the lens** is 25.4 mm, so
the beam fills the lens (`fill_factor = 1`). This matters: focal length and lens
diameter alone do not uniquely determine the focused THz spot unless the beam
size incident on the lens is known.

## What the simulation calculates

- TE and TM cutoff frequencies from the appropriate Bessel-function roots.
- Frequency-dependent propagation constant and relative group delay.
- Full vector eigenfields \(E_x,E_y,E_z,H_x,H_y,H_z\).
- Diffraction-limited Gaussian waist produced by the specified lens.
- Aperture clipping and Gaussian-to-mode power overlap.
- Copper-wall conductor loss from a surface-resistance perturbation.
- Modal power at the input and after 40 cm propagation.
- Coherent output field maps at 0.20, 0.50, 1.00, and 1.50 THz.
- CSV tables and publication-ready PNG figures.

For a perfectly centered, cylindrically symmetric, x-polarized Gaussian, only
one angular orientation of the \(m=1\) families is symmetry-accessible:
TE\(_{1n}\)-sin and TM\(_{1n}\)-cos. Other modes remain important in a real
experiment when the beam is offset, tilted, distorted, or coupled through a
non-ideal flange.

## Quick start

Python 3.10 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python run_simulation.py
```

The results are written to `results/`.
Numerical highlights from the supplied default run are collected in
[`RESULTS_SUMMARY.md`](RESULTS_SUMMARY.md).

Run the tests with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Change the experiment

All important quantities are command-line parameters:

```bash
python run_simulation.py \
  --diameter-mm 2.0 \
  --length-cm 40 \
  --focal-length-mm 100 \
  --lens-diameter-mm 50.8 \
  --fill-factor 0.80 \
  --minimum-thz 0.10 \
  --maximum-thz 2.00
```

If the measured 1/e THz field radius on the lens is \(w_L\), set

```text
fill_factor = w_L / 25.4 mm
```

For example, a measured 20 mm radius corresponds to `--fill-factor 0.7874`.

## Output files

| File | Contents |
|---|---|
| `cutoff_frequencies.csv` | Sorted TE/TM cutoff table |
| `modal_spectrum.csv` | Coupling, loss, output power, beta, and delay per mode |
| `total_spectrum.csv` | Total aperture, coupled, and output power fractions |
| `simulation_metadata.json` | Geometry, settings, and model assumptions |
| `01_mode_cutoffs.png` | Cutoffs of Gaussian-accessible modes |
| `02_coupling_and_output.png` | Aperture, coupled, and transmitted power |
| `03_output_modal_power.png` | Frequency-resolved output modal content |
| `04_copper_loss.png` | Copper loss over the 40 cm guide |
| `05_output_field_maps.png` | Coherent output fields at four THz frequencies |

## Model scope

This code is intended for fast physical interpretation and parameter sweeps. It
is **not** a full 3D finite-element model. The eigenfields assume a perfectly
circular, air-filled bore with ideal geometry. Copper loss is evaluated using
the good-conductor surface resistance

\[
R_s=\sqrt{\frac{\pi f\mu_0}{\sigma}}
\]

and first-order perturbation of the perfect-conductor fields. The model omits:

- input/output flange reflections and impedance matching;
- finite wall thickness, surface roughness, oxidation, and seams;
- waveguide bends, ellipticity, tapers, and mechanical tolerances;
- lens aberration and measured spatial phase;
- atmospheric absorption;
- offsets and angular misalignment.

Consequently, `copper_loss_db` is propagation loss from the wall model, whereas
the reduction between incident and coupled power is coupling/aperture loss.
Experimental insertion loss can be substantially larger than copper propagation
loss.

The equations and normalization are documented in
[`docs/model.md`](docs/model.md).

## Citation

If you use this project in research, cite the repository metadata in
`CITATION.cff` and describe any measured beam radius, alignment, surface
roughness, and flange geometry added to the default assumptions.

## License

MIT License.

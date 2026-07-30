# Physical model

## 1. Circular-waveguide eigenvalues

For bore radius \(a\), the transverse wavenumber is

\[
k_c=\frac{x_{mn}}{a}.
\]

For TE modes, \(x_{mn}\) is the \(n\)-th nonzero root of
\(J_m'(x)=0\). For TM modes it is the \(n\)-th root of \(J_m(x)=0\).
The cutoff frequency and propagation constant are

\[
f_c=\frac{c x_{mn}}{2\pi a},\qquad
\beta=\sqrt{k_0^2-k_c^2}.
\]

The code evaluates both cosine and sine angular orientations,
\(\cos(m\phi)\) and \(\sin(m\phi)\), when applicable.

## 2. Vector fields

The phasor convention is

\[
\exp(+j\omega t-j\beta z).
\]

For a TE mode, \(H_z=\psi\) and \(E_z=0\):

\[
\mathbf E_t=
\frac{j\omega\mu}{k_c^2}\hat{\mathbf z}\times\nabla_t\psi,\qquad
\mathbf H_t=-\frac{j\beta}{k_c^2}\nabla_t\psi.
\]

For a TM mode, \(E_z=\psi\) and \(H_z=0\):

\[
\mathbf E_t=-\frac{j\beta}{k_c^2}\nabla_t\psi,\qquad
\mathbf H_t=-\frac{j\omega\epsilon}{k_c^2}
\hat{\mathbf z}\times\nabla_t\psi.
\]

The forward modal power is calculated numerically:

\[
P_m=\frac{1}{2}\operatorname{Re}
\iint_A(\mathbf E_m\times\mathbf H_m^*)\cdot\hat{\mathbf z}\,dA.
\]

## 3. Lens and Gaussian field

If \(w_L\) is the 1/e field radius on a thin lens of focal length \(f_L\),
the diffraction-limited paraxial waist is

\[
w_0(f)=M^2\frac{\lambda f_L}{\pi w_L}.
\]

The default uses \(M^2=1\) and \(w_L=D/2=25.4\) mm. The aperture-plane
field is

\[
\mathbf E_{\mathrm{inc}}(r)
=\hat{\mathbf x}\exp\left(-\frac{r^2}{w_0^2}\right).
\]

The Gaussian power fraction intercepted by the radius-\(a\) bore is

\[
T_{\mathrm{ap}}=1-\exp\left(-\frac{2a^2}{w_0^2}\right).
\]

## 4. Modal coupling

The forward modal amplitude is estimated by power orthogonality:

\[
a_m=
\frac{\iint_A
(\mathbf E_{\mathrm{inc}}\times\mathbf H_m^*)\cdot\hat{\mathbf z}\,dA}
{2P_m}.
\]

The launched modal power is \(|a_m|^2P_m\). A small energy guard removes
finite-grid quadrature excess if the sum of projected modal power is slightly
larger than the geometrically admitted Gaussian power.

This is an aperture mode-matching estimate. A rigorous discontinuity solution
would include reflected free-space fields and evanescent waveguide modes.

## 5. Copper conductor loss

For conductivity \(\sigma\), the good-conductor surface resistance is

\[
R_s=\sqrt{\frac{\pi f\mu_0}{\sigma}}.
\]

The wall power dissipated per unit propagation length is

\[
P_{\mathrm{wall}}=
\frac{R_s}{2}\oint_{\partial A}|\mathbf H_t|^2\,ds.
\]

For field dependence \(\exp(-\alpha z)\),

\[
\alpha=\frac{P_{\mathrm{wall}}}{2P_m}.
\]

The power after guide length \(L\) is

\[
P_m(L)=P_m(0)\exp(-2\alpha L),
\]

and the corresponding power loss is

\[
\mathcal L_{\mathrm{dB}}=
\frac{20}{\ln 10}\alpha L.
\]

## 6. Dispersion

For an ideal air-filled guide,

\[
v_g=c\sqrt{1-\left(\frac{f_c}{f}\right)^2}.
\]

The reported delay is relative to free-space propagation over the same length:

\[
\Delta\tau=L\left(\frac{1}{v_g}-\frac{1}{c}\right).
\]

## 7. Coherent output field

At each selected frequency the output phasor is

\[
\mathbf E_{\mathrm{out}}(x,y;f)=
\sum_m a_m\mathbf E_m(x,y;f)
\exp[-\alpha_mL-j\beta_mL].
\]

Because distinct modes retain their propagation phase, the transverse output
pattern includes coherent multimode interference.

## Recommended interpretation

- `coupling_fraction_full_beam` includes aperture clipping and mode overlap.
- `coupling_fraction_aperture` is referenced only to power entering the bore.
- `copper_loss_db` is wall propagation loss, not total experimental insertion
  loss.
- `output_power_fraction` is referenced to the full Gaussian before the
  waveguide aperture.

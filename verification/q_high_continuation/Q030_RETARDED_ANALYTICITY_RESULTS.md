# q=0.30 retarded-response validation and analyticity barrier

**2026-11-05. Source: user-provided exact finite-80-shell Einstein–Vlasov
background plus newly generated q=.30 action-angle spectral nodes.
This is NOT a global q=.30 QNM or a gravitational Green-function residue.**

Physical input from the frozen reproducibility bundle:
`Einstein_Vlasov_q_continuation_20261105.zip`,
`validated/q_0.30/spectral_nodes.npz`,
SHA256 `f325b1be21e8a4698c2b1a880a02f67c1991413503529365caa1ff9e760cb75a`.

- Exact original massless distribution: `F=C(E0-E)_+(L-L0)_+`.
- 80 independently glued finite-width Vlasov shells.
- New gravitational mass ratio `q=0.30000000000316`.
- 87,840 radial nodes in 1,440 branches and 4,880
  carrier nodes in 80 branches.
- Original 61-vs-81 orbital quadrature upper-plane C7
  source-response discrepancy `1.27699%`, so input is
  acceptable for *initial susceptibility diagnostics*.

## An independently evaluated genuinely new q=.30 retarded susceptibility

At the **previously published q=.20 matter-mode frequency**

```
M omega_seed = 0.3797013686231943 - 0.0239854837986382 i
```

a separate reproduction of the original
`hbh_phase_xxic8_landau_continuation.py`
piecewise-linear retarded prescription gives:

```
chi_J(q=.30, omega_seed) = 0.0568142175282 + 1.5361612191563 i
chi_X(q=.30, omega_seed) = 1.9446338658585 + 1.1304170547293 i
```

Algebraically stabilized midpoint integration of the
SAME piecewise-linear C8 spectral integrals gives:

```
chi_J = 0.0568153514152 + 1.5361705409685 i
chi_X = 1.9446462831814 + 1.1304229382401 i
```

The relative differences in values are only ~6.11e-6,
caused by cancellation in the original sum of logarithms
and linear terms. The original C8 spectrum `q=.10`
is independently reproduced with
`chi_J = 0.0904420457063 + 1.484044177180 i`
at the published q=.10 root, near the frozen production
`0.0884816309420 + 1.484577616806 i`
(0.137% relative difference).

**These chi numbers are NOT a q=.30 complex eigenfrequency.
The old q=.20 omega is a fixed evaluation point / root seed.**

## Critical issue: C8's lower-sheet continuation is not holomorphic across knots

C8 first approximates the real-frequency density on each
branch by separate linear segments. It then applies the
second-sheet subtraction
`-2*pi*i*rho_segment(z)`, where the chosen segment is
selected with `searchsorted(Re(z))`. At a real
spectral knot `x_k` the adjacent linear densities
agree only at real `z=x_k`:
`rho_left(x_k)=rho_right(x_k)`. For complex z
with `Re(z)=x_k` and `Im(z)!=0`, their analytic
linear extrapolations generally differ by
`i*Im(z)*(slope_right-slope_left)`.
Therefore the prescribed lower-sheet quantity
has **complex-valued jumps along vertical knot
lines** and cannot be treated as a holomorphic
function on neighborhoods crossing those lines.

This is NOT just a floating-point concern; it is
a property of the chosen numerical continuation.

### Measured q=.30 failure

At `Re(M omega_seed)=0.379701368623`:

- **38 real spectral knots** fall inside ±2e-5,
  **18 knots** inside ±1e-5.
- Nearest knot separation ~9.03e-7.
- At the nearest knot the continued susceptibility
  jumps by roughly
  `Delta chi_J=-1.85e-4`,
  `Delta chi_X=+2.25e-4 i`.
- At finite-difference real step `h=1e-5`,
  dchi_J/domega differs from the within-cell
  analytic/imaginary-direction derivative by
  **~19.6× its magnitude**, with dchi_X/domega
  discrepancy **~18.3×**.
- After applying a stable midpoint formula, the
  imaginary-direction finite-difference derivative
  matches the exact within-cell analytic formula
  to `<1e-8` relative at h=1e-5.
- At `h=2e-7`, even real-direction finite
  differencing agrees to `<1e-9`, but such
  an extremely narrow cell is not a coherent
  physical analytic domain for tracking a
  QNM across changes in Re(omega).

### Crucial original-case control

For the archived reconstructed q=.10 spectrum at
its verified published q=.10 root, there are
**zero spectral knots within ±2e-5**. The
real/imag derivative comparisons agree.
The newly added CI unit test using the public
official q=.10 80-shell source reports
`CR_relative_disagreement=1.7973e-8`
and `PASS`.

Thus the q>.20 problem is NOT evidence that
the already verified q=.10 calculation used
the same unsafe local derivative regime.

## Scientific decision and enforced code behavior

The existing published C9 source-residue method
numerically differences an unnormalized Wronskian
at h=2e-5 and h=1e-5. The *new* q=.30
retarded susceptibility cannot safely support
that derivative or a global analytic
argument-principle root count in its present
piecewise-Re(z) interpolation.

The downstream `run_new_q.py --solve` now
**fails closed and writes a blocked JSON** if
the analyticity preflight finds this
frequency-direction discrepancy, rather than
publishing a misleading q=.30 pole or
source-normalized gravitational residue.

**What fixes it:** construct and validate a
genuinely holomorphic retarded continuation of
the *underlying kinetic action-angle integral*
on the frequency domain used for the global QNM.
This may require an analytic continuation of
the orbital frequency/vertex maps or an
independently certified spectral-density
representation; merely reducing finite-difference
step h is insufficient to follow a mode across
many spectral knots.

Before publishing a new mode, require:
1. Uniform complex-frequency convergence of
   the kinetic susceptibility under orbital
   quadrature refinement and choice of analytic
   representation.
2. A complex-analytic derivative agreeing along
   both real and imaginary directions on a
   finite neighborhood of the root.
3. Independent root continuation/winding
   checks with the same retarded prescription.
4. Normalized and unnormalized W derivatives,
   identically normalized exterior point and
   smooth source residues with stable q-grid
   and radial boundary/resolution controls.

No gravitational-wave detector strain was
examined in this study. This issue must be
resolved before treating the new q=.30
spectrum as a physically validated matter
QNM or ringdown detectability forecast.

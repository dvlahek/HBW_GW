# q=0.25 — controlled 80-shell Einstein–Vlasov QNM and exterior-source residue **candidate**

**Date:** 2026-11-05. **Evidence:** independent local C1–C7 physical
construction and GitHub C9 numerical execution.
**CI:** https://github.com/dvlahek/HBW_GW/actions/runs/42085159420
**JSON artifacts:**
https://github.com/dvlahek/HBW_GW/actions/runs/42085159420/artifacts/13443000996

This is a **reproducible new numerical candidate in the effective C9 model**,
NOT yet a physically certified QNM of the exact 80-shell background and
NOT an observation of matter in detector data.

## Physical source reconstruction

From the original massless Vlasov shell distribution
`F=C(E0-E)_+(L-L0)_+`, the uploaded original C1/C2/C4/C6/C7
scripts and the same previously frozen **scaled-delta path** used
for q=.30, recompute an **N=80 exact stationary shell staircase**
at q_target=.25:

| Metric or source check | q=.25 result |
|---|---:|
| Constructed inner-to-outer gravitational mass ratio | 0.2500000000009104 |
| Per-shell mass ratio | 0.9828205985452511 |
| Finite-shell regulator delta | 0.00005550417320229213 |
| Mass equation relative L2 residual | 3.90803125145e-7 |
| Lapse equation relative L2 residual | 3.06014291093e-9 |
| Max massless stress-energy trace residual | 5.68301807837e-16 |
| Minimum vacuum gap between physical shell supports | 0.00324452090493 |

These certify a **finite N=80 stationary background** within the
specified numerical tolerances. This does not prove physical
stability, N→∞ convergence, or astrophysical plausibility.
The pre-existing published q=.20 root used as a *seed*
was computed on a different N=56/regulator path. Hence
branch identity is not fully established solely by a small
frequency shift.

At the upper-half-plane C7 reference frequency,
independently regenerated real-frequency spectral nodes:

| nxi | radial nodes | carrier nodes | physical NPZ SHA256 |
|---:|---:|---:|---|
| 81 | 116,640 | 6,480 | ba232d44e02f474abbd6c4ba353d4e22894d7322b726514a8d5c17307d9bfdf3 |
| 101 | 145,440 | 8,080 | b3ffd96191ea6387f64d2b6c7da28efbdbd07dc00b0b3e03198bd9a795d3eb96 |
| 641 | 923,040 | 51,280 | e1bdfb360f75343c56df9a559c4b28eae11c4e472b5b010bbb8d32994fde85f3 |
| 801 | 1,153,440 | 64,080 | 998804e476d7208d315f94337298068c0f3476fe2bd5d6a70dd2654a97c87298 |

Between nxi=641 and nxi=801, upper-plane
chiJ shifts 0.0553% and chiX shifts 0.0213% at
the fixed reference point, not a proof of uniform
convergence in complex frequency.

## Locked analytic local C8 second-sheet approximations

From exactly the NEW high-resolution q=.25 input arrays
at the previously published q=.20 pole
`omega0=.3797013686231943-.0239854837986382i`,
construct three complex Taylor-12 sheets with frozen
crossing anchor:
`n641_a0, n801_a0, n801_am0005`.
Their common verified local validity disk radius is
`|omega-omega0|<=.004`. The original C8 Cauchy
integrals are retained; only the second-sheet
discrete spectral-branch assignment is frozen.
The normalized heldout Taylor-vs-exact-discrete-sheet
error is <5e-16.

GitHub source `q025_taylor_regenerated.py` stores
compressed payload SHA256
`58a0fbb91aa32fdb8263e362c3884e0d5c8353b5b97ef3949be8b6da2e4a103a`
and decompressed JSON SHA256
`26d15ca0261491415041d35612e720446d9135605f0cfd92cc4a6b0b24965441`.

At omega0:
- nxi641 vs 801: chiJ value 0.2633%, derivative 1.2387%;
  chiX value 0.04365%, derivative 0.11085%.
- nxi801 anchor at -0.0005 vs nominal:
  chiJ value 0.6167%, derivative 3.644%;
  chiX value 0.1831%, derivative 1.2503%.

The fixed discrete sheets are analytic locally.
**Uniqueness of the actual continuum retarded
Einstein–Vlasov sheet is NOT established.**

## New bounded global effective C9 mode candidate

The original 7×4 conservation-complete **effective-ocean**
model, using C9 Schwarzschild ECS boundary matching,
is evaluated with unchanged q=.20 mode and lambda seeds.
Each real/imag search displacement is confined to
[-.0024,+.0024] around omega0; acceptance
requires interior root and normalized |D|<2e-6.

| Source | Local numerical M omega | normalized |D| |
|---|---|---:|
| n641_a0 | 0.38099074548973594 - 0.025873483312466364i | 4.43531e-10 |
| n801_a0 | 0.3810200143912758 - 0.025862906485263668i | 2.96219e-10 |
| n801_am0005 | 0.3809538387456802 - 0.025896337962632018i | 9.94221e-10 |

All three are interior roots satisfying the
predeclared *numerical* cut, with maximum
pairwise complex spread from nxi801 reference
`7.41409451587e-5`, under the frozen
`2e-4` cross-configuration tolerance.

Two circular net-argument-winding checks
`r=.003,N=32` and `r=.0035,N=48`
give **+1 in all six configurations**.
This is a numerical *zeros minus possible
poles* net winding of the Riccati/ECS
representation. It is not an independent
pole-free analytic full-state certificate.

## NEW exterior source Green-function residue audit

Using the original fixed source families
(six exterior point sources and three
normalized smooth sources), the original
C9 effective-ocean Green-function
construction was evaluated at each of the
three local numerical candidates.

| Quantity | nxi641 | nxi801 (reference) | nxi801 shifted anchor |
|---|---:|---:|---:|
| Two-step dW/domega relative difference | 4.23554e-7 | 4.23081e-7 | 4.17675e-7 |
| Point-source residue modulus range | 0.105070–0.117924 | 0.105156–0.118022 | 0.105092–0.117941 |
| Smooth-source residue modulus range | 0.069710–0.223951 | 0.069778–0.224124 | 0.069710–0.224021 |

Maximum relative changes in *residue modulus*,
against the nxi801 nominal case, across all
point/smooth source locations:
- six exterior point sources: **0.000836675**
  = **0.0836675%**;
- three exterior smooth sources: **0.000976127**
  = **0.0976127%**.

These confirm **numerical residue stability in
the *effective C9* model**, not the physical
radiative metric/strain amplitude of the exact
finite-shell Einstein–Vlasov model. They are
not vacuum-normalized excitation ratios or
predicted merger mode amplitudes.

## Scientific boundaries and next step

1. The C9 effective-ocean radial response
   has not been derived as the exact linearized
   Einstein–Vlasov dynamics of this new N=80
   multi-shell stationary background.
2. Matching analytic discrete C8 sheets
   is necessary but insufficient to uniquely
   define the causal continuum retarded sheet.
3. The old q=.20 seed comes from a
   separately tuned 56-shell background.
   Establish a common N/regulator q-grid
   q=.20,.225,.25,.275,.30 and compare
   root continuation and Green source couplings.
4. The **earlier q=.30 negative search**
   was centered on the OLD q=.20 frequency
   and sampled only disks <=.0035.
   The q=.25 frequency shift implies the
   q=.30 branch might have traveled beyond
   that old region. Therefore the q=.30
   result does **not** imply the pole
   disappeared. Recenter q=.30 around
   q=.25 as an explicitly exploratory
   new predeclared study, with actual
   q=.30 641/801 C7 source regeneration.
5. Independently verify a full-state
   Wronskian root and physical exterior
   metric/strain response before claiming
   astrophysical matter-led ringdown detectability.

No GW detector strain is used, no
real-world event matter mode is claimed,
and the previously submitted CQG work
is not modified by this independent project.

# Recovered original C1–C7 pipeline: new Einstein–Vlasov q>0.20 backgrounds

**Research audit, 2026-11-05.** The original user-local files
`hbh_phase_xix_exact_vlasov_shell.py`,
`hbh_phase_xxic2_action_angle_audit.py`,
`hbh_phase_xxic4_smooth_polarization.py`,
`hbh_phase_xxic6_conservation_complete_response_v3.py`,
`hbh_phase_xxic7_spectral_response.py`,
`phase_xxic1_output.zip`, and `phase_xxic6_output.zip` were recovered.
The uploaded code and numerical artifacts have **not** been committed
into this public PR; their original bytes and resulting new-q spectra
are in a separate user-delivered reproducibility ZIP.

## New physics source: the C1 massless Einstein–Vlasov shell ODE

The original exact ODE uses
`F=C (E0-E)_+ (L-L0)_+` with `k=ell=1`, solves the
static `m'(r)=4pi r² rho` and lapse equation on the
compact support of each shell, and stops at its actual
outer matter-support zero. Each original C1 stage consists
of 80 self-similar shells at `q_total=.1`, each having
`q_step=0.9716279515771061`; the original shell CSV
contains masses, radial edges, lapse junctions, E0/L0
and phase-space normalizations.

Using that explicit recurrence and independently RETUNING
the true compact-support ODE's normalization C at each new
target `q_total`, we can reconstruct a genuine exact
**finite N=80** equilibrium for new q. This is **not**
an extrapolation or rename of the q=.1 spectral-node NPZ.
Each new q has its own regenerated full 87,840 radial
and 4,880 carrier C7 action–angle spectral nodes,
with an independently computed source SHA256.

One continuous exploratory physical-regularization path:
`delta(q)=[.005(14/80)^2] * [ln(1/q)/ln 10]^2`.
The original fixed delta at q=.1 is preserved. This newly
specified q-parameter path is distinct from automatically
assuming q alone uniquely fixes a collisionless atmosphere.
All future q-continuation claims must test its regulator
and N-convergence dependence.

Independent Einstein field residuals were computed by
differentiating high-resolution reconstructed metric
profiles and comparing with moments of the ORIGINAL
massless phase-space distribution at each radius.
For the new q=.3/.5/.7/.9/.95 backgrounds,
relative L2 mass-constraint errors are around
3.0–3.9e-7 and max massless tracelessness errors
around 5e-16. Independent radial hydrostatic
stress-energy residual at q=.3/.5/.9 is around
(4–7)e-11 in the original template units.
Vacuum shell junction mass/lapse continuity and
non-overlap pass. These are verified finite-N
backgrounds, NOT proven N→infinity limits.

## The critical recovered numerical bug: original C4 weights

Original `c4.smooth_weights` subtracts nearly equal
powers `E0^n-Emin^n`. Its effective `A4` and
`D5` energy moments vanish like
`t²`, `t=1-Emin/E0`, while its normalized energy
moment vanishes like `t³`. Naive differences at
q=.9/.95 suffer catastrophic cancellation.

For stable C7 builds only, the user-local source
recovery uses mathematically IDENTICAL factored
polynomials:

```
x = Emin/E0, t = 1-x
A4 = B² E0^7 t² (2x^5+4x^4+6x^3+8x²+10x+5)/35
D5 = E0^7 t² (6x^5+5x^4+4x³+3x²+2x+1)/42
Ered = (2BC/Omega_r) B² E0^7 t³
       (24x^4+37x³+39x²+30x+10)/420
```

The original C7 N80 q=.1 spectral reproduction passes
C6 independent upper-half-plane checks:
relative `chiJ` error **1.0784%** and `chiX`
error **1.2567%**.

**Do not interpret apparent suppression before this
algebraic fix.** At q=.9 the old unstable
`|chiJ|≈0.01238` becomes `|chiJ|≈0.64472`
with the factored weights at 61 orbital samples.
At q=.95 the old apparent `|chiJ|≈1.35e-5`
becomes `|chiJ|≈0.402`. These high-q
values still have inadequate independent orbital
quadrature convergence, so **no claim of
kinetic response survival or extinction is made**.

## Finite-q source evidence and convergence

All these are values of `chiJ, chiX` at
the old C2 **upper-half-plane** reference frequency
`0.370658573899+0.00590424942625i`,
**NOT** newly computed QNM poles or gravitational
residues.

| New finite-N total q | |chiJ| (nxi=61) | |chiX| | Max relative 61 vs 81 orbital quadrature error | Gate |
|---|---:|---:|---:|---:|---|
| 0.30 | 1.13074 | 1.67204 | **1.28%** | preliminary C7 spectral accuracy PASS |
| 0.50 | 1.12977 | 1.66642 | **1.30%** | preliminary C7 spectral accuracy PASS |
| 0.70 | 0.88857 | 1.20527 | 17.96% | FAIL — resample |
| 0.90 | 0.64472 | 0.87799 | 26.42% | FAIL — resample |
| 0.95 | 0.40203 | 0.38521 | **not fully evaluated** | FAIL — incomplete |

A separate physical q=0.298538 42-shell truncation
of the old q=.1 C1 chain gives
`chiJ=0.256275+1.095799i`,
near the independently constructed 80-shell
scaled-delta q=.3 case
`chiJ=0.264125+1.099464i`.
This is a useful cross-family consistency
control near q=.3, not a systematic N-limit.

**New q=.30 and q=.50 NPZ plus physical input
manifests exist in the user-delivered local ZIP.**
q=.70/.90/.95 full physically solved
backgrounds exist but their *spectral responses*
are not accepted as converged. This PR's validator
now requires an explicit independent C7
61-vs-81 quadrature check at max 3% and hard-stops
failed higher-q inputs before running the global
complex root search.

## What remains unsolved

The project now has an ACTUAL upstream finite-q
Einstein–Vlasov background and node generator,
so the earlier blanket statement that all q>.20
spectral inputs are unavailable is superseded
for **q=.3 and q=.5**.

However, **no new q>.20 global complex QNM pole,
mode residue or physical gravitational waveform
amplitude has yet been calculated** in this
work. The new q=.3 source manifold must enter
the retarded Landau response, global 7x4
Einstein–Vlasov ECS matching, and exterior
source Green-function derivative/residue,
with a NO-ROOT vs ROOT classification and
independent boundary/resolution audits. q=.5
cannot be jumped to directly from q=.3
under the current max 0.11 q-step root-identity
rule; build intermediate q values first.

No real gravitational-wave event or fresh H1/L1
detector noise was analyzed, and no observational
matter-led mode detection is claimed. This
remains a new research line separate from the
previously submitted manuscript.

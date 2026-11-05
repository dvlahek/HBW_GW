# q → 1: exact vacuum endpoint, measured residue audit and the missing continuation

**2026-11-05. Reproducible endpoint audit, NOT a numerical q=.3..999 Einstein–Vlasov pole continuation.**

Our task is to determine whether the long-lived Einstein–Vlasov
`ell=2` polar matter-led mode remains gravitationally observable
for *dilute environments*. The source mode is physically
self-consistent only on the nonrotating smooth Einstein–Vlasov
background, which is constructed at q=.05, .075, .10, .15, .20
in the present publication archive.

## Already measured finite-q evidence

The mass ratio `q=m_inner/M_outer` (NOT the merger binary mass
ratio) is the total gravitational mass ratio across the atmosphere.
For q=.05–.2 the gravitational mass difference is **80–95%
of the outer mass**. It is *not* a direct proper-rest-mass
measurement. Those rows are strongly self-gravitating and
are far from the q→1 dilute limit.

The archive contains exactly five self-consistent smooth-background
poles at q=.05/.075/.10/.15/.20, but a re-runnable complete
spectral-node catalogue is deposited only at q=.10.
The complete q-family archive and newly constructed q>.20
spectral-node inputs are absent. They CANNOT be synthesized
from interpolation of the known root table. Original
production code validates q=.10 and fails closed if a
finite-q case is requested without its real spectral
catalogue.

Same-exterior-source Green-function **residue amplitude ratios**:

| Source family | q=.10 | q=.20 |
|---|---:|---:|
| Six point sources | 0.96%–2.50% | 1.60%–3.98% |
| Three smooth sources | 1.17%–1.79% | 1.93%–2.89% |

These are **pole-residue ratios**, not predicted merger
waveform amplitudes or detector signals. The finite-duration
Gaussian-cosine burst source ratios are likewise
conditional on a *prescribed source profile*, not a
merger-excitation distribution. At q=.10 and width 5M,
the burst-end ratio ranges ~4.16%–6.39%;
q=.20 width 5M ~6.44%–9.75%.

**NO EXTRAPOLATION** of these ratios into q>.20,
especially not q=.90 or .99. The available trend
between .10 and .20 does not establish weak-environment
survival or observable excitation.

## A real, computable endpoint

In the original two-region global transfer construction
the dimensionless ocean coordinate interval is

```
L/M = 9 ln(1/q).
```

For example:
- q=.10 ⇒ L/M≈20.723
- q=.20 ⇒ L/M≈14.485
- q=.90 ⇒ L/M≈0.948
- q=.99 ⇒ L/M≈0.09045
- q=1 ⇒ **L=0 and exactly T=I**.

The original source code explicitly returns an exact
2×2 identity transfer matrix at q=1 without needing
a kinetic spectral catalogue. We can therefore compute
the **correct q=1 Schwarzschild vacuum Zerilli determinant**
using the existing horizon and exterior-complex-scaled
outgoing solutions. The audit tests the argument
winding around the original q=.10 matter root
`Mω=0.377514814460-0.020871654486i`
using a frozen radius .01 circle and 32 points, and
compares it with the known vacuum gravitational
fundamental
`Mω=0.373671684352-0.088962315642i`.

If the q=1 determinant winds **zero** on the q=.1
matter-pole disk, then that fixed complex region
contains no q=1 vacuum root (within the numerical
contour-method validity). This does NOT compute
what happens to the finite-q branch at q=.9.

The winding count is a numerical diagnostic, not
an interval-arithmetic theorem. We can change
contour sampling/Rmax as independent checks,
but never tune the contour or center to
force a desired zero count.

## Conditional mathematics of the q→1 limit

Suppose the exact *recomputed* finite-q transfer
`T_q(ω)` tends uniformly to `I` as q→1 on a
compact complex-frequency region with analytic
ingoing/outgoing bases, no branch cuts and
no Schwarzschild QNM zeros on the boundary.
Then its global determinant converges uniformly
to the vacuum Zerilli determinant. By Rouché's
theorem, the number of roots in that region
equals the vacuum count for sufficiently large q.

Thus a separately identifiable matter-led root
**cannot persist with a fixed nonzero gravitational
residue inside any vacuum-root-free compact
region** under those hypotheses. It may migrate
towards a vacuum QNM, the boundary or a branch cut,
or require singular/nonuniform susceptibility.
This is a **conditional implication**,
not proof that uniform convergence holds for
the actual cold kinetic limit.

One may NOT infer from `L→0` alone that
`T_q-I = O(L)`: the effective kinetic
generator could be unbounded near the
cold/dilute endpoint, precisely where a
memory-like susceptibility can be singular.
This needs an independent uniform-bound test.

## The decisive missing calculation

**To complete the user-requested actual physical q→1 continuation:**

1. Rebuild the background for q=.3/.5/.7/.9/.95/.99
   (and ideally .995/.999) with the same phase-space
   profile and stress-energy constraints.
2. Generate new action-angle radial and carrier
   spectral-node arrays for *each*, validating
   normalization and retarded continuation.
3. Track the global matter branch in complex
   ω, rechecking winding, Rmax, resolution and
   near-degenerate/cut behavior at each q.
4. Recompute `dW/dω` and pole residues against
   **identically normalized exterior point and
   smooth sources**, plus a finite-duration
   burst response, at every new q.
5. Compare source-normalized gravitational
   residue and excitation against a verified
   vacuum q=1 limit, not merely the pole
   frequency or late-time relative amplitude.
6. If the atmosphere becomes inconsistent,
   lacks stationary solutions or the kinetic
   response becomes singular, record that as a
   physical boundary of the model rather than
   forcing continuation.

**Decision:** No new astrophysical detector
forecast based on q>.20 should be claimed until
the **actual background + spectrum + source coupling**
exist. The existing Phase10 forecast is an honest
conditional phenomenological injectability result,
not a physical prediction of detectable
collisionless matter near a realistic merger.

This audit adds no new H1/L1 data, introduces no
new "detected" mode and does not modify
any submitted manuscript results.

# q=0.30 high-resolution C1–C7 / C8 / C9 global-pole test — bounded negative result

**Frozen 2026-11-05.** This is a separate matter-led ringdown
detectability research result, NOT a submitted-manuscript change.
CI: https://github.com/dvlahek/HBW_GW/actions/runs/42081748524
Raw CI JSON artifacts:
https://github.com/dvlahek/HBW_GW/actions/runs/42081748524/artifacts/13441620994

## Verified inputs

The original, physical, independently regenerated 80-shell
Einstein–Vlasov background at `q=.30000000000316` was processed
through original C1–C7 source machinery using stabilized,
algebraically equivalent C4 energy integrals.
Two independently regenerated higher orbital resolutions:

| nxi | radial spectral nodes | carrier nodes | NPZ SHA256 |
|---:|---:|---:|---|
| 641 | 923,040 | 51,280 | `92e63973d71474b4918fc2c2959227569324f69df7d98addad1838b39e730239` |
| 801 | 1,153,440 | 64,080 | `75a8a9a0a4a0ebd433343ba77a75343a3012ad9a7bb5ba1d1c4fd438222b2da0` |

The original initial compressed Taylor blob previously committed
on this branch was TRUNCATED. GitHub CI explicitly found
`zlib.error: incomplete or truncated stream`. It was
**deleted**, not reused or silently marked valid. A NEW
independently reconstructed Taylor source was generated
from the two just-regenerated physical NPZ files:
`q030_taylor_regenerated.py`, compressed payload SHA256
`c77f689c1f4606f065154175e1879ff6c2397ccb9398586e79d390d4b8120d49`.

Each dataset contains 12th-order scaled complex Taylor
coefficients for `chiJ`, `chiX` of the SAME original
C8 piecewise-linear Cauchy integral with a FROZEN
second-sheet spectral crossing anchor. The three cases:
`n641_a0`, `n801_a0`, `n801_am0005`.
They use the old published q=.20 seed
`M omega0 = .3797013686231943 - .0239854837986382 i`
and strict complex local validity radius **.004**.
The new payload passes a cryptographic checksum
and independently evaluated complex heldout interpolation
test to ~4e-16. This proves a surrogate reproduces
its CHOSEN discrete C8 sheet; it does NOT uniquely
define the physical retarded continuum sheet.

### Measured q=.30 susceptibility stability

At the same old q=.20 complex seed:

| Measurement | chiJ | chiX |
|---|---:|---:|
| nxi641 vs nxi801 relative VALUE difference | 0.5457% | 0.4319% |
| nxi641 vs nxi801 relative DERIVATIVE difference | 2.9146% | 3.0088% |
| nxi801 a0 vs a-0.0005 VALUE difference | 0.3708% | 0.2361% |
| nxi801 a0 vs a-0.0005 DERIVATIVE difference | 2.2281% | 1.6655% |

This is a bounded numerical stability observation, not
a uniform nxi→infinity or physical-sheet uniqueness theorem.

## Frozen bounded 7x4 global C9 / ECS determinant search

Use original C9 symbolic 7×4 conservation-complete
effective-ocean bulk, original vacuum Schwarzschild
ingoing/outgoing ECS boundary solver, `q=.3`,
`tau=1,edge_t=1,Rmax=80`, and published q=.20
lambda seed pair. All three new susceptibilities
use the exact same initialization and geometry.

A single constrained trust-region solve per case:
center `omega0`, each real and imaginary displacement
in `[-.0024,+.0024]`, so every trial remains
inside the Taylor validity disk. A valid root must
have normalized `|D| < 2e-6`, interior to the
box and converged. This is not a broad complex
plane QNM survey.

| Spectral source | Closest optimizer location | normalized |D| | Outcome |
|---|---|---:|---|
| n641_a0 | 0.382101368528 - 0.026385483799i | **0.01470305** | fails root acceptance |
| n801_a0 | 0.382101368602 - 0.026385483799i | **0.01459332** | fails |
| n801_am0005 | 0.382101368564 - 0.026385483799i | **0.01465042** | fails |

**All optimization endpoints are pinned to the search
box boundary.** They MUST NOT be presented as QNM
frequencies. The residual is ~7,300 times the
required numerical acceptance cutoff, and the
optimizer found NO isolated interior zero.

### Independent argument-phase diagnostic in the same disk

To check for possible missed roots, the same C9
determinant was evaluated around frozen concentric
complex contours centered on old q=.20 pole,
with `r=.003,N=32` and `r=.0035,N=48`.

For ALL THREE spectral sources, BOTH contours
give numerically **net winding 0** with well
sampled maximum phase increments <.295 rad.
Minimum normalized determinant magnitude on the
larger radius was:
- n641_a0: **.01337058**
- n801_a0: **.01346099**
- n801_am0005: **.01338297**

This is a stable **NET argument winding
(zeros MINUS possible poles)** of a
numerical C9 Riccati/ECS representation, not
a rigorous analytic pole-free disk certificate.
The smoothly chosen lambda solution branch and
an analytic full-state Wronskian still need
independent verification.

## Interpretation for detectable matter ringdown

**The old q=.20 matter mode is not recovered as a
numerically certified isolated q=.30 root within
the predeclared local complex neighborhood.**
This DOES NOT show that it vanishes anywhere
else in the complex plane. In particular, the
branch may migrate farther than the local radius,
hybridize with other modes, move to a different
Riemann sheet, or lose coupling, and the finite-N
new background follows a predeclared scaled-delta
family different from old q-family detail.

Crucial additional model limit: C9's
conservation-complete 7×4 *effective ocean*
operator remains the original idealized smooth
transfer model. The new C1 exact 80-shell
background and action-angle spectral nodes do
NOT by themselves prove that C9 is the exact
linearized radial Einstein–Vlasov operator
on this multi-shell geometry. No source-normalized
Green-function residue or detector waveform was
computed, and there is NO observational
matter-mode detection claim.

**Next step should be a physically controlled
q-grid continuation, not arbitrarily expanding
the root-search box until a root appears.**
Construct q=.25 from the same exact C1–C7
shell family, with 641/801 and frozen C8
sheet checks; compare old q=.20 mode against
q=.25. Independently derive the physical
retarded action-angle sheet and, if a pole
persists, its metric/strain residue at infinity.
Only then consider a new H1/L1 detectability
injection forecast.

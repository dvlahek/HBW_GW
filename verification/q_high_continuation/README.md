# High-q self-consistent Einstein–Vlasov continuation: executable intake

**Scope (2026-11-05):** a working, strict numerical *downstream*
continuation runner for newly generated spectral inputs. This is
**NOT** a reconstruction of the missing upstream stationary
Einstein–Vlasov background and action–angle generator, and does
**not** claim any new poles for q>.20.

## What is available

The publication repository contains:

- exact published q=.05/.075/.1/.15/.2 global roots;
- public, SHA256-verified full q=.1 spectral nodes, with expected
  1440 retarded radial branches and 80 carrier shells;
- published q=.1/.2 source-normalized point/smooth residues;
- production smooth kinetic response `code/hbh_phase_xxic8_landau_continuation.py`,
  global 7×4 Einstein–Vlasov bulk and ECS boundary matching
  `code/hbh_phase_xxic9_final_global_smooth_qnm_v4.py`;
- physical exterior Green-function residues from
  `code/hbh_np_source_excitation_test_v2.py`;
- a separately validated q=1 vacuum endpoint, with an analytic
  full-state Wronskian contour audit (draft PR #3).

It does **not** contain the upstream raw orbital/action–angle
spectral generator that can rebuild fresh equilibrium
`q>.20` arrays. The old complete q-family zip has
published results only up to q=.20; none of its caches may
be renamed, interpolated, copied or rescaled as a new q
background.

## First required upstream input

At each new q, starting at **q=.30** and stepping slowly
toward 1, a genuinely NEW smooth equilibrium calculation
must export in one directory:

- `new_q_background_audit.json`: numeric `q_integrated_from_background`,
  `shell_count`, `ocean_length_M=9*ln(1/q)`,
  `relative_mass_balance_error`, `einstein_constraint_max`,
  `vlasov_stationarity_max`, `stress_energy_conservation_max`,
  plus booleans `positive_phase_space_density`,
  `physical_orbit_support`,
  `stress_energy_from_same_distribution`.
  These diagnostics must be computed *by* the upstream
  physical construction, **not manually asserted**.
  The runner can check consistency/tolerances but cannot
  independently reproduce the Einstein–Vlasov equations
  from these scalar numbers.
- `new_q_spectral_nodes.npz`: nine original
  physical kinetic action–angle fields
  `nu_J,q_JJ,q_JX,n_r,n_phi,shell_J,nu_H,q_HH,shell_H`
  with matching shell counts and harmonic truncation.
  Generated from this SAME new equilibrium.
- `manifest.json` with q, physical-family identity,
  upstream generator script SHA256, actual producer
  command, and SHA256 of the *two real generated files*.
  This proves immutability of the supplied bytes but
  does not prove that upstream physics is correct.

Example MANIFEST FIELD STRUCTURE ONLY (not valid real
data, and do not copy placeholder checksums):

```json
{
  "schema_version": 1,
  "q": 0.30,
  "model_family": "published_reference_smooth_kinetic_family",
  "new_self_consistent_equilibrium_constructed": true,
  "q010_spectral_reused": false,
  "generator": {
    "script_sha256": "<actual 64-hex SHA256>",
    "source_description": "full source name and the physical model version actually run",
    "command": "python YOUR_ORIGINAL_UPSTREAM_GENERATOR.py --q 0.30 ..."
  },
  "background": {
    "file": "new_q_background_audit.json",
    "sha256": "<actual SHA256 of file bytes>"
  },
  "spectral_nodes": {
    "file": "new_q_spectral_nodes.npz",
    "sha256": "<actual SHA256 of file bytes>"
  }
}
```

**The full original q0.10 spectral NPZ cannot substitute for q0.30**:
its SHA is explicitly blacklisted when q>.20 is requested.

## Validation and execution

From the public HBW_GW repository root:

```bash
python -m pip install -r requirements.txt
python verification/q_high_continuation/run_new_q.py --selftest-q010

# Only when ORIGINAL new-q physical outputs exist:
python verification/q_high_continuation/run_new_q.py \
    --manifest /local/path/q030/manifest.json --validate-only

python verification/q_high_continuation/run_new_q.py \
    --manifest /local/path/q030/manifest.json \
    --solve --output-dir results/new_q030

# Only after q030 has all numerical gates PASS:
python verification/q_high_continuation/run_new_q.py \
    --manifest /local/path/q040/manifest.json \
    --previous-result results/new_q030/q_0.300_numerical_result.json \
    --solve --output-dir results/new_q040
```

The run automatically selects the previously published q=.20
root and pair of kinetic radial eigenvalue seeds for q=.30.
For every subsequent q, the immediately previous independently
verified result is mandatory. Roots and residues are NOT
interpolated; a previous published root is only an initial
guess for the actual *new* retarded determinant solve.

On validated genuine new-q inputs, the code evaluates
the original smooth constitutive response, global
complex determinant, isolated-root winding,
finite-difference derivative convergence, point
and smooth exterior Green-function residues. It
refuses missing inputs, failed background checks,
bad NPZ schema, wrong SHA, q-step jumps,
large residuals, non-isolated roots and derivative
or source-coupling instability.

## What still needs independent physical validation

These gates do not replace checking the upstream
stationary Einstein–Vlasov equations and conservation
from raw phase-space samples. Strong final publication
claims also require resolution convergence, a separate
full-state boundary-value pole audit, duplicate
finite-q source-normalized vacuum comparison and
analytic/cut controls.

**No new q>.20 pole or source residue is currently
verified.** The CI only tests the q=.10 published
input and that missing/reused new-q input fails
closed. It does NOT run a fake high-q calculation
or claim the original physical generator exists.

## Scientific acceptance criterion

The eventual *physical* result is a table of
`q,omega(q),dW/domega,source-normalized exterior
Green-function residue(q)`, for matching nonrotating
Einstein–Vlasov solutions and identical sources,
extending toward q=.95/.99. It must show whether
the long-lived pole vanishes, migrates, merges,
hits the continuum, or becomes invisible at
infinity. Only then is an astrophysical
detector forecast scientifically justified.

The already submitted manuscript is unchanged;
this is a separate research branch.

# Blind multiple-shooting verification

This directory contains a separate blind numerical verification of the
`q = 0.1` matter-led global mode.

The verification was performed within the **Podatkovna Infrastruktura i
Proširena Inteligencija (PiPi)** institutional research project at the
Faculty of Organization and Informatics, University of Zagreb. The individual
verifier requested not to be named.

The verification received the physical equations, the publication-resolution
spectral input, and a broad complex-frequency search window. The target
matter-mode frequency, manuscript, and production solver were withheld until
the search and numerical audits were complete.

Anthropic Claude was used as a coding assistant during implementation. The
researchers defined the physical specification, blind protocol, search window,
and validation criteria. The code, failed searches, numerical outputs, and
audit files were inspected before unblinding.

## Numerical method

The vacuum problem uses real-axis Riccati multiple shooting with adaptive
Runge--Kutta integration. No exterior complex scaling or Chebyshev
Dirichlet-to-Neumann solver is used. The polar Einstein--Vlasov bulk matrix is
re-derived symbolically, and the frequency-dependent retarded constitutive
response is evaluated from the supplied spectral nodes at each trial
frequency.

The blind search found one stable simple pole in the prescribed window,

```text
M omega = 0.377514811497 - 0.020871653994 i
```

The value was compared with the production result only after the blind search
and audits were complete.

## Reproducing the verification

Install the Python dependencies from the repository root and place the
publication-resolution file

```text
blind_q010_spectral_nodes.npz
```

in this directory. The spectral file is distributed as Supplementary Data 1
with the manuscript because it is not kept in this minimal Git repository.
Then run

```bash
python run_blind_search.py
python winding_audit.py
```

The first command performs the full seed search and writes the root-search and
resolution tables. The second performs the argument-principle winding audit.

## Files

- `independent_solver.py` - global multiple-shooting condition.
- `bulk_blind.py` - separately derived Einstein--Vlasov bulk system.
- `spectral.py` - retarded spectral continuation.
- `run_blind_search.py` - blind frequency search.
- `winding_audit.py` - argument-principle audit.
- `numerical_method.md` - numerical details and independence statement.
- `schwarzschild_benchmark.csv` - vacuum benchmark.
- `all_frequency_seeds.csv` - all attempted seeds, including failures.
- `all_converged_roots.csv` - deduplicated roots.
- `resolution_convergence.csv` - outer-boundary stability.
- `bulk_lambda_audit.csv` - full 7x4 SVD residuals.
- `winding_audit.csv` and `winding_audit_result.txt` - winding-number data.
- `final_verdict.txt` - compact audit summary.

This directory is provided as an additional numerical-method verification. It
is not presented as an external replication by an independent research group.

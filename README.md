# HBW_GW

Numerical code and compact audit data for our study of polar perturbations of smooth, self-gravitating massless Einstein--Vlasov atmospheres around Schwarzschild black holes.

The repository is intentionally small. It contains publication-facing modules from the verified calculation pipeline, compact result tables, and separate numerical checks of the global `q = 0.1` mode. Development runs and manuscript files are not included.

## Main result

For the smooth atmosphere with mass ratio `q = 0.1`, the matter-led polar mode is

```text
M omega = 0.377514814460 - 0.020871654486 i
```

The mode is also recovered in independently constructed smooth atmospheres at

```text
q = 0.05, 0.075, 0.10, 0.15, 0.20
```

with unit argument-principle winding and stable global boundary-condition audits. The accepted family is listed in `results/q_family_results.csv`.

At `q = 0.10` and `q = 0.20`, localized polar sources placed entirely in the outer Schwarzschild vacuum have non-zero overlap with the matter-led pole and finite retarded Green-function residues. Compact results are in `results/source_excitation_summary.csv`, `results/point_source_residues.csv`, and `results/smooth_source_residues.csv`.

## Separate numerical checks

A Chebyshev pseudospectral boundary-value implementation reproduces the `q = 0.1` global pole without using the transfer-matrix, homotopy, or exterior-complex-scaling machinery of the primary solver. The two calculations agree in complex frequency to `1.4e-8`. The solver is in `code/independent_pseudospectral_qnm.py`.

A second blind verification was performed within the **Podatkovna Infrastruktura i Proširena Inteligencija (PiPi)** institutional research project at the Faculty of Organization and Informatics. The target matter-mode frequency, manuscript, and production solver were withheld until the search and numerical audits were complete. The real-axis multiple-shooting calculation recovered

```text
M omega = 0.377514811497 - 0.020871653994 i
```

which differs from the primary result by `3.0e-9`. Its winding number is one and the full bulk relative singular-value residuals are below `8e-15`.

The complete code, all attempted seeds, failed searches, boundary checks, and winding audit are in `verification/blind_multiple_shooting/`. Anthropic Claude was used as a coding assistant during implementation of this blind check. The researchers defined the physical specification, blind protocol, search window, and validation criteria and inspected the code and audit outputs before unblinding.

## Microscopic response-profile robustness

At `q = 0.1` the retarded spectral weights were redistributed smoothly across the publication shells using tilt, central, and two-lobe deformations. The shell weights span approximately `0.68` to `1.39`, while the same simple global pole remains present in all seven cases. The largest relative complex-frequency shift is `6.45e-6`; every case retains unit winding, full Einstein--Vlasov residual control, and stable outer-vacuum matching.

This is a microscopic response-profile robustness test with the equilibrium background held fixed. It is not presented as a second self-consistent stationary Einstein--Vlasov family. The reproducible script is `verification/blind_multiple_shooting/profile_shape_robustness.py`; compact results are in `results/profile_shape_robustness.csv` and `results/profile_shape_robustness_stability.csv`.

## Repository layout

```text
code/          publication-facing calculation modules
results/       compact numerical audit tables
verification/  separate numerical verification packages
```

The filenames in `code/` retain internal phase labels used during verification so that the modules can be traced to the corresponding successful numerical stages. Exploratory scans and failed development branches are not included in the primary calculation package.

## Environment

Python 3.11 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The smooth global determinant is implemented in `code/hbh_phase_xxic9_final_global_smooth_qnm_v4.py`. The external-source Green-function calculation is in `code/hbh_np_source_excitation_test_v2.py`. The independent Chebyshev reproduction is in `code/independent_pseudospectral_qnm.py`.

The publication-resolution `q = 0.1` spectral-node array `blind_q010_spectral_nodes.npz` is required for the blind multiple-shooting and profile-shape checks. It should be placed in `verification/blind_multiple_shooting/`. The expected SHA256 checksum is recorded there so the publication input can be verified byte-for-byte.

## Scope

The repository reports the publication-grade family through `q = 0.20`. Exploratory calculations outside that verified range are deliberately not included.

# HBW_GW

Numerical code and compact audit data for our study of polar perturbations of smooth, self-gravitating massless Einstein--Vlasov atmospheres around Schwarzschild black holes.

This repository is intentionally small. It contains publication-facing modules from the verified calculation pipeline together with compact tables supporting the main numerical claims. Development runs and manuscript files are not included.

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

## Independent numerical checks

A Chebyshev pseudospectral boundary-value implementation reproduces the `q = 0.1` global pole without using the transfer-matrix, homotopy, or exterior-complex-scaling machinery of the primary solver. The two calculations agree in complex frequency to `1.4e-8`. The solver is in `code/independent_pseudospectral_qnm.py`.

A second, deliberately blind from-scratch real-axis multiple-shooting implementation was not supplied the target matter-mode frequency or production solver. After searching the prescribed complex-frequency window it recovered

```text
M omega = 0.377514811497 - 0.020871653994 i
```

which differs from the primary result by `3.0e-9`; its winding number is one and the full bulk relative singular-value residuals are below `8e-15`. The compact summary is in `results/blind_multiple_shooting_reproduction.csv`. The complete blind source package and all attempted seeds/audits are supplied with the manuscript as Supplementary Code 1.

## Repository layout

```text
code/       publication-facing calculation modules
results/    compact numerical audit tables
data/       notes on publication data
```

The filenames in `code/` retain internal phase labels used during verification so that modules can be traced to the successful numerical stages. Exploratory scans and failed development branches are not included.

## Environment

Python 3.11 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The smooth global determinant is implemented in `code/hbh_phase_xxic9_final_global_smooth_qnm_v4.py`. The external-source Green-function calculation is in `code/hbh_np_source_excitation_test_v2.py`. The independent Chebyshev reproduction is in `code/independent_pseudospectral_qnm.py`.

The publication-resolution `q = 0.1` spectral-node arrays used for the retarded constitutive response and blind reproduction are supplied with the manuscript as Supplementary Data 1. An archival data DOI can be added with the published version.

## Scope

The repository reports the publication-grade family through `q = 0.20`. Exploratory calculations outside that verified range are deliberately not included.

# HBW_GW

Numerical code and compact audit data for our study of polar perturbations of smooth, self-gravitating massless Einstein--Vlasov atmospheres around Schwarzschild black holes.

This repository is intentionally small. It contains publication-facing modules distilled from the verified calculation pipeline together with the tables that support the main numerical claims. Development runs and manuscript files are not included.

## Main numerical result

For the smooth atmosphere with mass ratio `q = 0.1`, the matter-led polar mode is

```text
M omega = 0.377514814460 - 0.020871654486 i
```

The mode is also recovered in independently constructed smooth atmospheres at

```text
q = 0.05, 0.075, 0.10, 0.15, 0.20
```

with unit argument-principle winding and stable global boundary-condition audits. The accepted family is listed in `results/q_family_results.csv`.

At `q = 0.10` and `q = 0.20` we additionally tested direct source coupling. Localized polar sources placed entirely in the outer Schwarzschild vacuum have non-zero overlap with the matter-led pole and finite retarded Green-function residues. The compact results are in `results/source_excitation_summary.csv`, `results/point_source_residues.csv`, and `results/smooth_source_residues.csv`.

## Repository layout

```text
code/       publication-facing calculation modules
results/    compact numerical audit tables
data/       note on archived spectral caches
```

The filenames in `code/` retain the internal phase labels used during verification so that the modules can be traced to the corresponding successful numerical stages. Development orchestration, exploratory scans, and failed branches have been removed from the public version.

## Environment

Python 3.11 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The smooth global determinant is implemented in `code/hbh_phase_xxic9_final_global_smooth_qnm_v4.py`. The external-source Green-function calculation is in `code/hbh_np_source_excitation_test_v2.py`.

The large spectral-node caches used by these calculations are not committed to this minimal Git repository. They will be deposited with the archival dataset accompanying the paper. Small numerical outputs needed to inspect the reported results are included under `results/`.

## Scope

The repository reports the publication-grade family through `q = 0.20`. Exploratory calculations outside that verified range are deliberately not included.

A manuscript citation and archival data DOI will be added with the paper.
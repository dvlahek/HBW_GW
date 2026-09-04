# HBW_GW

Publication-facing numerical code and compact source data for the study **“Cold collisionless matter can retain hidden dynamics.”**

The physical question is broader than one black-hole mode: a collisionless state can become macroscopically cold in an equilibrium moment while retaining a finite derivative-weighted linear response. The black-hole calculation is a self-consistent strong-field realization of that cold-response mechanism in the leading radiative gravitational sector, `ell = 2`.

The repository is intentionally compact. It contains the publication-facing global-response code, source-data tables, independent numerical checks and the audit data needed to inspect the reported claims. Development scans are not included.

## Main strong-field result

For the smooth Einstein--Vlasov atmosphere with mass ratio

```text
q = m/M = 0.10
```

the matter-led polar mode is

```text
M omega = 0.377514814460 - 0.020871654486 i
```

The mode is recovered in independently rebuilt smooth atmospheres at

```text
q = 0.05, 0.075, 0.10, 0.15, 0.20
```

with stable global boundary-condition and simple-pole audits. The accepted family is in `results/q_family_results.csv`.

## Second self-consistent profile family

The manuscript also changes the stationary phase-space profile itself to

```text
F_alt(E,L) = C (E0-E)_+^2 (L-L0)_+^2.
```

Its high-resolution continuation path is selected using stationary proper-shell bracketing before any spectral or QNM result is evaluated. The independent family approaches the same radially cold endpoint while retaining finite kinetic susceptibility. Its final global pole is

```text
M omega_alt = 0.375287456824 - 0.021057162748 i
```

which is `0.591145%` from the reference `q = 0.10` pole in complex frequency. Compact cold-orbit, seed, stability and QNM audit tables are in `results/second_family_*`; the interpretation and procedure are summarized in `verification/second_family/README.md`.

This is a profile-robustness result, not a universality theorem over arbitrary collisionless distributions.

## Identical exterior forcing

The source-coupling calculation places perturbations entirely in the common outer Schwarzschild region. The ordinary Schwarzschild `ell = 2` fundamental is recomputed in the same Zerilli normalization,

```text
M omega_vac = 0.373671684352 - 0.088962315642 i.
```

For the same source, observer and master-function normalization, the matter-led / vacuum Green-function pole-residue ratio is:

```text
q = 0.10: 0.96--2.50% for six point sources
          1.17--1.79% for three smooth sources

q = 0.20: 1.60--3.98% for six point sources
          1.93--2.89% for three smooth sources
```

A real finite-duration Gaussian-cosine burst scan is also reported. These values compare the two isolated pole contributions; they are not a numerical-relativity merger waveform or a detector forecast. Source data are in `results/source_data/` and the compact post-processing reproduction is in `verification/vacuum_normalized_excitation/`.

## Independent global checks

A Chebyshev pseudospectral boundary-value implementation reproduces the `q = 0.10` global pole without the production transfer-matrix or exterior-complex-scaling construction. The complex-frequency agreement is approximately `1.4e-8`. The solver is `code/independent_pseudospectral_qnm.py`.

A separate blind real-axis multiple-shooting verification within the **Podatkovna Infrastruktura i Proširena Inteligencija (PiPi)** institutional research project recovered

```text
M omega = 0.377514811497 - 0.020871653994 i
```

before unblinding, differing from the primary result by about `3.0e-9`. Its code, search seeds, failed searches, winding audit and bulk checks are in `verification/blind_multiple_shooting/`. Anthropic Claude was used as a coding assistant during implementation of that blind check; the physical specification, blind protocol and validation criteria were defined independently and the outputs were inspected before unblinding.

## Source data

`results/source_data/` contains the compact numerical data underlying the four main figures, together with the vacuum-normalized source comparison and real-burst scan.

The publication-resolution `q = 0.10` spectral-node binary is not stored in this minimal Git repository. It is supplied with the manuscript as Supplementary Data 1. The expected checksum for the blind reproduction is recorded in `verification/blind_multiple_shooting/SPECTRAL_DATA_SHA256.txt`.

## Repository layout

```text
code/          publication-facing calculation modules
results/       compact numerical results and source data
verification/  independent and post-processing verification packages
```

The filenames in `code/` retain internal phase labels so successful numerical stages remain traceable.

## Environment

Python 3.11 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The smooth global determinant is implemented in `code/hbh_phase_xxic9_final_global_smooth_qnm_v4.py`. The exterior-source Green-function calculation is in `code/hbh_np_source_excitation_test_v2.py`. The independent Chebyshev reproduction is in `code/independent_pseudospectral_qnm.py`.

## Additional exploratory detector injection

The repository retains an earlier idealized Gaussian-noise injection/recovery study under `verification/detector_injection/` and `results/detector_*` for transparency. It is **not used to support the current manuscript’s source-response claim** and is not an astrophysical detectability forecast.

## Scope

The publication-grade Einstein--Vlasov mass-ratio family reported here is restricted to `q = 0.05--0.20`. Exploratory calculations outside that verified range are deliberately not used in the manuscript.

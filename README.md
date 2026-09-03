# HBW_GW

Numerical code and compact audit data for our study of polar perturbations of smooth, self-gravitating massless Einstein--Vlasov atmospheres around Schwarzschild black holes.

The repository is intentionally small. It contains the verified scripts used for the final global mode calculation, two cached spectral datasets used in the source-coupling check, and the numerical tables that support the main robustness statements. Development runs and manuscript files are not included.

## Main numerical result

For the smooth atmosphere with mass ratio `q = 0.1`, the matter-led polar mode is

```text
M omega = 0.377514814460 - 0.020871654486 i
```

The same mode is recovered in independently constructed smooth atmospheres at

```text
q = 0.05, 0.075, 0.10, 0.15, 0.20
```

with unit argument-principle winding and stable global boundary-condition audits. The tabulated values are in `results/q_family_results.csv`.

We also tested direct source coupling at `q = 0.10` and `q = 0.20`. Localized polar sources placed entirely in the outer Schwarzschild vacuum have non-zero overlap with the matter-led pole and finite retarded Green-function residues. These results are summarized in `results/source_excitation_summary.csv`.

## Repository layout

```text
code/       verified calculation scripts
data/       cached spectral inputs for q=0.10 and q=0.20
results/    compact publication-facing audit tables
```

The filenames in `code/` retain the internal phase labels used during verification. They are kept unchanged so the public files can be traced directly to the successful numerical runs.

## Environment

Python 3.11 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

## Reproducing the q = 0.1 global mode

Run

```bash
python reproduce_q010.py
```

This rebuilds the small C7/C8 input archives expected by the verified global solver and runs the final smooth QNM calculation from the cached `q=0.1` spectral data.

## Reproducing the source-coupling test

Run

```bash
python reproduce_source_coupling.py
```

This assembles a minimal two-case family archive from the cached `q=0.10` and `q=0.20` spectral data and executes the verified retarded Green-function source-coupling calculation.

The generated pole contribution is a linear source-coupling diagnostic. It is not a merger waveform or a detector-sensitivity forecast.

## Notes on scope

The repository reports the publication-grade family through `q=0.20`. Exploratory calculations outside that verified range are deliberately not included here.

A manuscript citation and archival data DOI will be added with the paper.

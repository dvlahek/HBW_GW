#!/usr/bin/env python3
"""Fast publication-input and manuscript-claim smoke test.

This is not a full QNM re-solve. It verifies the public spectral input, checks
the constitutive response used by the separately implemented pseudospectral
solver against the production evaluator, and verifies representative recorded
results against the deposited CSVs.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from publication_inputs import load_publication_case
import hbh_phase_xxic8_landau_continuation as landau
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as production

radial, carrier, row = load_publication_case(0.10)
assert len(radial) == 1440, f"Unexpected radial branch count: {len(radial)}"
assert len(carrier) == 80, f"Unexpected carrier branch count: {len(carrier)}"

omega = complex(float(row["root_re"]), float(row["root_im"]))
reference = 0.377514814460 - 0.020871654486j
assert abs(omega - reference) < 5e-7, "q-family root differs from published q010 reference"

check = landau.evaluate_chi(omega, radial, carrier, second_sheet=True)
expected = production.smooth_chi(omega, radial, carrier)
for key in ("PH", "PJJ", "PJX", "chiJ", "chiX"):
    delta = abs(check[key] - expected[key])
    scale = max(1., abs(expected[key]))
    assert delta / scale < 1e-11, f"Constitutive evaluator disagreement for {key}: {delta}"

for key, csv_re, csv_im in (
    ("chiJ", "chiJ_re", "chiJ_im"),
    ("chiX", "chiX_re", "chiX_im"),
):
    published = complex(float(row[csv_re]), float(row[csv_im]))
    assert abs(check[key] - published) < 2e-3, (
        f"{key} inconsistent with compact publication q-family table"
    )

qfamily = pd.read_csv(ROOT / "results" / "q_family_results.csv")
assert np.allclose(sorted(qfamily.q.tolist()), [0.05, 0.075, 0.1, 0.15, 0.2])
alt = pd.read_csv(ROOT / "results" / "second_family_qnm_compare.csv")
assert not alt.empty
point = pd.read_csv(ROOT / "results" / "source_data" / "Fig4_point_excitation_ratios.csv")
smooth = pd.read_csv(ROOT / "results" / "source_data" / "Fig4_smooth_excitation_ratios.csv")
for frame in (point, smooth):
    assert set(np.round(frame.q.unique(), 2)) == {0.1, 0.2}
    derived = frame.kinetic_residue_abs / frame.vacuum_residue_abs
    assert np.max(np.abs(derived - frame.kinetic_to_vacuum_residue_ratio)) < 5e-13

try:
    load_publication_case(0.20)
except FileNotFoundError as exc:
    assert "--family" in str(exc)
else:
    raise AssertionError("The q=0.20 full spectral input must not be represented as deposited")

print("PUBLICATION_READINESS_INPUTS_PASS")
print(f"q010: {len(radial)} radial branches, {len(carrier)} carrier shells")
print(f"q010 root in family table: {omega}")
print("Production / independent retarded susceptibility: consistent")
print("q020 full recomputation: explicitly requires external full-family ZIP")

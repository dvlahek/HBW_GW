#!/usr/bin/env python3
"""Locate and validate the publication inputs for the global QNM calculations.

The openly deposited q=0.10 spectral nodes are sufficient for q=0.10
recomputations. Other q values require the original, separately supplied
full-family archive; a compact CSV cannot replace missing spectral nodes.
"""

import hashlib
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import hbh_phase_xxic8_landau_continuation as landau

ROOT = Path(__file__).resolve().parents[1]
Q010_SPECTRAL = ROOT / "verification" / "blind_multiple_shooting" / "blind_q010_spectral_nodes.npz"
Q_FAMILY_TABLE = ROOT / "results" / "q_family_results.csv"
Q010_SHA256 = "04016df710591f61023338b20117b9b3fd046f4fc3ca886060f2b4a951d9114f"
SPECTRAL_KEYS = {"nu_J", "q_JJ", "q_JX", "n_r", "n_phi", "shell_J", "nu_H", "q_HH", "shell_H"}


def _read_arrays(raw):
    with np.load(io.BytesIO(raw), allow_pickle=False) as npz:
        if set(npz.files) != SPECTRAL_KEYS:
            raise ValueError("Unexpected spectral-node array schema: " + repr(npz.files))
        return {key: np.array(npz[key]) for key in npz.files}


def _row(table, q):
    matches = table[np.isclose(table.q.to_numpy(dtype=float), float(q), rtol=0, atol=1e-9)]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one published mass-ratio row for q={q}.")
    return matches.iloc[0].to_dict()


def load_publication_case(q, family_zip=None):
    """Return (radial catalogue, carrier catalogue, q-family result row).

    Without --family, only q=0.10 is supported, using the checksum-verified
    public NPZ and the accepted q-family CSV from this repository.
    With --family, read the legacy complete family ZIP (no private data are
    embedded here). In both cases the same spectral catalogue is constructed.
    """
    q = float(q)
    if family_zip is None:
        if not np.isclose(q, 0.10, rtol=0, atol=1e-9):
            raise FileNotFoundError(
                f"Full spectral data for q={q:g} are not in this repository. "
                "Use --family PATH_TO_FULL_FAMILY.zip for this mass ratio, "
                "or reproduce the q=0.10 calculation from the public NPZ."
            )
        raw = Q010_SPECTRAL.read_bytes()
        actual_sha = hashlib.sha256(raw).hexdigest()
        if actual_sha != Q010_SHA256:
            raise ValueError(f"Publication q=0.10 spectral SHA256 mismatch: {actual_sha}")
        row = _row(pd.read_csv(Q_FAMILY_TABLE), q)
    else:
        path = Path(family_zip)
        if not path.is_file():
            raise FileNotFoundError(f"Full-family spectral archive not found: {path}")
        with zipfile.ZipFile(path, "r") as zf:
            names = zf.namelist()
            result_names = [n for n in names if n.endswith("consistent_q_family_results.csv")]
            if len(result_names) != 1:
                raise ValueError("Expected exactly one consistent_q_family_results.csv in full-family ZIP.")
            row = _row(pd.read_csv(zf.open(result_names[0])), q)
            n_shell = int(row["N"])
            tag = f"{q:.3f}".replace(".", "p")
            cache_names = [n for n in names if n.endswith(f"spectral_cache/q_{tag}_N{n_shell}.npz")]
            if len(cache_names) != 1:
                raise ValueError(f"Missing unique spectral cache for q={q:g}, N={n_shell} in full-family ZIP.")
            raw = zf.read(cache_names[0])
    radial, carrier = landau.build_catalog(_read_arrays(raw))
    if not radial or not carrier:
        raise ValueError("Empty kinetic spectral catalogue.")
    return radial, carrier, row

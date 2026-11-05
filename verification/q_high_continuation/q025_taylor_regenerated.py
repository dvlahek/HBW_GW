#!/usr/bin/env python3
"""New q=.25 80-shell C1-C7 retarded C8 Taylor sources at nxi 641 and 801.

Generated from the q=.25 scaled-delta source background by independent
C7 regeneration. Validates the exact frozen bytes, source SHA, degree
and held-out error against the same fixed discrete C8 causal-sheet anchor.
Not a physically unique continuation or a QNM root.
"""
import base64,hashlib,json,zlib
import numpy as np
PAYLOAD_B64="""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"""
COMPRESSED_SHA256="58a0fbb91aa32fdb8263e362c3884e0d5c8353b5b97ef3949be8b6da2e4a103a"
JSON_SHA256="26d15ca0261491415041d35612e720446d9135605f0cfd92cc4a6b0b24965441"
PHYSICAL_NPZ_SHA={
 "n641_a0":"e1bdfb360f75343c56df9a559c4b28eae11c4e472b5b010bbb8d32994fde85f3",
 "n801_a0":"998804e476d7208d315f94337298068c0f3476fe2bd5d6a70dd2654a97c87298",
 "n801_am0005":"998804e476d7208d315f94337298068c0f3476fe2bd5d6a70dd2654a97c87298"
}
def cases():
    raw=base64.b64decode(PAYLOAD_B64,validate=True)
    if hashlib.sha256(raw).hexdigest()!=COMPRESSED_SHA256:
        raise ValueError("q025 compressed Taylor checksum mismatch")
    source=zlib.decompress(raw)
    if hashlib.sha256(source).hexdigest()!=JSON_SHA256:
        raise ValueError("q025 decompressed Taylor source checksum mismatch")
    d=json.loads(source)
    if set(d)!=set(PHYSICAL_NPZ_SHA):
        raise ValueError("Missing an independent q025 high-resolution configuration")
    for label,x in d.items():
        if (x["q"]!=.25 or x["physical_sheet_not_certified"] is not True
            or x["not_a_QNM"] is not True or
            x["nxi"]!=(641 if label=="n641_a0" else 801)
            or x["source_spectral_sha256"]!=PHYSICAL_NPZ_SHA[label]
            or x["taylor_order"]!=12 or x["strict_validity_radius"]!=.004
            or x["worst_heldout_exact_cauchy_relative"]>1e-9):
            raise ValueError("q025 physical source/analytic model claims invalid")
        for key in ("chiJ","chiX"):
            v=np.array([complex(*item) for item in x[key+"_coeff_scaled"]])
            if v.shape!=(13,) or not np.all(np.isfinite(v)):
                raise ValueError("Invalid complex susceptibility Taylor vector")
    return d
def evaluate(case,w):
    center=complex(*case["center_Momega"])
    dz=(complex(w)-center)/case["strict_validity_radius"]
    if abs(dz)>1:
        raise ValueError("Analytic continuation outside verified q025 disk")
    ret={}
    for key in ("chiJ","chiX"):
        vals=np.array([complex(*x) for x in case[key+"_coeff_scaled"]])
        ret[key]=complex(np.polynomial.polynomial.polyval(dz,vals))
    return ret

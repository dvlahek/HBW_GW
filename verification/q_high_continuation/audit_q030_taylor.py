#!/usr/bin/env python3
"""Strict frozen Taylor coefficient provenance and schema inspection.

This is a data-schema gate, not an extra-mode detection or a global root.
"""
from pathlib import Path
import json
import math
import numpy as np
from q030_taylor_coefficients import cases

def audit():
    data=cases()
    summary={}
    for name, c in sorted(data.items()):
        print("CASE_KEYS",name,sorted(c))
        for k,v in c.items():
            if isinstance(v,dict):
                print("CHILD_KEYS",name,k,sorted(v))
        summary[name]={
            "q":c["q"],"nxi":c["nxi"],
            "physical_spectral_sha256":c["physical_spectral_sha256"],
            "taylor_order":c["taylor_order"],
            "strict_validity_radius":c["strict_validity_radius"],
            "errors":c["errors"],
        }
    return summary

if __name__=="__main__":
    print("Q030_FROZEN_TAYLOR_SOURCE_AUDIT",json.dumps(audit(),separators=(",",":")))

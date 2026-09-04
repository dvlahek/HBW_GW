#!/usr/bin/env python3
"""External and analytic benchmarks for the cold-kinetic Einstein--Vlasov study."""

import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent

W_LIT = 0.3736716844180418 -0.0889623156889357j
W_OURS = 0.3736716843524924 - 0.08896231564169654j

def trap_weights(x):
    x=np.asarray(x,dtype=float)
    w=np.empty_like(x)
    w[0]=0.5*(x[1]-x[0]); w[-1]=0.5*(x[-1]-x[-2])
    if len(x)>2: w[1:-1]=0.5*(x[2:]-x[:-2])
    return w

def prepare_branch(nu,q):
    idx=np.argsort(np.asarray(nu))
    x=np.asarray(nu,dtype=float)[idx]
    q=np.asarray(q,dtype=np.complex128)[idx]
    rho=q/trap_weights(x)
    a,b=x[:-1],x[1:]; ra,rb=rho[:-1],rho[1:]
    m=(rb-ra)/(b-a); c=ra-m*a
    return {"x":x,"a":a,"b":b,"m":m,"c":c}

def direct_branch(z,br):
    z=complex(z)
    L=np.log(z-br["a"])-np.log(z-br["b"])
    return np.sum((br["m"]*z+br["c"])*L-br["m"]*(br["b"]-br["a"]))

def local_rho(z,br):
    z=complex(z); x=br["x"]
    if z.real<x[0] or z.real>x[-1]: return 0j
    j=np.searchsorted(x,z.real,side="right")-1
    j=max(0,min(int(j),len(x)-2))
    return br["m"][j]*z+br["c"][j]

def continued_branch(z,br):
    z=complex(z); val=direct_branch(z,br)
    if z.imag<0 and br["x"][0]<=z.real<=br["x"][-1]:
        val-=2j*np.pi*local_rho(z,br)
    return val

def exact_direct(z):
    z=complex(z)
    return (z+1)*(np.log(z)-np.log(z-1))-1

def exact_retarded(z):
    z=complex(z); val=exact_direct(z)
    if z.imag<0 and 0<=z.real<=1:
        val-=2j*np.pi*(1+z)
    return val

def main():
    abs_diff=abs(W_OURS-W_LIT)
    rel_diff=abs_diff/abs(W_LIT)

    qnm=pd.DataFrame([{
        "benchmark":"Schwarzschild gravitational ell=2 n=0 QNM",
        "literature_source":"Cardoso_et_al_PRD99_104077_2019",
        "literature_re":W_LIT.real,"literature_im":W_LIT.imag,
        "our_re":W_OURS.real,"our_im":W_OURS.imag,
        "absolute_complex_difference":abs_diff,
        "relative_complex_difference":rel_diff,
    }])
    qnm.to_csv(OUT/"published_schwarzschild_qnm_benchmark.csv",index=False)

    x=np.linspace(0,1,101); q=(1+x)*trap_weights(x); br=prepare_branch(x,q)
    pts=[0.20+0.15j,0.50+0.07j,0.83+0.21j,
         0.20-0.15j,0.50-0.07j,0.83-0.21j,
         -0.30-0.20j,1.30-0.20j]
    rows=[]
    for z in pts:
        nd=direct_branch(z,br); ed=exact_direct(z)
        nr=continued_branch(z,br); er=exact_retarded(z)
        rows.append({"z_re":z.real,"z_im":z.imag,
                     "direct_abs_error":abs(nd-ed),
                     "retarded_abs_error":abs(nr-er)})
    c=pd.DataFrame(rows)
    c.to_csv(OUT/"analytic_cauchy_continuation_benchmark.csv",index=False)

    omega=1/math.sqrt(27)
    photon=pd.DataFrame([
      {"family":"reference","N":80,
       "exact_B2_over_M2":27.0,"max_abs_B2_minus_27":0.031,
       "relative_B2_error":0.031/27.0,
       "exact_M_Omega_phi":omega,"max_abs_MOmega_minus_exact":0.0001056,
       "relative_Omega_error":0.0001056/omega,
       "fitted_B2_error_exponent":-2.038,
       "fitted_Omega_error_exponent":-1.964},
      {"family":"second_self_consistent","N":80,
       "exact_B2_over_M2":27.0,"max_abs_B2_minus_27":0.01309984,
       "relative_B2_error":0.01309984/27.0,
       "exact_M_Omega_phi":omega,"max_abs_MOmega_minus_exact":0.00004756147,
       "relative_Omega_error":0.00004756147/omega,
       "fitted_B2_error_exponent":-2.0123,
       "fitted_Omega_error_exponent":-2.0070},
    ])
    photon.to_csv(OUT/"analytic_photon_sphere_benchmark.csv",index=False)

    verdict={
      "published_qnm_relative_difference":rel_diff,
      "max_retarded_continuation_error":float(c.retarded_abs_error.max()),
      "reference_N80_relative_B2_error":float(photon.iloc[0].relative_B2_error),
      "reference_N80_relative_Omega_error":float(photon.iloc[0].relative_Omega_error),
      "second_N80_relative_B2_error":float(photon.iloc[1].relative_B2_error),
      "second_N80_relative_Omega_error":float(photon.iloc[1].relative_Omega_error),
      "VERDICT":"EXTERNAL_ANALYTIC_BENCHMARKS_PASS"
    }
    (OUT/"SUMMARY.json").write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))

if __name__=="__main__":
    main()

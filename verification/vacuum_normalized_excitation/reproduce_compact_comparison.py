#!/usr/bin/env python3
"""Reproduce compact vacuum-normalized source ratios and real-burst scan.

This script uses the publication tables in results/.  It does not rebuild the
Einstein--Vlasov spectral cache.  The kinetic residues are the verified outputs
of code/hbh_np_source_excitation_test_v2.py; the Schwarzschild reference
residues and root audit are stored in results/source_data/.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SD = ROOT / "results" / "source_data"

WVAC = 0.3736716843524924 - 0.08896231564169654j
WKIN = {
    0.10: 0.3775148159474828 - 0.020871633919164j,
    0.20: 0.3797013686231943 - 0.0239854837986382j,
}
CARRIERS = [0.34, 0.36, WVAC.real, 0.39, 0.42]
SIGMAS = [2.0, 5.0, 10.0]

def real_burst_factor(w, sigma_t, omega_c):
    T = 8.0 * sigma_t
    t = np.linspace(0.0, T, 16001)
    t0 = 4.0 * sigma_t
    S = np.exp(-0.5*((t-t0)/sigma_t)**2) * np.cos(omega_c*(t-t0))
    S /= math.sqrt(max(float(np.trapezoid(S*S,t)),1e-300))
    I = np.trapezoid(S*np.exp(1j*w*t),t)
    return np.exp(-1j*w*T)*I

def main():
    point = pd.read_csv(SD / "Fig4_point_excitation_ratios.csv")
    smooth = pd.read_csv(SD / "Fig4_smooth_excitation_ratios.csv")

    for df in (point, smooth):
        calc = df.kinetic_residue_abs / df.vacuum_residue_abs
        err = np.max(np.abs(calc - df.kinetic_to_vacuum_residue_ratio))
        assert err < 5e-13

    rows=[]
    for q in (0.10,0.20):
        gain=abs(WVAC.imag)-abs(WKIN[q].imag)
        sq=smooth[smooth.q==q]
        for wc in CARRIERS:
            for sig in SIGMAS:
                tr=abs(real_burst_factor(WKIN[q],sig,wc))/abs(real_burst_factor(WVAC,sig,wc))
                for _,r in sq.iterrows():
                    rend=float(r.kinetic_to_vacuum_residue_ratio)*tr
                    rows.append({
                        'q':q,'source_name':r.source_name,'carrier_Momega':wc,
                        'sigma_t_over_M':sig,'kinetic_to_vacuum_at_burst_end':rend,
                        'ratio_25M_after_burst':rend*math.exp(gain*25.0),
                    })
    calc=pd.DataFrame(rows)
    summary=[]
    for q in (0.10,0.20):
        for sig in SIGMAS:
            x=calc[(calc.q==q)&(calc.sigma_t_over_M==sig)]
            summary.append({
                'q':q,'sigma_t_over_M':sig,
                'burst_end_ratio_min':x.kinetic_to_vacuum_at_burst_end.min(),
                'burst_end_ratio_max':x.kinetic_to_vacuum_at_burst_end.max(),
                'ratio_25M_min':x.ratio_25M_after_burst.min(),
                'ratio_25M_max':x.ratio_25M_after_burst.max(),
            })
    summary=pd.DataFrame(summary)
    ref=pd.read_csv(SD / "ringdown_real_burst_summary.csv")
    m=summary.merge(ref,on=['q','sigma_t_over_M'],suffixes=('_calc','_ref'))
    cols=['burst_end_ratio_min','burst_end_ratio_max','ratio_25M_min','ratio_25M_max']
    err=max(np.max(np.abs(m[f'{c}_calc']-m[f'{c}_ref'])) for c in cols)
    assert err < 2e-12
    print('VACUUM_NORMALIZED_COMPACT_REPRODUCTION_PASS')
    print(f'max summary difference = {err:.3e}')
    print(summary.to_string(index=False))

if __name__=='__main__':
    main()

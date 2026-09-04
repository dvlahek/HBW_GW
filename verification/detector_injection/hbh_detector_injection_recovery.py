#!/usr/bin/env python3
"""
HBH matter-led ringdown: idealized detector-level injection/recovery study.

This is a controlled synthetic feasibility test, not an analysis of real GW data.
It injects the Schwarzschild l=2 vacuum fundamental plus the predicted q=0.1
matter-led mode into Gaussian noise shaped by an analytic Advanced-LIGO design
PSD proxy. It compares a vacuum-only ringdown model with a vacuum + fixed-
frequency matter-mode model and repeats the comparison at progressively later
ringdown start times.
"""
from __future__ import annotations
import argparse, csv, math, zipfile
from pathlib import Path
import numpy as np

W_VAC = 0.373671684 - 0.088962316j
W_KIN = 0.377514814460 - 0.020871654486j
MTSUN_SI = 4.925490947e-6
MASSES = [30.0, 60.0, 100.0, 180.0]
SNRS = [12.0, 20.0, 30.0, 40.0]
RATIOS = [0.0, 0.02, 0.05, 0.10, 0.20, 0.30]
START_M = [0.0, 10.0, 20.0]

def mode_physical(w, mass_msun):
    Msec=MTSUN_SI*mass_msun
    return w.real/(2*math.pi*Msec), Msec/abs(w.imag)

def aligo_design_psd(f):
    f=np.asarray(f,float); x=np.maximum(f,1e-9)/215.0
    with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
        s=1e-49*(x**(-4.14)-5*x**(-2)+111*(1-x*x+0.5*x**4)/(1+0.5*x*x))
    s=np.where(np.isfinite(s)&(s>1e-50),s,1e6)
    return np.where(f>=20,s,1e6)

def tukey_end_window(t,duration,alpha=0.12):
    w=np.ones_like(t); edge=alpha*duration; m=t>(duration-edge)
    if np.any(m):
        u=(t[m]-(duration-edge))/edge; w[m]=0.5*(1+np.cos(np.pi*u))
    return w

def ringdown_basis(t,f_hz,tau_s,window):
    env=np.exp(-t/tau_s)*window; a=2*np.pi*f_hz*t
    return np.column_stack([env*np.cos(a),env*np.sin(a)])

def whiten_columns(X,fs):
    n=X.shape[0]; freqs=np.fft.rfftfreq(n,1/fs); psd=aligo_design_psd(freqs); scale=np.sqrt(psd); scale[0]=np.inf
    out=np.empty_like(X,float)
    for j in range(X.shape[1]):
        y=np.fft.irfft(np.fft.rfft(X[:,j])/scale,n=n); out[:,j]=y-np.mean(y)
    g=np.sqrt(np.mean(out*out))
    if not np.isfinite(g) or g<=0: raise RuntimeError('invalid whitening normalization')
    return out/g

def amp_phase(pair):
    a,b=float(pair[0]),float(pair[1]); return math.hypot(a,b), math.atan2(-b,a)

def build_mass_cache(mass,fs,duration):
    n=int(round(fs*duration)); t=np.arange(n)/fs; window=tukey_end_window(t,duration)
    fv,tv=mode_physical(W_VAC,mass); fk,tk=mode_physical(W_KIN,mass)
    Xraw=np.column_stack([ringdown_basis(t,fv,tv,window),ringdown_basis(t,fk,tk,window)])
    Xfull=whiten_columns(Xraw,fs); Msec=MTSUN_SI*mass; fits={}
    for sm in START_M:
        i0=int(np.searchsorted(t,sm*Msec,side='left')); X=Xfull[i0:,:]; Xv=X[:,:2]
        fits[sm]={'i0':i0,'X':X,'Xv':Xv,'pinv_X':np.linalg.pinv(X),'pinv_Xv':np.linalg.pinv(Xv),'n':len(X)}
    return {'Xfull':Xfull,'fits':fits,'f_vac_hz':fv,'tau_vac_s':tv,'f_kin_hz':fk,'tau_kin_s':tk}

def one_injection(rng,cache,mass,snr,ratio):
    X=cache['Xfull']; n=len(X); phv=rng.uniform(0,2*np.pi); phk=rng.uniform(0,2*np.pi)
    coeff=np.array([math.cos(phv),-math.sin(phv),ratio*math.cos(phk),-ratio*math.sin(phk)])
    h0=X@coeff; coeff*=snr/max(np.linalg.norm(h0),1e-300); h=X@coeff; y=h+rng.normal(0,1,n)
    rows=[]
    for sm,fc in cache['fits'].items():
        yy=y[fc['i0']:]; X2,X1=fc['X'],fc['Xv']; b1=fc['pinv_Xv']@yy; b2=fc['pinv_X']@yy
        r1=yy-X1@b1; r2=yy-X2@b2; rss1=float(r1@r1); rss2=float(r2@r2); nfit=fc['n']
        dbic=(rss1+2*math.log(nfit))-(rss2+4*math.log(nfit)); av,_=amp_phase(b2[:2]); ak,_=amp_phase(b2[2:])
        rows.append({'mass_msun':mass,'snr_full_window':snr,'ratio_true_at_t0':ratio,'start_M':sm,
                     'start_ms':1e3*sm*MTSUN_SI*mass,'f_vac_hz':cache['f_vac_hz'],'tau_vac_ms':1e3*cache['tau_vac_s'],
                     'f_kin_hz':cache['f_kin_hz'],'tau_kin_ms':1e3*cache['tau_kin_s'],'delta_bic':dbic,
                     'strong_two_mode':dbic>10,'positive_two_mode':dbic>0,'ratio_rec':ak/max(av,1e-300),
                     'rss_vac':rss1,'rss_two':rss2,'full_injected_snr':float(np.linalg.norm(h)),'fit_samples':nfit})
    return rows

def summarize(rows):
    groups={}
    for r in rows: groups.setdefault((r['mass_msun'],r['snr_full_window'],r['ratio_true_at_t0'],r['start_M']),[]).append(r)
    out=[]
    for key,rr in sorted(groups.items()):
        mass,snr,ratio,sm=key; db=np.array([x['delta_bic'] for x in rr]); rec=np.array([x['ratio_rec'] for x in rr])
        out.append({'mass_msun':mass,'snr_full_window':snr,'ratio_true_at_t0':ratio,'start_M':sm,'trials':len(rr),
                    'strong_fraction_dBIC_gt_10':float(np.mean(db>10)),'positive_fraction_dBIC_gt_0':float(np.mean(db>0)),
                    'median_delta_bic':float(np.median(db)),'p10_delta_bic':float(np.quantile(db,.1)),
                    'p90_delta_bic':float(np.quantile(db,.9)),'median_ratio_rec':float(np.median(rec))})
    return out

def threshold_table(summary,target=.90):
    out=[]; keys=sorted({(r['mass_msun'],r['snr_full_window'],r['start_M']) for r in summary})
    for mass,snr,sm in keys:
        arr=sorted([r for r in summary if r['mass_msun']==mass and r['snr_full_window']==snr and r['start_M']==sm],key=lambda x:x['ratio_true_at_t0'])
        det=[r for r in arr if r['ratio_true_at_t0']>0 and r['strong_fraction_dBIC_gt_10']>=target]
        ctrl=next(r for r in arr if r['ratio_true_at_t0']==0)
        out.append({'mass_msun':mass,'snr_full_window':snr,'start_M':sm,
                    'ratio_threshold_for_90pct_strong':det[0]['ratio_true_at_t0'] if det else float('nan'),
                    'zero_injection_false_positive_fraction':ctrl['strong_fraction_dBIC_gt_10'],
                    'zero_injection_median_delta_bic':ctrl['median_delta_bic']})
    return out

def write_csv(path,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--trials',type=int,default=100); ap.add_argument('--seed',type=int,default=1701); ap.add_argument('--fs',type=float,default=4096); ap.add_argument('--duration',type=float,default=.30); a=ap.parse_args()
    caches={m:build_mass_cache(m,a.fs,a.duration) for m in MASSES}; rng=np.random.default_rng(a.seed); rows=[]
    for mass in MASSES:
        for snr in SNRS:
            for ratio in RATIOS:
                for _ in range(a.trials): rows.extend(one_injection(rng,caches[mass],mass,snr,ratio))
    summary=summarize(rows); thresholds=threshold_table(summary); here=Path(__file__).resolve().parent
    write_csv(here/'detector_injection_trials.csv',rows); write_csv(here/'detector_injection_summary.csv',summary); write_csv(here/'detector_thresholds.csv',thresholds)
    with zipfile.ZipFile(here/'HBH_DETECTOR_INJECTION_RESULTS.zip','w',zipfile.ZIP_DEFLATED) as z:
        for fn in ['detector_injection_trials.csv','detector_injection_summary.csv','detector_thresholds.csv']: z.write(here/fn,fn)
if __name__=='__main__': main()

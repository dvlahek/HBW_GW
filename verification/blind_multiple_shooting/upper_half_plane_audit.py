"""Physical-sheet upper-half-plane winding audit for CQG-117441.

Region: -0.05 <= Re(M omega) <= 1, 0.001 <= Im(M omega) <= 0.5.
For Im(omega)>0 spectral.retarded_transform uses the direct physical sheet.
Bulk roots are followed by residual-strict predictor continuation. If a target
step leaves the local least-squares basin, the omega step is recursively
bisected. Only the requested contour points enter the winding sum.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from independent_solver import global_condition

HERE=Path(__file__).resolve().parent
BASE_SEEDS=(-0.06+0.02j,-0.15-0.02j)
FALLBACKS=((-0.04+0.05j,-0.18-0.05j),(-0.08+0.08j,-0.14-0.08j),(0+0j,-0.22+0j))
ANCHOR_RE=0.4; MAX_RELSV=1e-7; MIN_SEP=1e-5; MAX_BISECT=12

def contour(xmin,xmax,ymin,ymax,n):
    width,height=xmax-xmin,ymax-ymin; per=2*(width+height)
    nb=max(2,round(n*width/per)); nr=max(2,round(n*height/per)); nt=max(2,round(n*width/per)); nl=max(2,n-nb-nr-nt)
    pts=([complex(x,ymin) for x in np.linspace(xmin,xmax,nb,endpoint=False)]+[complex(xmax,y) for y in np.linspace(ymin,ymax,nr,endpoint=False)]+[complex(x,ymax) for x in np.linspace(xmax,xmin,nt,endpoint=False)]+[complex(xmin,y) for y in np.linspace(ymax,ymin,nl,endpoint=False)])
    i0=min(range(len(pts)),key=lambda i:abs(pts[i].imag-ymin)+abs(pts[i].real-ANCHOR_RE)); pts=pts[i0:]+pts[:i0]
    pts=[pts[0]]+list(reversed(pts[1:])); return pts+[pts[0]]

def pair_distance(a,b): return min(abs(a[0]-b[0])+abs(a[1]-b[1]),abs(a[0]-b[1])+abs(a[1]-b[0]))
def ordered_like(lams,target):
    return tuple(lams) if abs(lams[0]-target[0])+abs(lams[1]-target[1])<=abs(lams[0]-target[1])+abs(lams[1]-target[0]) else (lams[1],lams[0])
def predictor(target_w,prev_w,prev_lams,prevprev_w=None,prevprev_lams=None):
    if prevprev_w is None or prevprev_lams is None or abs(prev_w-prevprev_w)<1e-15:return prev_lams
    r=(target_w-prev_w)/(prev_w-prevprev_w); return (prev_lams[0]+r*(prev_lams[0]-prevprev_lams[0]),prev_lams[1]+r*(prev_lams[1]-prevprev_lams[1]))

def solve_point(w,previous,pred,rmax):
    attempts=[]
    for s in (pred,previous,*FALLBACKS):
        if s is not None and not any(pair_distance(s,x)<1e-12 for x in attempts):attempts.append(s)
    valid=[]
    for seeds in attempts:
        try:
            D,lams,relsvs,chis=global_condition(w,seeds,Rmax=rmax); target=pred if pred is not None else previous; lams=ordered_like(tuple(lams),target)
            sep=abs(lams[0]-lams[1]); maxsv=max(map(float,relsvs))
            if np.isfinite(D.real) and np.isfinite(D.imag) and sep>MIN_SEP and maxsv<MAX_RELSV: valid.append((pair_distance(lams,target),maxsv,D,lams,tuple(relsvs),tuple(chis),sep))
        except Exception: pass
    if not valid: raise RuntimeError(f"no residual-valid pair at w={w}")
    valid.sort(key=lambda x:(x[0],x[1])); _,_,D,lams,relsvs,chis,sep=valid[0]; return D,lams,relsvs,chis,sep

def advance(target_w,prev_w,prev_lams,prevprev_w,prevprev_lams,rmax,depth=0):
    pred=predictor(target_w,prev_w,prev_lams,prevprev_w,prevprev_lams)
    try:
        return solve_point(target_w,prev_lams,pred,rmax),(prev_w,prev_lams),0
    except RuntimeError:
        if depth>=MAX_BISECT: raise
        mid=prev_w+0.5*(target_w-prev_w)
        # Recursively reach the midpoint first. This is the key basin-preserving step.
        midres,midbridge,n1=advance(mid,prev_w,prev_lams,prevprev_w,prevprev_lams,rmax,depth+1)
        _,midlams,_,_,_=midres
        # Then continue from the accepted midpoint to the original target.
        finalres,finalbridge,n2=advance(target_w,mid,midlams,prev_w,prev_lams,rmax,depth+1)
        return finalres,finalbridge,n1+n2+1

def run(label,n,bounds,rmax,outdir):
    pts=contour(*bounds,n); vals=[]; rows=[]; adaptive_total=0
    w0=pts[0]; D,lams,relsvs,chis,sep=solve_point(w0,BASE_SEEDS,BASE_SEEDS,rmax)
    prevprev_w=prevprev_lams=None; prev_w=w0; prev_lams=lams
    def add(k,w,D,lams,relsvs,chis,sep,sub):
        vals.append(D); rows.append({"k":k,"w_re":w.real,"w_im":w.imag,"D_re":D.real,"D_im":D.imag,"abs_D":abs(D),"arg_D":float(np.angle(D)),"lambda1_re":lams[0].real,"lambda1_im":lams[0].imag,"lambda2_re":lams[1].real,"lambda2_im":lams[1].imag,"lambda_sep":sep,"rel_sv1":float(relsvs[0]),"rel_sv2":float(relsvs[1]),"adaptive_substeps":sub,"chiJ_re":chis[0].real,"chiJ_im":chis[0].imag,"chiX_re":chis[1].real,"chiX_im":chis[1].imag})
    add(0,w0,D,lams,relsvs,chis,sep,0)
    for k,w in enumerate(pts[1:],1):
        result,bridge,sub=advance(w,prev_w,prev_lams,prevprev_w,prevprev_lams,rmax); D,lams,relsvs,chis,sep=result; adaptive_total+=sub
        if sub>0: prevprev_w,prevprev_lams=bridge
        else: prevprev_w,prevprev_lams=prev_w,prev_lams
        prev_w,prev_lams=w,lams; add(k,w,D,lams,relsvs,chis,sep,sub)
        if k%25==0 or k==len(pts)-1: print(f"[{label}] {k}/{len(pts)-1} w={w.real:+.6f}{w.imag:+.6f}i |D|={abs(D):.6e} maxSV={max(relsvs):.3e} sep={sep:.3e} sub={sub}",flush=True)
    z=np.asarray(vals,complex); steps=np.angle(z[1:]/z[:-1]); un=np.unwrap(np.angle(z)); winding=float(np.sum(steps)/(2*np.pi))
    for r,p in zip(rows,un):r["arg_D_unwrapped"]=float(p)
    cp=outdir/f"upper_half_plane_{label}.csv"
    with cp.open("w",newline="",encoding="utf-8") as f:wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    aa=np.abs(z);i=int(np.argmin(aa));mx=float(np.max(np.abs(steps)))
    return {"label":label,"physical_sheet":True,"bounds":{"re_min":bounds[0],"re_max":bounds[1],"im_min":bounds[2],"im_max":bounds[3]},"requested_contour_points":n,"evaluated_points_including_closure":len(pts),"rmax":rmax,"winding":winding,"winding_integer":int(round(winding)),"winding_unwrapped":float((un[-1]-un[0])/(2*np.pi)),"min_abs_D":float(aa[i]),"min_abs_D_at":{"re":rows[i]["w_re"],"im":rows[i]["w_im"]},"max_abs_phase_step_rad":mx,"phase_step_below_pi":bool(mx<np.pi),"closure_relative_D_mismatch":float(abs(z[-1]-z[0])/max(abs(z[0]),1e-300)),"max_full_bulk_relative_sv":float(max(max(r["rel_sv1"],r["rel_sv2"]) for r in rows)),"min_bulk_root_separation":float(min(r["lambda_sep"] for r in rows)),"adaptive_substeps_total":adaptive_total,"csv":cp.name}

def main():
    p=argparse.ArgumentParser();p.add_argument("--n-primary",type=int,default=384);p.add_argument("--n-dense",type=int,default=768);p.add_argument("--re-min",type=float,default=-0.05);p.add_argument("--re-max",type=float,default=1.0);p.add_argument("--im-min",type=float,default=0.001);p.add_argument("--im-max",type=float,default=0.5);p.add_argument("--rmax",type=float,default=40.0);p.add_argument("--outdir",type=Path,default=HERE/"upper_half_plane_audit_output");a=p.parse_args()
    if a.im_min<=0 or not(a.re_min<0<a.re_max):raise ValueError("invalid contour")
    a.outdir.mkdir(parents=True,exist_ok=True);b=(a.re_min,a.re_max,a.im_min,a.im_max);ss=[run("primary",a.n_primary,b,a.rmax,a.outdir),run("dense",a.n_dense,b,a.rmax,a.outdir)]
    agree=ss[0]["winding_integer"]==ss[1]["winding_integer"];zero=all(s["winding_integer"]==0 for s in ss);phase=all(s["phase_step_below_pi"] for s in ss);res=all(s["max_full_bulk_relative_sv"]<MAX_RELSV for s in ss);close=all(s["closure_relative_D_mismatch"]<1e-8 for s in ss)
    final={"audit":"physical-sheet upper-half-plane growing-mode audit","interpretation_scope":"finite-window linear spectral audit only","bulk_residual_gate":MAX_RELSV,"contours":ss,"integer_windings_agree":agree,"both_windings_zero":zero,"phase_continuity_gate":phase,"bulk_residual_gate_pass":res,"closure_gate_pass":close,"verdict":"PASS_NO_GROWING_ROOTS_IN_TESTED_REGION" if agree and zero and phase and res and close else "REVIEW_REQUIRED"}
    (a.outdir/"upper_half_plane_audit_summary.json").write_text(json.dumps(final,indent=2),encoding="utf-8");lines=["PHYSICAL-SHEET UPPER-HALF-PLANE AUDIT","======================================",f"Region: {a.re_min} <= Re(Mw) <= {a.re_max}, {a.im_min} <= Im(Mw) <= {a.im_max}",f"Bulk residual acceptance gate: {MAX_RELSV:.1e}",""]
    for s in ss:lines += [f"{s['label']}:",f"  points: {s['requested_contour_points']}",f"  winding: {s['winding']:.12f} (integer {s['winding_integer']})",f"  min |D|: {s['min_abs_D']:.12e} at {s['min_abs_D_at']['re']:+.9f}{s['min_abs_D_at']['im']:+.9f}i",f"  max |phase step|: {s['max_abs_phase_step_rad']:.6f} rad",f"  closure relative D mismatch: {s['closure_relative_D_mismatch']:.3e}",f"  max full-bulk relative SV: {s['max_full_bulk_relative_sv']:.3e}",f"  min bulk-root separation: {s['min_bulk_root_separation']:.6e}",f"  adaptive substeps: {s['adaptive_substeps_total']}",""]
    lines.append(f"VERDICT: {final['verdict']}");(a.outdir/"upper_half_plane_audit_summary.txt").write_text("\n".join(lines)+"\n",encoding="utf-8");print("\n".join(lines),flush=True)
    if final["verdict"]!="PASS_NO_GROWING_ROOTS_IN_TESTED_REGION":raise SystemExit(2)
if __name__=="__main__":main()

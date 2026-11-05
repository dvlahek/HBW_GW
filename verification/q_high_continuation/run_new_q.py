#!/usr/bin/env python3
"""Fail-closed high-q Einstein-Vlasov spectral intake and pole/residue runner.

No new q>.20 stationary background or action-angle spectral catalogue
is available in the publication repository. This code DOES NOT synthesize
one by interpolating the q=.05-.20 pole table or reusing the q=.10 NPZ.

With genuinely NEW, independently generated, self-consistent background
and action-angle nodes for a specific q>.20, this runner:
  - checks independently recorded provenance and frozen SHA-256 bytes
  - checks spectral array schema and shell/harmonic structure
  - evaluates the original retarded kinetic susceptibilities
  - solves the complex global polar Einstein-Vlasov determinant
  - checks winding, root, derivative stability, source gravitational residues
  - saves inputs and numerical errors without publishing a q-continuation
    unless ALL physical and numerical gates pass.

The first upstream background *generator* is still missing; manifest
constraint metrics are user-supplied diagnostics, not independently
recomputed physics. No numerical result for q>.20 can be claimed
until a valid upstream generator/output bundle is supplied.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic8_landau_continuation as landau
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9
from publication_inputs import Q010_SHA256,Q_FAMILY_TABLE,SPECTRAL_KEYS,load_publication_case

ALLOWED_SOURCE_Q=(.25,.3,.4,.5,.6,.7,.8,.9,.95,.97,.99,.995,.999)
MAX_NPZ_BYTES=256*1024*1024
MAX_MASS_ERROR=1e-5
MAX_CONSTRAINT_NORM=1e-5
EXPECTED_MODEL="published_reference_smooth_kinetic_family"
NUMERICAL_ROOT_CUTOFF=2e-6
CONTINUATION_MAX_Q_STEP=.11
ALLOWED_HARMONICS_R={-3,-2,-1,1,2,3}
ALLOWED_HARMONICS_P={-2,0,2}


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def finite_scalar(obj,key):
    if key not in obj:
        raise ValueError("Missing required background physical diagnostic: "+key)
    v=float(obj[key])
    if not math.isfinite(v):
        raise ValueError("Nonfinite diagnostic: "+key)
    return v


def within_directory(parent,relative):
    if (not isinstance(relative,str) or not relative or
        Path(relative).is_absolute() or ".." in Path(relative).parts):
        raise ValueError("Manifest may only reference local in-directory relative input files")
    path=(parent/relative).resolve()
    if not path.is_relative_to(parent.resolve()) or not path.is_file():
        raise FileNotFoundError("Required physical input file is missing: "+str(path))
    if path.is_symlink():
        raise ValueError("Input may not be symlinked")
    return path


def manifest_signature_payload(j):
    return json.dumps({k:v for k,v in j.items() if k!="manifest_sha256"},
                      sort_keys=True,separators=(",",":"))


def validate_manifest(path):
    path=Path(path).resolve()
    j=json.loads(path.read_text())
    if j.get("schema_version")!=1:
        raise ValueError("Expected a versioned high-q physical spectral archive manifest")
    q=float(j.get("q",-1))
    if not any(abs(q-v)<1e-10 for v in ALLOWED_SOURCE_Q):
        raise ValueError("Only prelisted q>.20 new backgrounds accepted")
    if j.get("model_family")!=EXPECTED_MODEL:
        raise ValueError("Wrong smooth equilibrium physical family")
    if j.get("new_self_consistent_equilibrium_constructed") is not True:
        raise ValueError("Cannot substitute interpolation for NEW stationary Einstein-Vlasov solution")
    if j.get("q010_spectral_reused") is not False:
        raise ValueError("Reused q=.10 response kernels are scientifically invalid")
    # New C1-C7 reconstruction uncovered catastrophic cancellation in the
    # original energy-integrated C4 weights near q=1. A finite SHA-locked
    # spectral array alone does not establish a resolved source response.
    if j.get("spectral_convergence_verified") is not True:
        raise ValueError("High-q spectral quadrature convergence must be independently verified")
    err=float(j.get("spectral_convergence_max_relative",float("inf")))
    if not math.isfinite(err) or err<0 or err>.03:
        raise ValueError("High-q source response failed frozen 3% 61-vs-81 quadrature convergence gate")
    generator=j.get("generator",{})
    if (not isinstance(generator,dict) or
        not isinstance(generator.get("script_sha256"),str) or
        len(generator["script_sha256"])!=64 or
        not isinstance(generator.get("source_description"),str) or
        len(generator["source_description"].strip())<16 or
        not isinstance(generator.get("command"),str) or
        len(generator["command"].strip())<5):
        raise ValueError("Missing physical generator code identity and command")
    outputs={}
    for key in ("background","spectral_nodes"):
        meta=j.get(key,{})
        if not isinstance(meta,dict) or len(str(meta.get("sha256","")))!=64:
            raise ValueError("Missing exact SHA lock for "+key)
        source=within_directory(path.parent,meta.get("file"))
        nbytes=source.stat().st_size
        if key=="spectral_nodes" and not 1000<nbytes<MAX_NPZ_BYTES:
            raise ValueError("Spectral nodes file size not credible/safe")
        if sha256(source)!=meta["sha256"]:
            raise ValueError("Physical input checksum does not match frozen manifest: "+key)
        outputs[key]=source
    if j["spectral_nodes"]["sha256"]==Q010_SHA256:
        raise ValueError("q>.20 cannot reuse published q=.10 spectral nodes")
    b=json.loads(outputs["background"].read_text())
    if not isinstance(b,dict):
        raise ValueError("Bad independent stationary background audit")
    if abs(finite_scalar(b,"q_integrated_from_background")-q)>1e-8:
        raise ValueError("Derived total gravitational mass ratio disagrees with target q")
    for k in ("relative_mass_balance_error","einstein_constraint_max",
              "vlasov_stationarity_max","stress_energy_conservation_max"):
        x=abs(finite_scalar(b,k))
        cap=MAX_MASS_ERROR if k=="relative_mass_balance_error" else MAX_CONSTRAINT_NORM
        if x>cap:
            raise ValueError(f"Physical background {k} exceeds gate {cap}")
    N=int(finite_scalar(b,"shell_count"))
    if N<8 or N>10000:
        raise ValueError("Invalid number of independent action-angle shells")
    for k in ("positive_phase_space_density","physical_orbit_support",
              "stress_energy_from_same_distribution"):
        if b.get(k) is not True:
            raise ValueError("Missing physical assertion "+k)
    if abs(finite_scalar(b,"ocean_length_M")-9*math.log(1/q))>1e-6:
        raise ValueError("Background has inconsistent q/coordinate ocean length")
    return q,j,b,outputs


def load_spectral_arrays(npz_path,shell_count):
    if Path(npz_path).suffix.lower()!=".npz":
        raise ValueError("Expected unpickled compressed NumPy spectral nodes")
    with np.load(npz_path,allow_pickle=False) as data:
        if set(data.files)!=SPECTRAL_KEYS:
            raise ValueError("Missing/extra kinetic branch fields")
        arrays={key:np.asarray(data[key]) for key in SPECTRAL_KEYS}
    for key in SPECTRAL_KEYS:
        if arrays[key].ndim!=1 or arrays[key].size==0:
            raise ValueError("Kinetic node field must be nonempty 1D: "+key)
        if not np.issubdtype(arrays[key].dtype,np.number):
            raise ValueError("Kinetic node field not numeric: "+key)
        if not np.all(np.isfinite(arrays[key])):
            raise ValueError("Kinetic node field nonfinite: "+key)
    rkeys=("nu_J","q_JJ","q_JX","n_r","n_phi","shell_J")
    if len(set(arrays[k].shape for k in rkeys))!=1:
        raise ValueError("Radial kinetic node arrays have incompatible shapes")
    hkeys=("nu_H","q_HH","shell_H")
    if len(set(arrays[k].shape for k in hkeys))!=1:
        raise ValueError("Carrier kinetic node arrays have incompatible shapes")
    if arrays["nu_J"].size<10*shell_count or arrays["nu_H"].size<4*shell_count:
        raise ValueError("Too few kinetic integration quadrature nodes")
    for name in ("shell_J","shell_H","n_r","n_phi"):
        values=arrays[name]
        if not np.all(values==np.rint(values)):
            raise ValueError("Discrete shell/harmonic index is fractional: "+name)
    ns=set(np.asarray(arrays["shell_J"],int))
    nh=set(np.asarray(arrays["shell_H"],int))
    if len(ns)!=shell_count or len(nh)!=shell_count or ns!=nh:
        raise ValueError("Discrepant spectral shell catalogues or shell count")
    if not set(np.asarray(arrays["n_r"],int)).issubset(ALLOWED_HARMONICS_R):
        raise ValueError("Nonphysical radial harmonics in published truncation")
    if not set(np.asarray(arrays["n_phi"],int)).issubset(ALLOWED_HARMONICS_P):
        raise ValueError("Nonphysical azimuthal harmonics in published truncation")
    return arrays


def inspect_reference_q010():
    radial,carrier,row=load_publication_case(.10)
    spectral=ROOT/"verification/blind_multiple_shooting/blind_q010_spectral_nodes.npz"
    N=int(row["N"])
    arrays=load_spectral_arrays(spectral,N)
    if len(carrier)!=80 or len(radial)!=1440 or len(arrays["nu_J"])<1000:
        raise RuntimeError("Original public q=.10 spectral catalogue failed dimensional smoke gate")
    return {"q":.1,"N":N,"sha256":sha256(spectral),
            "radial_branches":len(radial),"carrier_branches":len(carrier),
            "radial_node_count":int(len(arrays["nu_J"])),
            "carrier_node_count":int(len(arrays["nu_H"])),
            "status":"PUBLIC Q=.10 REFERENCE ONLY, NOT q>.20 data"}


def inspect_new_q(path):
    q,manifest,b,outputs=validate_manifest(path)
    arrays=load_spectral_arrays(outputs["spectral_nodes"],int(b["shell_count"]))
    radial,carrier=landau.build_catalog(arrays)
    if len(carrier)!=int(b["shell_count"]) or len(radial)<int(b["shell_count"]):
        raise ValueError("Retarded kinetic catalogue does not cover independent shells")
    return {"q":q,"manifest":manifest,
            "background_asserted_diagnostics":b,
            "spectral_SHA256":manifest["spectral_nodes"]["sha256"],
            "spectral_path":str(outputs["spectral_nodes"]),
            "radial_branch_count":len(radial),"carrier_branch_count":len(carrier),
            "radial_nodes":int(len(arrays["nu_J"])),
            "carrier_nodes":int(len(arrays["nu_H"])),
            "radial":radial,"carrier":carrier}


def seed_from_predecessor(q,path):
    """Mandatory adjacent q publication root/accepted prior-q continuation."""
    if path is None:
        if not (.2<q<=.3+1e-10):
            raise ValueError("First new q must start near published q=.20; chain later q values with --previous-result")
        with Q_FAMILY_TABLE.open(newline="") as f:
            rows=list(csv.DictReader(f))
        row=next(r for r in rows if abs(float(r["q"])-.2)<1e-10)
        return (.2,complex(float(row["root_re"]),float(row["root_im"])),
                (complex(float(row["lambda1_re"]),float(row["lambda1_im"])),
                 complex(float(row["lambda2_re"]),float(row["lambda2_im"]))))
    doc=json.loads(Path(path).read_text())
    qold=float(doc["q"])
    if (doc.get("status")!="DOWNSTREAM_NUMERICAL_PASS_UPSTREAM_UNVERIFIED" or
        not .2<qold<q or q-qold>CONTINUATION_MAX_Q_STEP or
        doc.get("source_input_sha_and_declared_physics_passed") is not True):
        raise ValueError("Predecessor must be a close ACCEPTED q>.20 independent physical continuation")
    w=complex(*doc["root_Momega"])
    lams=tuple(complex(*x) for x in doc["lambda_pair"])
    return qold,w,lams


def run_solver(case,previous_result,outdir):
    q=case["q"]
    previous_q,seed,lams=seed_from_predecessor(q,previous_result)
    # New q=.30 node spectra can cross dozens of real-frequency spline
    # knots over the steps previously used for dW/domega. The original
    # lower-half C8 interpolation is nonholomorphic across these knots:
    # block a falsely precise QNM or gravitational residue at the outset.
    from retarded_analyticity_gate import pack_catalog,derivative_analyticity_gate
    spectral_gate=derivative_analyticity_gate(
        pack_catalog(case["radial"],case["carrier"]),seed)
    if not spectral_gate["global_QNM_residue_gate_pass"]:
        target=Path(outdir)
        target.mkdir(parents=True,exist_ok=True)
        (target/f"q_{q:.3f}_BLOCKED_nonanalytic_retarded.json").write_text(
            json.dumps({
                "status":"BLOCKED_NONANALYTIC_RETARDED_C8",
                "q":q,
                "new_finite_q_root_computed":False,
                "new_gravitational_residue_computed":False,
                "physical_input_SHA256":case["spectral_SHA256"],
                "diagnostics":spectral_gate,
                "required_fix":"Validated globally analytic retarded kinetic continuation, with q>.2 quadrature and dW/domega convergence"
            },indent=2)+"\n")
        raise RuntimeError("Physical q data valid but retarded C8 is nonanalytic across spectral knots: global mode/residue claim BLOCKED")
    xivc=c9.load_xivc()
    phase_ix=xivc.load_phase_ix()
    Mfun,row_names=c9.build_generic_bulk_matrix(phase_ix)
    c9.ACTIVE_LAMBDA_SEED_PAIR=lams
    det=c9.SmoothGlobalDet(xivc,Mfun,row_names,case["radial"],case["carrier"],q=q)
    root=c9.solve_global(det,seed,max_nfev=160)
    ad,dr,w,ans,info=root
    if not np.isfinite(ad) or ad>NUMERICAL_ROOT_CUTOFF:
        raise RuntimeError(f"No accepted global root at q={q} (|D|={ad})")
    if abs(w-seed)>.05:
        raise RuntimeError("Unexpected pole jump; manual branch-identity audit required")
    wind=c9.winding_number(det,w,radius=4e-4,n=40)
    if abs(wind-1.)>.1:
        raise RuntimeError("Global root cannot be certified as simple isolated pole by local winding")
    d1=c9.finite_difference_derivative(det,w,2e-5)
    d2=c9.finite_difference_derivative(det,w,1e-5)
    drel=float(abs(d1-d2)/max(abs(d2),1e-300))
    if drel>.025 or abs(d2)<1e-7:
        raise RuntimeError(f"Unstable W derivative or degenerate root: {drel}")
    if info is None or not np.all(np.isfinite(info["lambdas"])):
        raise RuntimeError("No validated kinetic lambda pair at new q")
    import hbh_np_source_excitation_test_v2 as source
    summary,point,smooth=source.run_case(
        q,case["radial"],case["carrier"],w,
        tuple(complex(x) for x in info["lambdas"]),
        xivc=xivc,Mfun=Mfun,row_names=row_names)
    if summary["dW_derivative_relative_difference"]>.025:
        raise RuntimeError("Unstable unnormalized physical source Green-function residue")
    if not all(np.isfinite(x["residue_abs"]) and x["residue_abs"]>0
               for x in point+smooth):
        raise RuntimeError("Nonfinite source-normalized physical Green-function residues")
    result={
        "status":"DOWNSTREAM_NUMERICAL_PASS_UPSTREAM_UNVERIFIED",
        "q":q,"previous_q":previous_q,
        "source_input_sha_and_declared_physics_passed":True,
        "upstream_equilibrium_constraints_independently_recomputed":False,
        "input_manifest":case["manifest"],
        "source_spectral_SHA256":case["spectral_SHA256"],
        "root_Momega":[float(w.real),float(w.imag)],
        "lambda_pair":[[float(complex(x).real),float(complex(x).imag)]
                       for x in info["lambdas"]],
        "normalized_global_det_abs":float(ad),
        "root_step_from_seed_abs":float(dr),
        "winding":float(wind),
        "normalized_derivative_h2e5":[float(d1.real),float(d1.imag)],
        "normalized_derivative_h1e5":[float(d2.real),float(d2.imag)],
        "normalized_derivative_relative_difference":drel,
        "source_summary":summary,"point_source_residues":point,
        "smooth_source_residues":smooth,
        "caveats":[
            "Upstream generator's physical Einstein-Vlasov constraints are declared by source manifest; this code does not independently rederive upstream equilibrium",
            "This single q candidate is not a completed q->1 branch continuation",
            "A pole and external-source residue are not a predicted astrophysical merger excitation or detector strain",
            "Global winding here is a numerically sampled Riccati-type determinant and should receive an independent full-state root-count check before publication",
            "Residues must be independently normalized against same vacuum source before reporting detection-related amplitude ratios"
        ]}
    dest=Path(outdir);dest.mkdir(parents=True,exist_ok=True)
    (dest/f"q_{q:.3f}_numerical_result.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--selftest-q010",action="store_true")
    parser.add_argument("--manifest",type=Path)
    parser.add_argument("--validate-only",action="store_true")
    parser.add_argument("--solve",action="store_true")
    parser.add_argument("--previous-result",type=Path)
    parser.add_argument("--output-dir",type=Path,default=ROOT/"results/q_high_continuation")
    args=parser.parse_args()
    if args.selftest_q010:
        if args.manifest or args.solve:
            parser.error("--selftest-q010 does not generate or solve a new physical q")
        x=inspect_reference_q010()
        print("Q_HIGH_REFERENCE_Q010_INPUT_GATE_PASS",json.dumps(x,separators=(",",":")))
        return
    if args.manifest is None:
        parser.error("NEW q>.20 stationary-background + action-angle spectral manifest is REQUIRED; cannot interpolate archived roots")
    if args.solve == args.validate_only:
        parser.error("Choose exactly one of --validate-only or --solve")
    case=inspect_new_q(args.manifest)
    public={k:v for k,v in case.items() if k not in ("radial","carrier")}
    if args.validate_only:
        print("Q_HIGH_PHYSICAL_INPUT_INTAKE_PASS",json.dumps(public,separators=(",",":")))
        return
    result=run_solver(case,args.previous_result,args.output_dir)
    print("Q_HIGH_NEW_PHYSICAL_POLE_RESIDUE_PASS",json.dumps({
        "q":result["q"],"root":result["root_Momega"],
        "winding":result["winding"],
        "source_spectral_SHA256":result["source_spectral_SHA256"]
    },separators=(",",":")))


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Re-fit both frozen crossed-origin CVs and calculate paired origin contrasts.

Post-result diagnostic: water geometry is undirected; WPI port-count product is
not vessel OD. Matched predictive gains do not identify transport mechanisms.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math,subprocess,sys,zipfile
from pathlib import Path
import numpy as np

ZIP_SHA={
 14296307100:"c00508ea0384efe73841a82a7ad5e73d3df8a13eed36ecea23b4546238d0d2b4",
 14296891019:"36dbee9c52bebf2c8110fd894bd915f11b748c45057c3576dbd67df65bb8a8ff",
 14296766522:"6dcd3987d79f117a782ea919c6961f2ede8fac7d93af9228ffa62fb9a072186b",
 14296567644:"3380b9ecca7f2b399e9a5da2856d8d4ec9753b476ab8aec40d858ed13f353780",
}
P1="scripts/audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117.py"
P2="scripts/audit_WCVP126_native_alien_double_unseen_islands_20261117.py"
AREA="results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv"
BLOBS={P1:"0c06c1f013d58624b83c1b03da8e433f77e3e409",
       P2:"afd22b2cad2243da45a56be21ccf59467461a1d1",
       AREA:"53b618fa51a1a418ccc6e4363f4c866c389c90cc"}
F1="WCVP126_native_alien_crossed_held_island_prediction_audit.json"
F2="original_WCVP126_native_alien_double_unseen_5fold.json"
MODELS=("B","G","P","GP")
METRICS=("G_vs_B_pct","P_vs_B_pct","G_added_after_P_pct","P_added_after_G_pct")
CONTRASTS=("native_minus_introduced_G_vs_B","introduced_minus_native_P_vs_B",
           "native_minus_introduced_G_after_P","introduced_minus_native_P_after_G")
def need(ok,why):
    if not ok:raise RuntimeError(why)
def digest(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def compare(old,new,path="",state=None):
    if state is None:state={"leaves":0,"numeric_leaves":0,"max_abs_numeric_difference":0.}
    if isinstance(old,dict):
        need(isinstance(new,dict) and set(old)==set(new),"keys "+path)
        for k in old:compare(old[k],new[k],path+"/"+k,state)
    elif isinstance(old,list):
        need(isinstance(new,list) and len(old)==len(new),"list "+path)
        for i,(a,b) in enumerate(zip(old,new)):compare(a,b,path+"/"+str(i),state)
    elif isinstance(old,(int,float)) and not isinstance(old,bool):
        need(isinstance(new,(int,float)) and not isinstance(new,bool) and math.isclose(old,new,rel_tol=1e-8,abs_tol=1e-10),f"numeric {path}: {old} vs {new}")
        state["numeric_leaves"]+=1
        state["max_abs_numeric_difference"]=max(state["max_abs_numeric_difference"],abs(old-new))
        state["leaves"]+=1
    else:
        need(old==new,f"exact {path}: {old} vs {new}")
        state["leaves"]+=1
    return state
def gains(x):
    b,g,p,gp=np.moveaxis(x,-1,0)
    return np.stack([100*(b-g)/b,100*(b-p)/b,100*(p-gp)/p,100*(g-gp)/g],axis=-1)
def direct(d):
    return np.stack([d[...,0,0]-d[...,1,0],d[...,1,1]-d[...,0,1],
                     d[...,0,2]-d[...,1,2],d[...,1,3]-d[...,0,3]],axis=-1)
def paired(report,strict):
    k="folds" if strict else "per_island_loss"
    origin="alien_flora" if strict else "introduced_flora"
    key="fold" if strict else "code"
    native=report["native_flora"][k]; introduced=report[origin][k]
    ids=[x[key] for x in native]
    need(ids==[x[key] for x in introduced] and len(ids)==len(set(ids))==(5 if strict else 104),"paired fold roster")
    for x,y in zip(native,introduced):
        if strict:
            for field in ("held_island_codes","n_unseen_islands","n_train_pairs","n_test_double_unseen_pairs"):
                need(x[field]==y[field],"paired strict cohort "+field)
        else:need(x["n_test"]==y["n_test"],"paired incident count")
    loss=np.asarray([[[r["MSE_"+m] for m in MODELS] for r in rows] for rows in [native,introduced]]).transpose(1,0,2)
    need(np.isfinite(loss).all() and (loss>0).all(),"loss validity")
    n=len(ids)
    point=gains(loss.mean(0))
    field="macro_equal_fold_relative_mse_improvement_percent" if strict else "relative_delta_MSE_percent_positive_is_improvement"
    for i,lab in enumerate(("native_flora",origin)):
        expected=[report[lab][field][m] for m in METRICS]
        need(np.allclose(point[i],expected,rtol=1e-8,atol=1e-10),"archived point gain mismatch")
    idx=np.asarray(list(itertools.product(range(n),repeat=n)),dtype=int) if strict else np.random.default_rng(20261125).integers(0,n,(9999,n))
    # Same draw indices for BOTH origins; each gain retains its own MSE denominator.
    boot=direct(gains(loss[idx].mean(1)))
    observed=direct(point)
    loo=np.array([direct(gains(np.delete(loss,i,axis=0).mean(0))) for i in range(n)])
    intervals=np.quantile(boot,[.025,.975],axis=0).T
    return {"paired_units":n,"unit":"five held-out fold groups" if strict else "104 held-out islands",
      "resamples":len(idx),"method":"all 5^5 ordered fold resamples" if strict else "9999 same-island draws; RNG seed 20261125",
      "native_gains_pct":dict(zip(METRICS,map(float,point[0]))),
      "introduced_gains_pct":dict(zip(METRICS,map(float,point[1]))),
      "direct_contrasts":{name:{"estimate_percentage_points":float(observed[i]),
        "descriptive_q025_q975":[float(x) for x in intervals[i]],
        "fraction_resamples_positive":float(np.mean(boot[:,i]>0)),
        "leave_one_unit_out_range":[float(loo[:,i].min()),float(loo[:,i].max())]}
        for i,name in enumerate(CONTRASTS)},
      "scope":"Descriptive paired prediction-loss resampling; dyads share endpoints and training folds overlap. Not independent causal inference or a calibrated hypothesis-test p-value."}
def main():
    p=argparse.ArgumentParser()
    for k in ("downloads","repo","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    for name,h in BLOBS.items():need(blob((a.repo/name).read_bytes())==h,"frozen source git blob "+name)
    unpack=a.out/"_source";unpack.mkdir()
    for ident,h in ZIP_SHA.items():
        path=a.downloads/f"{ident}.zip"
        need(digest(path.read_bytes())==h,"artifact digest "+str(ident))
        with zipfile.ZipFile(path) as z:
            need(z.testzip() is None,"CRC")
            for n in z.namelist():
                need(Path(n).name==n and not (unpack/n).exists(),"unexpected or duplicate member "+n)
                (unpack/n).write_bytes(z.read(n))
    flags=["--plants",str(unpack/"WCVP126_flora_similarity_strict_WPI2020_geographic_pairs.csv"),
           "--geometry",str(unpack/"WCVP126_GSHHG_source_coastal_anchors_and_sea_geometry_pairs.csv"),
           "--island-area",str(a.repo/AREA),
           "--sea-audit",str(unpack/"GSHHG126_sea_geometry_source_coverage_audit.json"),
           "--shipping-audit",str(unpack/"WCVP126_geodesic_port_shipment_source_coverage_audit.json"),
           "--out-dir",str(a.out)]
    for name in (P1,P2):subprocess.run([sys.executable,str(a.repo/name),*flags],check=True)
    one=json.loads((a.out/F1).read_text());two=json.loads((a.out/F2).read_text())
    need(one["n_common_anchored_source_islands"]==104 and one["n_common_valid_island_pairs"]==5355,"original matched source population")
    need(two["source_common_pairs"]==5355 and two["native_flora"]["n_double_unseen_test_pairs"]==1089,"strict sample")
    verified={}
    results={}
    for lab,name,doc,strict in (("single_held_island",F1,one,False),("double_unseen",F2,two,True)):
        original=json.loads((unpack/name).read_text())
        verified[lab]=compare(original,doc)
        fresh=paired(doc,strict)
        compare(paired(original,strict),fresh)
        results[lab]=fresh
    receipt={"schema":"chocho_WCVP104_reciprocal_origin_refit_and_paired_contrasts_v1",
      "status":"EXACT_SOURCE_REFIT_AND_PAIRED_DESCRIPTIVE_CONTRASTS",
      "source_artifact_sha256":ZIP_SHA,"original_code_and_area_git_blobs":BLOBS,
      "matched_islands":104,"matched_dyads":5355,"strict_double_unseen_test_dyads":1089,
      "source_report_reproduction":verified,"paired_predictive_origin_contrasts":results,
      "limits":["Same-origin relative MSE percent, then contrast in percentage points; distinct target richness margins and error scales.",
        "No causal double dissociation established by separate significance or prediction performance.",
        "Water predictor is undirected GSHHG shortest-path detour, not current direction; WPI is port-count product, not cargo or ship routes.",
        "Post-result diagnostic; no new source or outcomes selected. No raw WCVP or underlying coastline-graph re-acquisition."]}
    (a.out/"RECIPROCAL_ORIGIN_PAIRED_REPLAY.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    # No large input tables reuploaded; upstream source archives remain pinned.
    import shutil
    shutil.rmtree(unpack)
    print("RECIPROCAL_ORIGIN_REPLAY_PASS",json.dumps(receipt),flush=True)
if __name__=="__main__":main()

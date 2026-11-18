#!/usr/bin/env python3
"""Fixed 104-island WCVP whole-introduced-flora vs potential-butterfly-resource bridge.

Retrospective conditional prediction, NOT independent consumer biology. Use
original WCVP/HOSTS 239×369 source matrix and independently generated WCVP
plant-pair Jaccard. Five pre-existing source-hash island folds have no overlap
between held test island endpoints and training island endpoints.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 as src
import audit_WCVP126_native_alien_double_unseen_islands_20261117 as strict

SOURCE_SHA="ec6a5530d65ac0b04f230e0c51fd68bf543555641a17b94795996bf11ce257d0"
REPORT_SHA="ccda0bc669fb5ffe7326175fcc0501dbb1447790925c06add92501d1d30f255b"
N_REGIONS=104
N_SOURCE_PAIRS=5355
N_SOURCE_TEST_PAIRS=1089
N_DRAW=4999
SEED=20261118
ALPHA=20.0

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def csvrows(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
def jac(a,b):
    u=(a|b).bit_count()
    return (a&b).bit_count()/u if u else None
def tied_ranks(x):
    x=np.asarray(x,dtype=float)
    order=np.argsort(x,kind="mergesort");rank=np.empty(len(x),float)
    k=0
    while k<len(x):
        z=k+1
        while z<len(x) and x[order[z]]==x[order[k]]:
            z+=1
        rank[order[k:z]]=(k+z-1)/2.
        k=z
    return rank
def corr_spearman(x,y):
    a,b=tied_ranks(x),tied_ranks(y)
    return float(np.corrcoef(a,b)[0,1])
def setup(source_root,plants,geometry,area):
    srcdir=source_root/"provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031"
    matrix=srcdir/"beta_baseline_v3.matrix.csv.gz"
    receipt=srcdir/"beta_baseline_v3.json"
    if sha(matrix)!=SOURCE_SHA or sha(receipt)!=REPORT_SHA:
        raise RuntimeError("Original 239×369 source binary matrix or receipt hash changed")
    sys.path.insert(0,str(source_root/"scripts"))
    import analyze_global_resource_gain_fixed_margins as frozen
    regions,taxa,native,added,current,units=frozen.read_matrix(matrix,receipt)
    if len(regions)!=369 or len(taxa)!=239 or units!={"native":26530,"added":14553,"current":41083}:
        raise RuntimeError("Source resource counts/identity unexpectedly changed")
    if any((n&a) or (n|a)!=c for n,a,c in zip(native,added,current)):
        raise RuntimeError("Original resource nesting violated")
    region_idx={code:i for i,code in enumerate(regions)}
    areas={r["code"]:float(r["area_km2"]) for r in csvrows(area)}
    if len(areas)!=126:raise RuntimeError("Original 126 WCVP area roster changed")
    entries,missing=src.input_pair(SimpleNamespace(plants=plants,geometry=geometry))
    if len(entries)!=N_SOURCE_PAIRS or len(missing)!=1:
        raise RuntimeError("Original matched WCVP source pair cohort changed")
    codes={e["island1"] for e in entries}|{e["island2"] for e in entries}
    if len(codes)!=N_REGIONS or not codes.issubset(region_idx):
        raise RuntimeError("Original 104 coastal WCVP source identities mismatch 239×369 matrix")
    mapped=[]
    n_undefined_current=0
    n_undefined_native=0
    for e in entries:
        i,j=region_idx[e["island1"]],region_idx[e["island2"]]
        nres=jac(native[i],native[j])
        cres=jac(current[i],current[j])
        if nres is None:n_undefined_native+=1
        if cres is None:
            n_undefined_current+=1
            continue
        mapped.append({
            **e,"response":cres,
            "native_resource_jaccard":nres if nres is not None else 0.,
            "native_resource_missing":int(nres is None),
            "native_count_A":native[i].bit_count(),
            "native_count_B":native[j].bit_count(),
            "current_count_A":current[i].bit_count(),
            "current_count_B":current[j].bit_count(),
            "alien_flora_Jaccard":e["alien"],
            "native_flora_Jaccard":e["native"],
        })
    return mapped,areas,{
        "source_pairs":len(entries),"source_islands":len(codes),
        "native_resource_undefined_pairs":n_undefined_native,
        "current_resource_undefined_pairs":n_undefined_current,
        "eligible_source_pair_count":len(mapped),
        "native_zero_resources_islands":[r for r in sorted(codes) if native[region_idx[r]]==0],
        "current_zero_resources_islands":[r for r in sorted(codes) if current[region_idx[r]]==0],
        "source_taxa":len(taxa),"source_units":units
    }

def features(entries,area,*,include_current_margins,include_alien_similarity):
    arr=[]
    for z in entries:
        na,nb=z["native_count_A"],z["native_count_B"]
        ca,cb=z["current_count_A"],z["current_count_B"]
        a,b=area[z["island1"]],area[z["island2"]]
        row=[
          z["log_geo"],z["log_geo"]**2,z["same_level1"],z["abs_lat"],
          math.log1p(min(a,b)),math.log1p(max(a,b)),
          z["port_pair_proxy"],
          z["native_flora_Jaccard"],z["native_resource_jaccard"],
          float(z["native_resource_missing"]),
          math.log1p(min(na,nb)),math.log1p(max(na,nb)),
          math.log1p(min(z["n1"],z["n2"])),math.log1p(max(z["n1"],z["n2"])),
          math.log1p(min(z["a1"],z["a2"])),math.log1p(max(z["a1"],z["a2"])),
        ]
        if include_current_margins:
            row += [math.log1p(min(ca,cb)),math.log1p(max(ca,cb))]
        if include_alien_similarity:
            row += [z["alien_flora_Jaccard"]]
        arr.append(row)
    out=np.array(arr,dtype=float)
    if not np.isfinite(out).all():
        raise RuntimeError("Nonfinite plant/resource pair predictor")
    return out

def fold_cohort(entries,*,native_only):
    group={}
    for e in entries:
        for v in [e["island1"],e["island2"]]:
            group.setdefault(strict.label(v),set()).add(v)
    if len(group)!=5:raise RuntimeError("Original five source-hash island bins changed")
    idx_byfold=[]
    for fold in range(5):
        held=group[fold]
        tr=np.array([k for k,e in enumerate(entries)
             if e["island1"] not in held and e["island2"] not in held and
             (not native_only or not e["native_resource_missing"])],int)
        te=np.array([k for k,e in enumerate(entries)
             if e["island1"] in held and e["island2"] in held and
             (not native_only or not e["native_resource_missing"])],int)
        if len(tr)<1500 or len(te)<35:
            raise RuntimeError(f"Insufficient frozen fold support: {fold}, {len(tr)}, {len(te)}")
        # Exclude cross-boundary pairs for each fold, exactly as original source.
        assert all(entries[k]["island1"] not in held and entries[k]["island2"] not in held for k in tr)
        assert all(entries[k]["island1"] in held and entries[k]["island2"] in held for k in te)
        idx_byfold.append((fold,tr,te))
    return idx_byfold

def evaluate(entries,area,*,native_only=False):
    pairs=fold_cohort(entries,native_only=native_only)
    y=np.array([e["response"] for e in entries],dtype=float)
    fits={}
    for label,has_current in [("E",False),("S",True)]:
        X0=features(entries,area,include_current_margins=has_current,include_alien_similarity=False)
        X1=features(entries,area,include_current_margins=has_current,include_alien_similarity=True)
        allfolds=[]
        for fold,tr,te in pairs:
            pred0=src.ridge(X0[tr],y[tr],X0[te])
            pred1=src.ridge(X1[tr],y[tr],X1[te])
            mse0=float(np.mean((y[te]-pred0)**2))
            mse1=float(np.mean((y[te]-pred1)**2))
            allfolds.append({"fold":fold,"n_train_pairs":len(tr),"n_test_pairs":len(te),
                "MSE_baseline":mse0,"MSE_plus_alien_flora_jaccard":mse1,
                "percent_MSE_gain":100*(mse0-mse1)/mse0,
                "positive_target_fraction":float(np.mean(y[te]>0))})
        macro0=float(np.mean([r["MSE_baseline"] for r in allfolds]))
        macro1=float(np.mean([r["MSE_plus_alien_flora_jaccard"] for r in allfolds]))
        gain=100*(macro0-macro1)/macro0
        ww=np.array([r["n_test_pairs"] for r in allfolds])
        w0=np.average([r["MSE_baseline"] for r in allfolds],weights=ww)
        w1=np.average([r["MSE_plus_alien_flora_jaccard"] for r in allfolds],weights=ww)
        rng=np.random.default_rng(SEED+(0 if label=="E" else 1)+(9 if native_only else 0))
        mat=np.array([[z["MSE_baseline"],z["MSE_plus_alien_flora_jaccard"]] for z in allfolds])
        bootstrap=[]
        for z in range(N_DRAW):
            ix=rng.integers(0,len(mat),len(mat))
            m0,m1=mat[ix].mean(axis=0)
            bootstrap.append(100*(m0-m1)/m0)
        fits[label]={
            "baseline_features":X0.shape[1],
            "alien_jaccard_added_features":X1.shape[1]-X0.shape[1],
            "n_test_pairs":int(ww.sum()),
            "mean_source_y":float(np.mean(y)),
            "test_y_zero_fraction":float(np.mean(y[[k for _,_,te in pairs for k in te]]==0)),
            "macro_equal_fold_MSE_baseline":macro0,
            "macro_equal_fold_MSE_plus_alien_flora_jaccard":macro1,
            "macro_equal_fold_relative_gain_pct":gain,
            "pair_weighted_relative_gain_pct":float(100*(w0-w1)/w0),
            "n_folds_adding_alien_flora_helped":sum(r["percent_MSE_gain"]>0 for r in allfolds),
            "descriptive_fold_bootstrap_95":list(map(float,np.quantile(bootstrap,[.025,.975]))),
            "folds":allfolds,
        }
    return fits

def main():
    p=argparse.ArgumentParser()
    for name in ["archived-root","plants","geometry","island-area","out-dir"]:
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    if abs(src.RIDGE_ALPHA-ALPHA)>1e-12:
        raise RuntimeError("Original fixed ridge penalty drift")
    e,area,source=setup(a.archived_root,a.plants,a.geometry,a.island_area)
    # Exact original 104 island, 5,355 WCVP pairs and original 1,089 held source pair IDs.
    source_mask=np.array([strict.label(z["island1"])==strict.label(z["island2"]) for z in e])
    if int(source_mask.sum())+source["current_resource_undefined_pairs"] < N_SOURCE_TEST_PAIRS-3:
        raise RuntimeError("Unexpectedly lost source strict WCVP pair support")
    scores=evaluate(e,area,native_only=False)
    scores_native=evaluate(e,area,native_only=True)
    summary={
      "schema":"chocho_source_locked_104_WCVP_island_whole_alien_flora_vs_potential_239_butterfly_resource_similarity_v1",
      "status":"EXPLORATORY_OBSERVATIONAL_SHARED_WCVP_SOURCE_NO_BUTTERFLY_OCCUPANCY",
      "protocol":"results/island104_alien_flora_to_butterfly_resources_20261118/PROTOCOL.md",
      "source_original_matrix_commit":"bb4d8f00fbdb97fa2d90f606a8c9e7e0620646f2",
      "source_original_matrix_sha256":SOURCE_SHA,
      "source_original_WCVP_pair_CI":44010307551,
      "source_original_GSHHG_coast_CI":44010673580,
      "source_original_bit_counts":source,
      "source_original_ridge_alpha":ALPHA,
      "source_hashed_fivefold_label":"chocho-WCVP126-double-unseen-v1|TDWG3 code modulo5",
      "source_valid_paired_whole_flora_and_resource_pairs":len(e),
      "source_binary_bit_resource_current_Jaccard_mean_on_pairs":float(np.mean([z["response"] for z in e])),
      "raw_correlations_nonindependent":{
        "Spearman_alien_flora_jaccard_to_current_butterfly_potential_jaccard":corr_spearman(
            [z["alien_flora_Jaccard"] for z in e],[z["response"] for z in e]),
        "Spearman_native_flora_jaccard_to_current_butterfly_potential_jaccard":corr_spearman(
            [z["native_flora_Jaccard"] for z in e],[z["response"] for z in e]),
      },
      "all_source_pairs":{"no_current_resource_margin_baseline_E":scores["E"],
                          "current_resource_margin_conditioned_baseline_S":scores["S"]},
      "both_islands_nonempty_native_resource_sensitivity":{
                          "no_current_resource_margin_baseline_E":scores_native["E"],
                          "current_resource_margin_conditioned_baseline_S":scores_native["S"]},
      "limits":[
        "Both whole-alien-plant predictor and butterfly potential-resource response are constructed from overlapping contemporary WCVP species distributions and HOSTS interaction records; raw correlation and nested predictive gain are not independent trophic observations.",
        "The source does not contain true butterfly occupancy, site-level larval use, plant arrival dates, external independent butterfly community assemblies or direct ship/current trajectories.",
        "The model conditions on known source botanical species richness, and the stricter S additionally conditions on held-region contemporary butterfly potential-resource marginal richness; S is retrospective, not a new-region predictive forecast.",
        "All island dyads share endpoints; only five pre-defined code-hash folds, overlapping training and descriptive resampling uncertainty.",
        "Accepted alien whole-flora Jaccard includes documented HOSTS host plant species; a positive increment can reflect construction of the response rather than independent host-specific use.",
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    out=a.out_dir/"source_104_flora_to_butterfly_resource_strict_transfer.json"
    out.write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("WCVP104_ALIEN_FLORA_TO_POTENTIAL_BUTTERFLY_RESOURCE_BRIDGE",json.dumps({
      "status":summary["status"],"source":source,
      "raw":summary["raw_correlations_nonindependent"],
      "scores":{k:{m:{"gain_pct":v["macro_equal_fold_relative_gain_pct"],
                      "strict_test_pairs":v["n_test_pairs"],
                      "interval":v["descriptive_fold_bootstrap_95"],
                      "positive_folds":v["n_folds_adding_alien_flora_helped"]}
                  for m,v in q.items()}
                for k,q in [("all",scores),("native_supported",scores_native)]}
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

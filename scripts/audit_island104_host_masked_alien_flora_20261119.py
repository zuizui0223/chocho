#!/usr/bin/env python3
"""Host-masked alien-flora negative control for the frozen WCVP104 island resource model.

All flora and potential consumer resources are derived from WCVP/HOSTS, so the
analysis measures structural geography, not independent biology or colonization.
Original 104-island WCVP flora Jaccard is reproduced from individual accepted
plant IDs BEFORE splitting known HOSTS plants and unlisted WCVP plants.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import audit_island104_alien_flora_to_butterfly_resource_20261118 as bridge

SOURCE_RAW_SHA={
    "original_insect_host":"0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9",
    "native_allflora":"b66350a5d965e707cbfc5cb1de1d6fb72a62b6a01e5632eed2267b80cb02d091",
    "current_allflora":"bdf3ad670c435f570a715b23cb31ca7ab6add6e3241bee2105dc101fc675ebb7"
}
EXPECTED_SOURCE={
    "native_rows":1202784,"current_rows":1338575,
    "fixed_239_interaction_edges":2596,"known_239_accepted_host_plants":1706,
    "coastal_islands":104,"WCVP_shared_pairs":5355,
    "resource_defined_pairs":5353,"double_unseen_test_pairs":1088,
}
BOOT_SEED=20261119
N_BOOT=4999
MODELS=("B","W","H","O","H_O","W_H","W_O")
LEVELS=("E","S")

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        for row in csv.DictReader(f):
            yield row

def get_host_ids(path,panel):
    edges=set()
    for r in rows(path):
        insect=r["insect_species"].strip()
        if insect not in panel:continue
        p=r["accepted_plant_name_id"].strip()
        fam=r["family"].strip()
        if p and fam:edges.add((insect,p))
    if len(edges)!=EXPECTED_SOURCE["fixed_239_interaction_edges"]:
        raise RuntimeError(f"Original 239 HOSTS accepted interaction edges changed: {len(edges)}")
    if len({b for b,p in edges})!=239:
        raise RuntimeError("239 butterfly insect species not all represented")
    hosts={p for b,p in edges}
    if len(hosts)!=EXPECTED_SOURCE["known_239_accepted_host_plants"]:
        raise RuntimeError(f"Original 239 HOSTS accepted plant count changed: {len(hosts)}")
    return hosts

def read_region_plant_sets(native_file,current_file,codes):
    n={k:set() for k in codes}
    c={k:set() for k in codes}
    stats={}
    for which,p,bag in [("native",native_file,n),("current",current_file,c)]:
        nrows=0
        for r in rows(p):
            nrows+=1
            code=r["area_code_l3"].strip()
            if code in bag:
                pid=r["accepted_plant_name_id"].strip()
                if not pid:raise RuntimeError("Missing accepted source taxon ID")
                if pid in bag[code]:
                    raise RuntimeError(f"Duplicate accepted species×region: {which}/{code}/{pid}")
                bag[code].add(pid)
        stats[which+"_rows"]=nrows
    if stats["native_rows"]!=EXPECTED_SOURCE["native_rows"] or stats["current_rows"]!=EXPECTED_SOURCE["current_rows"]:
        raise RuntimeError(f"Original all-flora accepted plant sidecar rows differ: {stats}")
    for code in codes:
        if not n[code].issubset(c[code]):
            raise RuntimeError(f"Original native accepted WCVP plant IDs not nested in current for {code}")
    alien={code:c[code]-n[code] for code in codes}
    return n,c,alien,stats

def overlap(a,b):
    union=len(a|b)
    shared=len(a&b)
    return shared,union,(shared/union if union>0 else None)

def relative_gain(base,model):
    if base<=0:raise RuntimeError("Nonpositive held-fold MSE")
    return float(100*(base-model)/base)

def eval_models(entries,areas,H,O):
    folded=bridge.fold_cohort(entries,native_only=False)
    if len(folded)!=5 or sum(len(te) for _,tr,te in folded)!=1088:
        raise RuntimeError("Original held-source fold population changed")
    y=np.array([e["response"] for e in entries],dtype=float)
    h=np.array([e["host_alien_jaccard"] for e in entries],float)
    o=np.array([e["other_alien_jaccard"] for e in entries],float)
    w=np.array([e["alien_flora_Jaccard"] for e in entries],float)
    hmiss=np.array([e["host_undefined"] for e in entries],float)
    omiss=np.array([e["other_undefined"] for e in entries],float)
    result={}
    for level in LEVELS:
        current=(level=="S")
        # Reproduce exact previous source bridge E/S before introducing flags.
        b_original=bridge.features(entries,areas,include_current_margins=current,include_alien_similarity=False)
        w_original=bridge.features(entries,areas,include_current_margins=current,include_alien_similarity=True)
        baseline=bridge_model_score(b_original,w_original,y,folded)
        target={"E":26.36323764850138,"S":33.76782049915366}[level]
        if abs(baseline["gain_pct"]-target)>1e-8:
            raise RuntimeError(f"Original 104-island WCVP flora bridge not reproduced {level} {baseline['gain_pct']}")
        common=np.column_stack([b_original,hmiss,omiss])
        X={
            "B":common,
            "W":np.column_stack([common,w]),
            "H":np.column_stack([common,h]),
            "O":np.column_stack([common,o]),
            "H_O":np.column_stack([common,h,o]),
            "W_H":np.column_stack([common,w,h]),
            "W_O":np.column_stack([common,w,o]),
        }
        fold_out=[]
        for fold,tr,te in folded:
            predictions={}
            for m,x in X.items():
                pred=bridge.src.ridge(x[tr],y[tr],x[te])
                predictions[m]=float(np.mean((y[te]-pred)**2))
            fold_out.append({"fold":fold,"n_train_pairs":len(tr),"n_test_pairs":len(te),**{"MSE_"+m:v for m,v in predictions.items()}})
        equal={m:float(np.mean([z["MSE_"+m] for z in fold_out])) for m in MODELS}
        if any(v<=0 for v in equal.values()):raise RuntimeError("Nonpositive model MSE")
        gains={m:relative_gain(equal["B"],equal[m]) for m in MODELS[1:]}
        cvals={
            "O_to_O_plus_H":relative_gain(equal["O"],equal["H_O"]),
            "H_to_H_plus_O":relative_gain(equal["H"],equal["H_O"]),
            "W_to_W_plus_H":relative_gain(equal["W"],equal["W_H"]),
            "W_to_W_plus_O":relative_gain(equal["W"],equal["W_O"])
        }
        # Source fold-only descriptive Monte Carlo intervals; same fold draws for all models.
        rng=np.random.default_rng(BOOT_SEED+(0 if level=="E" else 1))
        mat=np.array([[z["MSE_"+m] for m in MODELS] for z in fold_out],float)
        draw={k:[] for k in list(gains)+list(cvals)}
        for _ in range(N_BOOT):
            idx=rng.integers(0,5,5)
            vals=dict(zip(MODELS,mat[idx].mean(axis=0)))
            for m in gains:draw[m].append(relative_gain(vals["B"],vals[m]))
            draw["O_to_O_plus_H"].append(relative_gain(vals["O"],vals["H_O"]))
            draw["H_to_H_plus_O"].append(relative_gain(vals["H"],vals["H_O"]))
            draw["W_to_W_plus_H"].append(relative_gain(vals["W"],vals["W_H"]))
            draw["W_to_W_plus_O"].append(relative_gain(vals["W"],vals["W_O"]))
        weight=np.array([f["n_test_pairs"] for f in fold_out])
        weighted={m:float(np.average([f["MSE_"+m] for f in fold_out],weights=weight)) for m in MODELS}
        result[level]={
            "source_original_whole_alien_gain_without_flags_pct":baseline["gain_pct"],
            "source_original_whole_alien_descriptive_interval":baseline["interval"],
            "n_test_pairs":int(weight.sum()),
            "base_with_both_empty_flags_MSE":equal["B"],
            "MSE_by_model_equal_fold":equal,
            "gain_pct_vs_common_flagged_base":gains,
            "incremental_gain_pct":cvals,
            "descriptive_fivefold_bootstrap95":{
                k:list(map(float,np.quantile(v,[.025,.975]))) for k,v in draw.items()
            },
            "n_improving_folds_vs_common_flagged_base":{
                m:sum(z["MSE_"+m]<z["MSE_B"] for z in fold_out) for m in MODELS[1:]
            },
            "pair_weighted_gain_vs_common_flagged_base_pct":{
                m:relative_gain(weighted["B"],weighted[m]) for m in MODELS[1:]
            },
            "source_n_folds":5,"folds":fold_out,
        }
    return result

def bridge_model_score(X0,X1,y,folds):
    bb=[];pp=[]
    for _,tr,te in folds:
        pred0=bridge.src.ridge(X0[tr],y[tr],X0[te])
        pred1=bridge.src.ridge(X1[tr],y[tr],X1[te])
        bb.append(float(np.mean((y[te]-pred0)**2)))
        pp.append(float(np.mean((y[te]-pred1)**2)))
    b=float(np.mean(bb));p=float(np.mean(pp))
    return {"gain_pct":relative_gain(b,p),"interval":None}

def main():
    ap=argparse.ArgumentParser()
    for k in ("archived-root","plants","geometry","island-area","allflora-native","allflora-current","insect-hosts","out-dir"):
        ap.add_argument("--"+k,required=True,type=Path)
    a=ap.parse_args()
    source=a.archived_root/"provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/beta_baseline_v3.json"
    if bridge.sha(source)!=bridge.REPORT_SHA:
        raise RuntimeError("Frozen source model identity changed")
    q=json.loads(source.read_text(encoding="utf-8"))
    order=q["resource_matrix"]["species_order"]
    if len(order)!=239 or len(set(order))!=239:
        raise RuntimeError("Frozen original butterfly panel species order changed")
    h_sha=digest(a.insect_hosts);n_sha=digest(a.allflora_native);c_sha=digest(a.allflora_current)
    for kind,actual,expected in [
        ("original_insect_host",h_sha,SOURCE_RAW_SHA["original_insect_host"]),
        ("native_allflora",n_sha,SOURCE_RAW_SHA["native_allflora"]),
        ("current_allflora",c_sha,SOURCE_RAW_SHA["current_allflora"])
    ]:
        if actual!=expected:
            raise RuntimeError(f"Original 2026-10-31 accepted source sidecar sha drift: {kind} {actual} != {expected}")
    hosts=get_host_ids(a.insect_hosts,set(order))
    entries,area,audit=bridge.setup(a.archived_root,a.plants,a.geometry,a.island_area)
    if audit["source_pairs"]!=5355 or len(entries)!=5353 or len(set([e["island1"] for e in entries]+[e["island2"] for e in entries]))!=104:
        raise RuntimeError("Prior 104 island botanical/resource cohort not recreated")
    # The plant CSV includes more botanical region pairs than the original 104-island
    # coastal geometry cohort. Select the 5,355 matched pairs using the FROZEN
    # source geometry join, not the full source CSV nor response-defined entries.
    raw_plant_rows=list(rows(a.plants))
    source_pair={tuple(sorted((r["source_region"],r["destination_region"]))):r for r in raw_plant_rows}
    if len(source_pair)!=len(raw_plant_rows):
        raise RuntimeError("Original full plant CSV has duplicate unordered pair keys")
    original_matched_pairs,unmatched=bridge.src.input_pair(
        bridge.SimpleNamespace(plants=a.plants,geometry=a.geometry))
    original_matched_keys=[tuple(sorted((e["island1"],e["island2"]))) for e in original_matched_pairs]
    if len(original_matched_keys)!=5355 or len(set(original_matched_keys))!=5355 or len(unmatched)!=1:
        raise RuntimeError("Original 5,355 unique matched coastal flora dyads not recreated")
    if any(k not in source_pair for k in original_matched_keys):
        raise RuntimeError("Original matched coastal flora dyad absent from plant source CSV")
    source_pair_rows=[source_pair[k] for k in original_matched_keys]
    src_codes={e["island1"] for e in entries}|{e["island2"] for e in entries}
    if {r["source_region"] for r in source_pair_rows}|{r["destination_region"] for r in source_pair_rows}!=src_codes:
        raise RuntimeError("Original matched coastal flora island identities differ from resource cohort")
    # Validate ALL matched plant dyads, including two whose resource response is undefined.
    native,current,alien,flora_rows=read_region_plant_sets(a.allflora_native,a.allflora_current,src_codes)
    count_compared=0
    failures=[]
    for pair in source_pair_rows:
        ia,ib=pair["source_region"],pair["destination_region"]
        x,y=alien[ia],alien[ib]
        n_shared,n_union,n_j=overlap(x,y)
        expected=(int(pair["source_alien_WCVP_accepted"]),int(pair["dest_alien_WCVP_accepted"]))
        if (len(x),len(y))!=expected or n_shared!=int(pair["alien_shared_accepted_taxa"]) or n_union!=int(pair["alien_union_accepted_taxa"]) or n_j is None or abs(n_j-float(pair["alien_Jaccard"]))>1e-12:
            failures.append({"source":ia,"dest":ib,"computed":[len(x),len(y),n_shared,n_union,n_j],"original":[*expected,pair["alien_shared_accepted_taxa"],pair["alien_union_accepted_taxa"],pair["alien_Jaccard"]]})
            if len(failures)>5:break
        count_compared+=1
    if failures:
        raise RuntimeError("Original accepted WCVP alien Jaccard cannot be reconstructed from pinned original raw sidecars: "+json.dumps(failures[:5]))
    if count_compared!=5355:raise RuntimeError("Host-mask original full WCVP source cohort not matched")
    # Model only the 5,353 original dyads with defined butterfly-resource Jaccard.
    for e in entries:
        ia,ib=e["island1"],e["island2"]
        x,y=alien[ia],alien[ib]
        H1,H2=x&hosts,y&hosts
        O1,O2=x-hosts,y-hosts
        sh,un,jh=overlap(H1,H2)
        so,uo,jo=overlap(O1,O2)
        e.update(host_alien_jaccard=(jh if jh is not None else 0.0),
                 other_alien_jaccard=(jo if jo is not None else 0.0),
                 host_undefined=int(jh is None),other_undefined=int(jo is None),
                 host_union_size=un,other_union_size=uo)
    result=eval_models(entries,area,None,None)
    info={
        "schema":"chocho_104_WCVP_island_HOSTS239_host_masked_other_alien_flora_20261119_v1",
        "status":"SOURCE_MATCHED_EXPLORATORY_HOST_MASK_VS_UNLISTED_ALIEN_PLANT_COMPOSITION",
        "source_original_WCVP_accepted_sidecars_sha256":SOURCE_RAW_SHA,
        "source_original_239_butterflies":239,
        "source_original_239_resolved_HOSTS_edges":2596,
        "source_original_239_distinct_host_accepted_plant_IDs":len(hosts),
        "source_original_frozen_matrix_sha256":bridge.SOURCE_SHA,
        "source_original_allflora_region_rows":flora_rows,
        "source_104_island_pairs_WCVP_origin_recomputed_exactly":count_compared,
        "source_resource_jaccard_defined_pairs":len(entries),
        "source_strict_both_endpoints_unseen_test_pairs":1088,
        "source_104_island_alien_taxa_counts":{
            "distinct_all_alien_accepted_IDs":len(set().union(*alien.values())),
            "distinct_alien_239HOSTS_IDs":len(set().union(*(a&hosts for a in alien.values()))),
            "distinct_alien_other_IDs":len(set().union(*(a-hosts for a in alien.values()))),
            "total_alien_species_x_region":sum(len(v) for v in alien.values()),
            "host_known_alien_species_x_region":sum(len(v&hosts) for v in alien.values()),
            "other_alien_species_x_region":sum(len(v-hosts) for v in alien.values()),
        },
        "source_plant_mask_pair_undefined_union_counts":{
            "host":sum(e["host_undefined"] for e in entries),
            "other":sum(e["other_undefined"] for e in entries),
        },
        "model_likelihood":{
            "E_non_target_richness":result["E"],
            "S_with_target_derived_current_resource_richness":result["S"]
        },
        "limits":[
            "Other alien species are only NOT LISTED as hosts for this selected 239-butterfly HOSTS catalogue; this is not proved biological nonhost status.",
            "Response of potential 239-butterfly host resource geography and source-vetted host predictor are structurally coupled by original WCVP+HOSTS data.",
            "Original source WCVP accepted distribution status lacks date of introduction and sites of actual larval feeding.",
            "Only five dependent source-hash geographic folds; confidence intervals are descriptive and not independent source replication.",
            "Contemporary butterfly resource richness is source-target-derived, so S is retrospective.",
            "Additional H/O predictors and missingness indicators may be correlated; a positive host incremental gain does not identify independent host use.",
        ],
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    target=a.out_dir/"island104_host_masked_alien_flora_bridge_source_scores.json"
    target.write_text(json.dumps(info,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("WCVP104_SOURCE_HOST_MASKED_ALIEN_FLORA",json.dumps({
        "source_validated":count_compared,
        "plant_counts":info["source_104_island_alien_taxa_counts"],
        "missing":info["source_plant_mask_pair_undefined_union_counts"],
        "results":{k:{
            "baseline_E_or_S_original_whole_gain_pct":v["source_original_whole_alien_gain_without_flags_pct"],
            "gain_by_mask":v["gain_pct_vs_common_flagged_base"],
            "conditional_H_and_O":v["incremental_gain_pct"],
            "sample_n":v["n_test_pairs"]
        } for k,v in result.items()}
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

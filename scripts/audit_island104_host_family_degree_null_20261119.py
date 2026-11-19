#!/usr/bin/env python3
"""Family × island-occupancy matched pseudo-HOSTS null; exploratory SAME-source diagnostic."""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np

import audit_island104_host_masked_alien_flora_20261119 as mask
import audit_island104_alien_flora_to_butterfly_resource_20261118 as bridge

SEED = 20261119
N_RANDOM = 199
FOLD_EXPECTED = [405, 231, 190, 91, 171]
ORIGINAL_GAINS = {
    "E": {"W":26.363237648501535, "H":33.8214798093749,
          "O":19.53859900994786, "H_O":33.707207863415725},
    "S": {"W":33.76782049915363, "H":41.80431031928996,
          "O":26.111822764290107, "H_O":41.58352027778727},
}
MODELS = ("B", "H", "O", "H_O")
STATS = ("H_given_O", "O_given_H", "B_to_H", "B_to_O", "H_minus_O")

def degree_bin(deg):
    if deg < 1 or deg > 104:
        raise RuntimeError("Alien taxon outside 1-104 island range")
    return min(6, deg.bit_length() - 1)

def source_family_map(path, all_plant_ids):
    result = {}
    for r in mask.rows(path):
        p = r["accepted_plant_name_id"].strip()
        if not p or p not in all_plant_ids:
            continue
        family = r["family"].strip()
        if not family:
            family = "UNRESOLVED_WCVP_FAMILY"
        if p in result and result[p] != family:
            raise RuntimeError(f"Accepted WCVP taxon family conflict for {p}")
        result[p] = family
    missing = all_plant_ids - result.keys()
    if missing:
        raise RuntimeError(f"WCVP accepted-family map missing {len(missing)} alien taxa")
    return result

def build_vectors(bitmask, island_bits, codes, entry_pairs, fullmask):
    by_code_h = [int(x & bitmask) for x in island_bits]
    by_code_o = [int(x & (fullmask ^ bitmask)) for x in island_bits]
    n_h = np.array([x.bit_count() for x in by_code_h], dtype=int)
    n_o = np.array([x.bit_count() for x in by_code_o], dtype=int)
    h = np.zeros(len(entry_pairs), dtype=float)
    o = np.zeros(len(entry_pairs), dtype=float)
    hm = np.zeros(len(entry_pairs), dtype=float)
    om = np.zeros(len(entry_pairs), dtype=float)
    margins = np.empty((len(entry_pairs),4),dtype=float)
    for k,(i,j) in enumerate(entry_pairs):
        ha,hb = by_code_h[i],by_code_h[j]
        oa,ob = by_code_o[i],by_code_o[j]
        hu=(ha|hb).bit_count()
        ou=(oa|ob).bit_count()
        if hu:
            h[k] = (ha&hb).bit_count()/hu
        else:
            hm[k] = 1.
        if ou:
            o[k] = (oa&ob).bit_count()/ou
        else:
            om[k] = 1.
        margins[k] = [math.log1p(min(n_h[i],n_h[j])),
                      math.log1p(max(n_h[i],n_h[j])),
                      math.log1p(min(n_o[i],n_o[j])),
                      math.log1p(max(n_o[i],n_o[j]))]
    return h,o,hm,om,margins,n_h,n_o

def measure(vectors, E, S, y, folds):
    h,o,hm,om,margins,nh,no=vectors
    output={}
    for level,base in (("E",E),("S",S)):
        common=np.column_stack([base,hm,om,margins])
        X={
            "B":common,
            "H":np.column_stack([common,h]),
            "O":np.column_stack([common,o]),
            "H_O":np.column_stack([common,h,o]),
        }
        fold_results=[]
        for fold,tr,te in folds:
            rec={"fold":int(fold),"n_test_pairs":int(len(te))}
            for model in MODELS:
                pred=bridge.src.ridge(X[model][tr],y[tr],X[model][te])
                rec["MSE_"+model]=float(np.mean((y[te]-pred)**2))
            fold_results.append(rec)
        equal={model:float(np.mean([f["MSE_"+model] for f in fold_results]))
               for model in MODELS}
        weights=np.array([f["n_test_pairs"] for f in fold_results],dtype=float)
        weighted={model:float(np.average([f["MSE_"+model] for f in fold_results],
                                         weights=weights)) for model in MODELS}
        def stats(mse):
            B,H,O,HO = [mse[m] for m in MODELS]
            return {"H_given_O":mask.relative_gain(O,HO),
                    "O_given_H":mask.relative_gain(H,HO),
                    "B_to_H":mask.relative_gain(B,H),
                    "B_to_O":mask.relative_gain(B,O),
                    "H_minus_O":mask.relative_gain(B,H)-mask.relative_gain(B,O)}
        output[level]={
            "equal_fold":stats(equal),
            "pair_weighted":stats(weighted),
            "equal_fold_MSE":equal,
            "folds":fold_results,
            "null_union_undefined_H":int(hm.sum()),
            "null_union_undefined_O":int(om.sum()),
            "host_region_richness_sum":int(nh.sum()),
            "other_region_richness_sum":int(no.sum()),
        }
    return output

def main():
    ap=argparse.ArgumentParser()
    for name in ("archived-root","plants","geometry","island-area",
                 "allflora-native","allflora-current","insect-hosts",
                 "families","out-dir"):
        ap.add_argument("--"+name,required=True,type=Path)
    a=ap.parse_args()
    for kind,path in (("original_insect_host",a.insect_hosts),
                      ("native_allflora",a.allflora_native),
                      ("current_allflora",a.allflora_current)):
        if mask.digest(path)!=mask.SOURCE_RAW_SHA[kind]:
            raise RuntimeError(f"Frozen WCVP/HOSTS source drift: {kind}")
    receipt=a.archived_root/"provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/beta_baseline_v3.json"
    if mask.digest(receipt)!=bridge.REPORT_SHA:
        raise RuntimeError("Frozen HOSTS resource matrix receipt changed")
    species_order=json.loads(receipt.read_text())["resource_matrix"]["species_order"]
    if len(species_order)!=239 or len(set(species_order))!=239:
        raise RuntimeError("Frozen butterfly species panel changed")
    host_ids=mask.get_host_ids(a.insect_hosts,set(species_order))
    entries,areas,source_audit=bridge.setup(a.archived_root,a.plants,a.geometry,a.island_area)
    if (source_audit["source_pairs"]!=5355 or len(entries)!=5353):
        raise RuntimeError("Coastal WCVP/resource cohort changed")
    codes=sorted({e["island1"] for e in entries}|{e["island2"] for e in entries})
    if len(codes)!=104:
        raise RuntimeError("Original 104 islands not reproduced")
    native,current,alien,_=mask.read_region_plant_sets(a.allflora_native,a.allflora_current,set(codes))
    plants=sorted(set().union(*alien.values()))
    plant_set=set(plants)
    true_host_set=plant_set & host_ids
    if len(plants)!=7887 or len(true_host_set)!=668 or len(plant_set-true_host_set)!=7219:
        raise RuntimeError("Original WCVP alien species host-mask split changed")
    if sum(len(v&host_ids) for v in alien.values())!=7726:
        raise RuntimeError("Original WCVP alien host island incidence changed")

    # Source audit: original coastal geography JOIN, all 5,355 source dyads.
    source_plant_rows={tuple(sorted((r["source_region"],r["destination_region"]))):r for r in mask.rows(a.plants)}
    full_pairs,unmatched=bridge.src.input_pair(bridge.SimpleNamespace(plants=a.plants,geometry=a.geometry))
    if len(full_pairs)!=5355 or len(unmatched)!=1:
        raise RuntimeError("Original 5,355 source coastal plant dyads changed")
    for e in full_pairs:
        i,j=e["island1"],e["island2"]
        r=source_plant_rows[tuple(sorted((i,j)))]
        sh,un,jac=mask.overlap(alien[i],alien[j])
        src1=int(r["source_alien_WCVP_accepted"])
        src2=int(r["dest_alien_WCVP_accepted"])
        if (len(alien[r["source_region"]])!=src1 or len(alien[r["destination_region"]])!=src2
            or sh!=int(r["alien_shared_accepted_taxa"])
            or un!=int(r["alien_union_accepted_taxa"])
            or jac is None or abs(jac-float(r["alien_Jaccard"]))>1e-12):
            raise RuntimeError(f"Original WCVP Jaccard failed: {i}/{j}")

    family=source_family_map(a.families,plant_set)
    # Exact accepted WCVP family agreement with the pinned 239 HOSTS links.
    disagree=0
    for r in mask.rows(a.insect_hosts):
        if r["insect_species"].strip() in species_order:
            p=r["accepted_plant_name_id"].strip()
            if p in true_host_set and family[p]!=r["family"].strip():
                disagree+=1
    if disagree:
        raise RuntimeError(f"WCVP family disagrees with frozen HOSTS on {disagree} rows")

    code_idx={c:i for i,c in enumerate(codes)}
    plant_idx={p:i for i,p in enumerate(plants)}
    island_bits=[]
    for code in codes:
        z=0
        for p in alien[code]:
            z|=1<<plant_idx[p]
        island_bits.append(z)
    fullmask=(1<<len(plants))-1
    degree={p:sum(p in alien[code] for code in codes) for p in plants}
    strata=defaultdict(list)
    for p in plants:
        strata[(family[p],degree_bin(degree[p]))].append(plant_idx[p])
    original_hmask=0
    for p in true_host_set:
        original_hmask|=1<<plant_idx[p]
    keys=sorted(strata)
    groups=[(k,sorted(strata[k]),sum(plants[i] in true_host_set for i in strata[k])) for k in keys]
    fixed_host=sum(nh for _,members,nh in groups if nh==len(members))
    mutable_host=668-fixed_host
    if mutable_host<=0:
        raise RuntimeError("No host label variation under predeclared family × occupancy null")

    pair_indices=[(code_idx[e["island1"]],code_idx[e["island2"]]) for e in entries]
    y=np.array([e["response"] for e in entries],dtype=float)
    folded=bridge.fold_cohort(entries,native_only=False)
    sizes=[int(len(te)) for fold,tr,te in folded]
    if sizes!=FOLD_EXPECTED:
        raise RuntimeError(f"Frozen strict held-out folds changed: {sizes}")
    E=bridge.features(entries,areas,include_current_margins=False,include_alien_similarity=False)
    S=bridge.features(entries,areas,include_current_margins=True,include_alien_similarity=False)
    original_vectors=build_vectors(original_hmask,island_bits,codes,pair_indices,fullmask)
    for k,e in enumerate(entries):
        e["host_alien_jaccard"]=float(original_vectors[0][k])
        e["other_alien_jaccard"]=float(original_vectors[1][k])
        e["host_undefined"]=int(original_vectors[2][k])
        e["other_undefined"]=int(original_vectors[3][k])
    old=mask.eval_models(entries,areas,None,None)
    for level in ("E","S"):
        for model,expected in ORIGINAL_GAINS[level].items():
            got=old[level]["gain_pct_vs_common_flagged_base"][model]
            if abs(got-expected)>1e-8:
                raise RuntimeError(f"Previous host-mask published score did not reproduce: {level}/{model}")

    observed=measure(original_vectors,E,S,y,folded)
    rng=np.random.default_rng(SEED)
    null={level:{stat:[] for stat in STATS} for level in ("E","S")}
    null_weighted={level:{stat:[] for stat in STATS} for level in ("E","S")}
    missing_union_counts=[]
    for draw in range(N_RANDOM):
        bits=0
        for _,group,k in groups:
            if not k:
                continue
            if k==len(group):
                chosen=group
            else:
                chosen=rng.choice(np.array(group,dtype=np.int32),size=k,replace=False)
            for idx in chosen:
                bits|=1<<int(idx)
        if bits.bit_count()!=668:
            raise RuntimeError("Permutation did not preserve original 668 host labels")
        vectors=build_vectors(bits,island_bits,codes,pair_indices,fullmask)
        v=measure(vectors,E,S,y,folded)
        missing_union_counts.append([int(v["E"]["null_union_undefined_H"]),int(v["E"]["null_union_undefined_O"])])
        for level in ("E","S"):
            for stat in STATS:
                null[level][stat].append(v[level]["equal_fold"][stat])
                null_weighted[level][stat].append(v[level]["pair_weighted"][stat])

    summary={}
    for level in ("E","S"):
        summary[level]={
            "observed":observed[level],
            "null_family_occupancy_199":{
                stat:{"mean":float(np.mean(vals)),
                      "median":float(np.median(vals)),
                      "q025_q975":[float(x) for x in np.quantile(vals,[.025,.975])],
                      "min":float(min(vals)),"max":float(max(vals)),
                      "n_at_least_observed":int(sum(x>=observed[level]["equal_fold"][stat] for x in vals)),
                      "descriptive_upper_tail_add_one":float((1+sum(x>=observed[level]["equal_fold"][stat] for x in vals))/(N_RANDOM+1))}
                for stat,vals in null[level].items()},
            "null_pair_weighted":{
                stat:{"median":float(np.median(vals)),
                      "q025_q975":[float(x) for x in np.quantile(vals,[.025,.975])]}
                for stat,vals in null_weighted[level].items()}
        }
    out={
        "schema":"chocho_WCVP104_HOSTS239_family_degree_null_20261119_v1",
        "status":"SOURCE_LOCKED_EXPLORATORY_STRUCTURAL_HOST_IDENTITY_NULL",
        "original_plant_source_pairs_exact":5355,
        "original_resource_defined_pairs":5353,
        "test_both_endpoints_unseen_pairs":1088,
        "test_fold_sizes":sizes,
        "source_butterflies":239,
        "source_accepted_plant_IDs":7887,
        "source_host_accepted_IDs":668,
        "source_nonlisted_accepted_IDs":7219,
        "null_seed":SEED,"null_draws":N_RANDOM,
        "family_occupancy_strata":{
            "strata_total":len(groups),
            "strata_with_both_host_and_other":sum(0<nh<len(members) for _,members,nh in groups),
            "strata_fully_host":sum(nh==len(members) for _,members,nh in groups),
            "host_labels_immovable_due_full_host_strata":fixed_host,
            "host_labels_in_mixed_strata":mutable_host,
            "missing_WCVP_family_taxa":sum(family[p]=="UNRESOLVED_WCVP_FAMILY" for p in plants),
            "occupancy_bin_rule":"floor(log2(number_of_104_islands_present)), capped at 6",
        },
        "null_undefined_union_summary":{"H_min":min(x[0] for x in missing_union_counts),
           "H_max":max(x[0] for x in missing_union_counts),
           "O_min":min(x[1] for x in missing_union_counts),
           "O_max":max(x[1] for x in missing_union_counts)},
        "scores":summary,
        "interpretation_guardrails":[
            "A result above this null is still structurally linked because HOSTS host names define the potential butterfly-resource response.",
            "Pseudo host label assignment preserves family and broad global island occupancy-bin counts, not exact island-specific host margins; those are modeled explicitly.",
            "Null draws share one island/plant/butterfly population, so upper-tail ranks are descriptive, not independent biological significance.",
            "O and pseudo-O mean plants not assigned HOSTS host labels, not independently experimentally validated nonhosts.",
            "S conditions on response-derived contemporary butterfly-resource marginal richness and is retrospective.",
            "No observed butterfly survival, host preference, abundance, or new establishment is measured."
        ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    path=a.out_dir/"family_occupancy_matched_host_label_null_20261119.json"
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("FAMILY_OCCUPANCY_MATCHED_HOST_NULL",json.dumps({
        "source_pairs":5355,"held_pairs":1088,"families":len({family[p] for p in plants}),
        "mixed_strata":out["family_occupancy_strata"]["strata_with_both_host_and_other"],
        "immovable_hosts":fixed_host,
        "E_observed":observed["E"]["equal_fold"],
        "E_null_H_given_O":summary["E"]["null_family_occupancy_199"]["H_given_O"],
        "S_observed":observed["S"]["equal_fold"],
        "S_null_H_given_O":summary["S"]["null_family_occupancy_199"]["H_given_O"],
    }),flush=True)

if __name__=="__main__":
    main()

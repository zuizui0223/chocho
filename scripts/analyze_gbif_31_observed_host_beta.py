#!/usr/bin/env python3
"""Independent OBSERVED butterfly-assemblage pilot versus GloNAF larval hosts.

Only 30 species in GBIF's archived 31-species postgate sample were matched to
the frozen 239-butterfly resource panel. These records are an independent
OUTCOME, unlike HOSTS-derived potential-resource Jaccard.

No GBIF non-record is a confirmed species absence. Modern GloNAF/WCVP flora
is not colonial-time plant chronology. Region-held-out Ridge comparison is
spatial transport of contemporaneous record similarity, not a future forecast.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
from shapely.geometry import Point

from audit_gbif_31_butterfly_wgsrpd3_source import (
    SOURCE as GBIF_ZIP, GEOM, PANEL, REGIONS, EXPECTED_RETURNED,
    find_occurrence_member, geometry, normalize_species, read_csv as load_csv,
    sha256 as file_sha,
)
from analyze_documented_host_specificity import (
    HOSTS, GLO_NAF, TDWG, binomial, glonaf_index, pair_jaccard,
)
from analyze_matched_glonaf_consumer_history import produce_matched_regions
from analyze_colonial_legacy_resource_beta import (
    load as load_beta_sources, orthogonal_controls, residual, correlation,
)

ROOT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
LOCK=ROOT/"gbif_31_observed_beta_model_lock_v1.json"
GATE=ROOT/"gbif_31_observed_region_source_gate.json"
REGION_COUNTS=ROOT/"gbif_31_observed_region_counts.csv"
OUT=ROOT/"gbif_31_observed_host_beta_pilot_v1.json"
PAIR_CSV=ROOT/"gbif_31_observed_host_beta_pairs_v1.csv"


def canonical_exact_hosts(selected_species,hosts=HOSTS):
    host_by_butterfly=defaultdict(set)
    with hosts.open(encoding="utf8-sig",newline="") as handle:
        for r in csv.DictReader(handle):
            b=(r.get("Insect Genus","").strip()+" "+
               r.get("Insect Species","").strip()).casefold()
            if b not in selected_species:
                continue
            name=binomial(r.get("Hostplant Genus"),r.get("Hostplant Species"))
            if name:
                host_by_butterfly[b].add(name)
    if len(host_by_butterfly)<20:
        raise RuntimeError("Too few original HOSTS butterflies in GBIF source")
    return host_by_butterfly


def observed_assemblages(source,source_gate,region_gate,region_codes):
    if file_sha(source)!=source_gate["source"]["archive_sha256"]:
        raise RuntimeError("Independent GBIF DOI download bytes changed")
    tree,shapes,codes=geometry()
    allowed=set(region_codes)
    candidates={r["species"].strip().casefold() for r in load_csv(PANEL)}
    species_by_region=defaultdict(set)
    early_by_region=defaultdict(set)
    late_by_region=defaultdict(set)
    counts=Counter()
    gbif_species=set()
    for k in region_codes:
        species_by_region[k]=set()
    with zipfile.ZipFile(source) as z:
        (member,delim,_),_members=find_occurrence_member(z)
        with z.open(member) as f:
            reader=csv.DictReader(io.TextIOWrapper(f,encoding="utf-8-sig",newline=""),delimiter=delim)
            for original in reader:
                counts["rows"]+=1
                row={k.casefold():v for k,v in original.items() if k}
                name=normalize_species(row)
                if name not in candidates:
                    continue
                gbif_species.add(name)
                try:
                    lat=float(row.get("decimallatitude") or "")
                    lon=float(row.get("decimallongitude") or "")
                except ValueError:
                    continue
                if not (-90<=lat<=90 and -180<=lon<=180 and
                        math.isfinite(lat) and math.isfinite(lon)):
                    continue
                pt=Point(lon,lat)
                hits=[codes[int(j)] for j in tree.query(pt) if shapes[int(j)].covers(pt)]
                if not hits:
                    continue
                code=sorted(hits)[0]
                if code not in allowed:
                    continue
                species_by_region[code].add(name)
                counts["matched_region_records"]+=1
                try:
                    year=int(row.get("year") or "")
                except ValueError:
                    continue
                if 2010<=year<=2017:
                    early_by_region[code].add(name)
                if 2018<=year<=2026:
                    late_by_region[code].add(name)
    if counts["rows"]!=EXPECTED_RETURNED:
        raise RuntimeError("Archived GBIF records differed from source gate")
    if len(gbif_species)!=source_gate["coverage"]["source_taxa_matching_frozen_239"]:
        raise RuntimeError("GBIF species retained count inconsistent with original source gate")
    for row in region_gate:
        key=row["wgsrpd3_code"]
        if key not in allowed:
            continue
        if len(species_by_region[key])!=int(row["indexed_butterfly_species"]):
            raise RuntimeError("Independent reconstruction disagrees with source-gate species richness at "+key)
    return species_by_region,early_by_region,late_by_region,gbif_species,dict(counts)


def jaccard_from_sets(assemblages,ids):
    # Species order fixed lexically; the finite 30-species window does not
    # represent the full 239 species and is never treated as true absence.
    names=sorted(set().union(*assemblages.values()))
    indexing={taxon:i for i,taxon in enumerate(names)}
    masks=[]
    for key in ids:
        p=0
        for name in assemblages.get(key,()):
            p |= 1<<indexing[name]
        masks.append(p)
    return masks,names


def ji(masks,ia,ib):
    arr,missing=pair_jaccard(masks,ia,ib,undefined_zero=True)
    return arr,int(missing)


def fit_ridge_predict(train,test,y,X,penalty=10.):
    """No outcome-derived species richness; standardize on training pairs."""
    xtr=X[train]
    means=xtr.mean(axis=0)
    std=xtr.std(axis=0)
    std=np.where(std>1e-10,std,1.)
    a=(xtr-means)/std
    b=(X[test]-means)/std
    a=np.column_stack([np.ones(len(a)),a])
    b=np.column_stack([np.ones(len(b)),b])
    gram=a.T@a
    diag=np.eye(gram.shape[0])*penalty
    diag[0,0]=0
    beta=np.linalg.solve(gram+diag,a.T@y[train])
    return np.clip(b@beta,0,1)


def evaluate_loio(y,features,ia,ib,node_count,lam):
    out={}
    for model,feature_arrays in features.items():
        design=np.column_stack(feature_arrays)
        by_node=[]
        for idx in range(node_count):
            held=(ia==idx)|(ib==idx)
            if held.sum()!=node_count-1:
                raise RuntimeError("Pair sample is not the complete nodal cohort")
            preds=fit_ridge_predict(~held,held,y,design,lam)
            error=preds-y[held]
            by_node.append({
                "node":idx,
                "n_heldout_pairs":int(held.sum()),
                "mse":float(np.mean(error**2)),
                "mae":float(np.mean(np.abs(error))),
            })
        out[model]={
            "island_or_region_equal_mse":float(np.mean([x["mse"] for x in by_node])),
            "region_equal_mae":float(np.mean([x["mae"] for x in by_node])),
            "per_region":by_node,
        }
    comparisons={}
    for target in ("M1_known_host","M2_other_aliens","M3_both"):
        a=out["M0_general_flora_and_effort"]["per_region"]
        b=out[target]["per_region"]
        delta=[x["mse"]-y["mse"] for x,y in zip(a,b)]
        comparisons[target]={
            "mse_reduction_percent":100*float(np.mean(delta))/out["M0_general_flora_and_effort"]["island_or_region_equal_mse"],
            "mean_region_equal_mse_difference":float(np.mean(delta)),
            "regions_improved":sum(x>0 for x in delta),
            "regions_worsened":sum(x<0 for x in delta),
            "median_improvement_absolute":float(np.median(delta)),
            "minimum_improvement":min(delta),
            "maximum_improvement":max(delta),
        }
    summaries={k:{kk:v for kk,v in x.items() if kk!="per_region"} for k,x in out.items()}
    return summaries,comparisons


def analyse():
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock["schema"]!="chocho_independent_GBIF31_observed_butterfly_beta_284source_lock_v1":
        raise RuntimeError("Predetermined observed-butterfly test design absent")
    source_gate=json.loads(GATE.read_text(encoding="utf8"))
    if source_gate["pilot_gate"]!="REGION_BETA_FEASIBLE_FOR_LIMITED_31_SPECIES":
        raise RuntimeError("Original independently observed species coverage gate rejected sample")
    original_regions=load_csv(REGION_COUNTS)
    good=[row for row in original_regions
          if int(row["indexed_butterfly_species"])>=5 and
             int(row["gbif_provider_days"])>=5]
    if len(good)!=111:
        raise RuntimeError("Prelocked 111 record-eligible WGSRPD3 regions changed")
    good_codes={r["wgsrpd3_code"] for r in good}
    region_qa={r["wgsrpd3_code"]:r for r in good}
    beta=load_beta_sources()
    if len(beta["ids"])!=284:
        raise RuntimeError("Frozen colonial/geographic beta cohort changed")
    selected_mask=np.array([code in good_codes for code in beta["ids"]])
    if sum(selected_mask)!=111:
        raise RuntimeError("Source GBIF/colonial WGSRPD3 mapping inconsistent")
    selected_codes=[k for k in beta["ids"] if k in good_codes]
    q=(selected_mask[beta["ia"]] & selected_mask[beta["ib"]])
    ia0,ib0=beta["ia"][q],beta["ib"][q]
    rebase={j:k for k,j in enumerate(np.flatnonzero(selected_mask))}
    ia=np.asarray([rebase[int(j)] for j in ia0],int)
    ib=np.asarray([rebase[int(j)] for j in ib0],int)
    if len(ia)!=6105 or len(ib)!=6105:
        raise RuntimeError("GBIF eligibility cannot produce expected 111 choose2 pairs")

    species_all,early,late,selected_taxa,counts=observed_assemblages(
        GBIF_ZIP,source_gate,original_regions,selected_codes)
    observations,stored_species=jaccard_from_sets(species_all,selected_codes)
    first,first_species=jaccard_from_sets(early,selected_codes)
    second,second_species=jaccard_from_sets(late,selected_codes)
    observed_j,obs_missing=ji(observations,ia,ib)
    first_j,early_missing=ji(first,ia,ib)
    second_j,late_missing=ji(second,ia,ib)
    if obs_missing:
        raise RuntimeError("Every study region must have at least 5 recorded butterflies")
    if len(stored_species)<20:
        raise RuntimeError("Matched observed species too sparse to evaluate")
    source_masks,_,source_audit=produce_matched_regions(GLO_NAF,TDWG,(.05,))
    aliases,glonaf_audit=glonaf_index(GLO_NAF)
    host_by_consumer=canonical_exact_hosts(selected_taxa,HOSTS)
    native_host_names=set().union(*host_by_consumer.values())
    matched_aliases=native_host_names.intersection(aliases)
    known_idx={j for name in matched_aliases for j in aliases[name]}
    if len(known_idx)<20:
        raise RuntimeError("Too few documented GloNAF hosts for observed GBIF species")
    hostmask=sum(1<<j for j in sorted(known_idx))
    allflora=[source_masks[.05][k] for k in selected_codes]
    hosts=[m&hostmask for m in allflora]
    nondocumented=[m&~hostmask for m in allflora]
    host_j,host_zero=ji(hosts,ia,ib)
    other_j,other_zero=ji(nondocumented,ia,ib)
    allalien_j,all_zero=ji(allflora,ia,ib)
    if all_zero:
        raise RuntimeError("Good GBIF regions had an undefined all-alien flora pair")

    flora_rich=np.asarray([m.bit_count() for m in allflora],float)
    eff=np.asarray([float(region_qa[k]["gbif_provider_days"]) for k in selected_codes])
    data_source=np.asarray([float(region_qa[k]["gbif_dataset_keys"]) for k in selected_codes])
    common=[x[q] for x in beta["control_sets"]["colonial_vs_resource_beta_after_flora"]]
    covars={
       "log_geographic_distance":common[0],
       "same_global_level1":common[1],
       "temperature_difference":common[2],
       "aridity_difference":common[3],
       "native_WCVP_flora_Jaccard":common[4],
       "native_239_potential_resource_Jaccard":common[5],
       "native_WCVP_flora_richness_min":common[6],
       "native_WCVP_flora_richness_max":common[7],
       "change_in_all_WCVP_flora_Jaccard":common[8],
       "four_empire_shared_duration":beta["exposure"][q],
       "GloNAF_all_alien_Jaccard":allalien_j,
       "GloNAF_all_alien_richness_min":np.log1p(np.minimum(flora_rich[ia],flora_rich[ib])),
       "GloNAF_all_alien_richness_max":np.log1p(np.maximum(flora_rich[ia],flora_rich[ib])),
       "GBIF_indexed_provider_days_min":np.log1p(np.minimum(eff[ia],eff[ib])),
       "GBIF_indexed_provider_days_max":np.log1p(np.maximum(eff[ia],eff[ib])),
       "GBIF_indexed_dataset_count_min":np.log1p(np.minimum(data_source[ia],data_source[ib])),
       "GBIF_indexed_dataset_count_max":np.log1p(np.maximum(data_source[ia],data_source[ib])),
    }
    base=list(covars.values())
    basis_host=orthogonal_controls(base+[other_j])
    basis_other=orthogonal_controls(base+[host_j])
    host_partial=correlation(residual(host_j,basis_host),residual(observed_j,basis_host))
    other_partial=correlation(residual(other_j,basis_other),residual(observed_j,basis_other))
    if host_partial is None or other_partial is None:
        raise RuntimeError("No variation in host or other-alien predictor on independent GBIF outcome")
    models={
       "M0_general_flora_and_effort":base,
       "M1_known_host":base+[host_j],
       "M2_other_aliens":base+[other_j],
       "M3_both":base+[host_j,other_j],
    }
    scores,comparisons=evaluate_loio(observed_j,models,ia,ib,len(selected_codes),
                                   float(lock["evaluation"]["linear_ridge_lambda"]))
    world={code:level for code,level in
           json.loads(Path("data/frozen/wgsrpd3_level1_lookup_v0.1.json").read_text(encoding="utf8"))["units"].items()}
    level=Counter(world[code] for code in selected_codes)
    report={
        "schema":"chocho_GBIF31_observed_assemblage_HOSTS_GloNAF_same111regions_v1",
        "status":"INDEPENDENT_BUTTERFLY_GBIF_INDEXED_ASSEMBLAGE_PILOT_NOT_TRUE_OCCUPANCY",
        "method_lock_sha256":file_sha(LOCK),
        "source":{
            "GBIF_download_DOI":"10.15468/dl.pp5nc9",
            "GBIF_archive_sha256":file_sha(GBIF_ZIP),
            "GloNAF_Zenodo_MD5":"16f5585d2500676951a4b664bbf0ffde",
            "HOSTS_commit":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "official_TDWG_commit":"52da7828aba9d461dd133c27b3bd7a4407161f54",
            "original_239_frozen_butterfly_DOI_species_count":30,
        },
        "cohort":{
            "eligible_wgsrpd3_regions":len(selected_codes),
            "nonindependent_region_pairs":len(ia),
            "actual_recorded_butterfly_species":len(stored_species),
            "all_original_source_taxa":len(selected_taxa),
            "world_division_region_counts":dict(level),
            "GBIF_region_record_counts":counts,
            "unique_recorded_taxa_early_2010_2017":len(first_species),
            "unique_recorded_taxa_late_2018_2026":len(second_species),
            "pairs_both_no_early_records":early_missing,
            "pairs_both_no_late_records":late_missing,
            "original_documented_focal_host_names":len(native_host_names),
            "host_names_matched_GloNAF":len(matched_aliases),
            "host_raw_GloNAF_taxon_IDs":len(known_idx),
            "pairs_both_no_documented_host_record":host_zero,
            "pairs_both_no_other_alien_record":other_zero,
        },
        "description":{
            "mean_recorded_butterfly_community_Jaccard":float(observed_j.mean()),
            "median_recorded_butterfly_community_Jaccard":float(np.median(observed_j)),
            "mean_known_host_GloNAF_Jaccard":float(host_j.mean()),
            "mean_other_GloNAF_Jaccard":float(other_j.mean()),
            "mean_all_GloNAF_Jaccard":float(allalien_j.mean()),
            "mean_earlier_2010_2017_recorded_community_Jaccard":float(first_j.mean()),
            "mean_later_2018_2026_recorded_community_Jaccard":float(second_j.mean()),
        },
        "adjusted_partial_rank":{
            "known_focal_larval_host_Jaccard_given_all_flora_geography_climate_and_effort":float(host_partial),
            "other_not_documented_aliens_Jaccard_given_known_hosts_and_controls":float(other_partial),
            "not_a_causal_or_iid_pair_p_value":True,
        },
        "region_held_out_ridge":{
            "definition":"train 110x109/2=5995 pairs and score 110 pairs sharing the unseen 111th region; equal weight per 111 regions",
            "ridge_lambda":float(lock["evaluation"]["linear_ridge_lambda"]),
            "models":scores,
            "incremental_vs_general_flora_and_effort":comparisons,
        },
        "limits":[
            "Only original GBIF occurrence DOI's 30 exact matches, not globally representative 239 or 225 butterfly species.",
            "Most 111 eligible regions are in the Americas; no claims of worldwide butterfly-community generality.",
            "GBIF presence is an indexed detection and species absence is not observed.",
            "Eligibility determined from all-period GBIF records, so leave-region-out validation is cross-sectional, not a prospective future test.",
            "Plant and host records are undated contemporary inventories; no causal host introduction or colonial-era temporal mechanism.",
            "GBIF provider/day counts partially proxy effort but cannot correct all observer biases or species detection probability.",
            "GBIF sampling is taxonomically selected and the region-pair observations are dependent; do not report naïve independent-pair p-values.",
            "HOSTS taxonomic matches remain source lower bounds, non-documented taxa are not proven biological nonhosts.",
        ],
    }
    # Numeric pair-level aggregate permits source review without exposing
    # original GBIF coordinates, gbifID or observer identity.
    pair_rows=[]
    for k,(a,b) in enumerate(zip(ia,ib)):
        pair_rows.append({
            "wgsrpd3_a":selected_codes[int(a)],"wgsrpd3_b":selected_codes[int(b)],
            "observed_butterfly_Jaccard":float(observed_j[k]),
            "general_GloNAF_alien_Jaccard":float(allalien_j[k]),
            "known_focal_HOSTS_Jaccard":float(host_j[k]),
            "other_alien_Jaccard":float(other_j[k]),
            "gbif_earlier_butterfly_Jaccard":float(first_j[k]),
            "gbif_later_butterfly_Jaccard":float(second_j[k]),
        })
    return report,pair_rows


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT)
    p.add_argument("--pairs",type=Path,default=PAIR_CSV)
    args=p.parse_args()
    out,pairrows=analyse()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with args.pairs.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(pairrows[0]))
        w.writeheader();w.writerows(pairrows)
    print(json.dumps({
        "status":out["status"],
        "source":out["source"],
        "cohort":out["cohort"],
        "partial":out["adjusted_partial_rank"],
        "heldout":out["region_held_out_ridge"],
    },indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

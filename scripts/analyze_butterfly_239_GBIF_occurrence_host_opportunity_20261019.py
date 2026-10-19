#!/usr/bin/env python3
"""Source-conditional global observed butterfly records vs host-resource geography.

GBIF responses are bounded, systematically paged and not representative
downloads. Comparison is association between *sampled presence records* and
potential larval-host resources under other-butterfly sampling-effort proxies,
butterfly-species and continent, plus climate-band sensitivity controls.

No simulated butterfly observations and no inference of true absence,
butterfly colonisation, egg production or invasive-host fitness.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np

from analyze_bce_augmented_plant_network_climate_null_20261019 import (
    source_augmented_model,analog_mask,
)
from analyze_butterfly_climate_matched_plant_range_null import (
    climate_strata,
)

WINDOWS=("2010_2017","2018_2025","pooled")
SOURCES=("all_accepted_bases","human_and_preserved_only")
EFFORT_GATES=(1,3,5)
QS=(.05,.10)
TAXON_SOURCES={"HUMAN_OBSERVATION","PRESERVED_SPECIMEN"}
PRIMARY_ARGS=("all_accepted_bases","pooled",3,.05,"level1")
SEED=20261019


def source_observation_counts(path:Path,manifest:dict):
    windows={(s,w):Counter() for s in SOURCES for w in WINDOWS}
    known={r["species"]:r for r in manifest["species_receipts"]}
    assert len(known)==manifest["species_requested"]
    assert set(WINDOWS[:2])=={k for k,_ in manifest["queries"]["windows"]}
    with path.open(encoding="utf-8",newline="") as f:
        rows=csv.DictReader(f)
        required={"species","window","wgsrpd3","basis","sampled_record_count"}
        assert set(rows.fieldnames or ())>=required
        seen=set()
        for r in rows:
            sp=r["species"];w=r["window"];geo=r["wgsrpd3"];basis=r["basis"]
            if sp not in known or w not in WINDOWS[:2]:
                raise ValueError("GBIF source metadata and aggregates disagree")
            if known[sp]["status"]!="COMPLETE":
                raise ValueError("Incomplete sampled species incorrectly supplied as observed")
            key=(sp,w,geo,basis)
            if key in seen:
                raise ValueError("Duplicate original aggregated row")
            seen.add(key)
            count=int(r["sampled_record_count"])
            if count<=0:raise ValueError("Nonpositive occurrence count")
            for source in SOURCES:
                if source=="all_accepted_bases" or basis in TAXON_SOURCES:
                    windows[(source,w)][sp,geo]+=count
                    windows[(source,"pooled")][sp,geo]+=count
    return windows


def climate_and_plant_candidates(model,regions,climate):
    assert len(model["species"])==239
    masks,eligible=analog_mask(model,regions,climate)
    assert len(eligible)==235
    index={code:i for i,code in enumerate(regions)}
    bylevel=model["level1"]
    sys_climate={
        "bio1":{i:float(climate[regions[i]][0]) for i in range(len(regions))},
        "bio12":{i:float(climate[regions[i]][2]) for i in range(len(regions))},
    }
    classes,metadata=climate_strata(
        climate,regions,bylevel,"temperature_precipitation_2x2")
    assert len(classes)==365
    introduced={
        h:sum(1<<k for k in cells)
        for h,cells in model["introduced_sets"].items()
    }
    added=[]
    for i,hosts in enumerate(model["species_hosts"]):
        union=0
        for h in hosts:
            union |= introduced.get(h,0)
        added.append(union & ~model["native_species"][i])
    assert sum(bits.bit_count() for bits in added)==14965
    by_sp={model["species"][i]:{
        "i":i,"genus":model["species"][i].split()[0],
        "native":model["native_species"][i],"added":added[i],
        "masks":{q:masks[q].get(i) for q in QS}
    } for i in range(239)}
    return by_sp,classes,index,metadata


def group_results(model,regions,manifest,counts,*,nb=999,main_n=9999,
                  target_species=None, effort_reference_species=None):
    records={x["species"]:x for x in manifest["species_receipts"]}
    species={s for s,r in records.items() if r["status"]=="COMPLETE"}
    # The source-selected target cohort and observer-effort background are
    # different populations. A recovered species may be evaluated as a
    # target without being silently added to the older 135-species effort
    # proxy; conversely, its own presence must still enter the target outcome.
    target=set(species) if target_species is None else set(target_species)
    effort=set(species) if effort_reference_species is None else set(effort_reference_species)
    if not target<=species or not effort<=species:
        raise ValueError("Original target/background GBIF source cohort is incomplete")
    bysp,classes,index,clmeta=climate_and_plant_candidates(
        model["model"],regions,model["climate"])
    species&=set(bysp)
    target &= species
    effort &= species
    level1=model["model"]["level1"]
    eligible={s for s in target if bysp[s]["masks"][.05] is not None}
    results={}
    for source in SOURCES:
        for window in WINDOWS:
            occ=counts[(source,window)]
            observed_by_region=defaultdict(set)
            effort_by_region=defaultdict(set)
            for (sp,code),c in occ.items():
                if sp in species and c>0 and code in index:
                    observed_by_region[code].add(sp)
                    if sp in effort:
                        effort_by_region[code].add(sp)
            for gate in EFFORT_GATES:
                for q in QS:
                    for stratification in ("level1","level1_bio1_bio12_2x2"):
                        block=defaultdict(lambda:{"N":0,"K":0,"A":0,"O":0})
                        eligible_species=set()
                        mapped_species={sp for sp in eligible if bysp[sp]["masks"][q] is not None}
                        # Species with no GBIF positive still contribute candidate
                        # cells but not a phantom observed presence.
                        for sp in sorted(mapped_species):
                            s=bysp[sp];native=s["native"];add=s["added"]
                            m=s["masks"][q]
                            for j,code in enumerate(regions):
                                if native>>j&1 or not m>>j&1:continue
                                other_effort=len(effort_by_region.get(code,set())-{sp})
                                if other_effort<gate:continue
                                g=str(level1[j])
                                if stratification!="level1":
                                    g+="|"+str(classes[j])
                                group=block[(sp,g)]
                                group["N"]+=1
                                group["K"]+=int(sp in observed_by_region.get(code,set()))
                                group["A"]+=int(bool(add>>j&1))
                                group["O"]+=int(bool(add>>j&1) and
                                                sp in observed_by_region.get(code,set()))
                        valid=[]
                        for (sp,g),b in sorted(block.items()):
                            if not (0<b["A"]<b["N"]):continue
                            if not 0<=b["O"]<=min(b["A"],b["K"]):
                                raise ValueError("Invalid observed overlap")
                            item={
                                "species":sp,"genus":bysp[sp]["genus"],
                                "stratum":g,**b,
                                "expected":b["A"]*b["K"]/b["N"],
                            }
                            valid.append(item)
                            eligible_species.add(sp)
                        totalN=sum(x["N"] for x in valid)
                        totalK=sum(x["K"] for x in valid)
                        totalA=sum(x["A"] for x in valid)
                        totalO=sum(x["O"] for x in valid)
                        expect=sum(x["expected"] for x in valid)
                        informative=sum(x["K"]>0 and x["A"]>0
                                        and x["A"]<x["N"] for x in valid)
                        enough=(len(eligible_species)>=25 and
                                informative>=10 and totalK>=25 and
                                totalA>=25)
                        label=f"{source}|{window}|E{gate}|q{q:.2f}|{stratification}"
                        out={
                            "status":("EXPLORATORY_SAMPLED_PRESENCE_OVERLAP" if enough
                                      else "NOT_IDENTIFIABLE_SAMPLED_SPECIES_OR_REGIONS"),
                            "source_basis":source,"window":window,
                            "minimum_other_observed_butterfly_species":gate,
                            "climate_analogue_lower_quantile":q,
                            "strata":stratification,
                            "source_species_complete":len(species),
                            "focal_target_species_eligible_for_source":len(target),
                            "effort_reference_species_eligible_for_source":len(effort),
                            "species_native_climate_eligible":len(mapped_species),
                            "species_with_effort_supported_two_resource_classes":len(eligible_species),
                            "species_with_sampled_positive_in_informative_strata":
                                len({x["species"] for x in valid if x["K"]>0}),
                            "eligible_species_level1_strata":len(valid),
                            "informative_positive_strata":informative,
                            "candidate_region_species_n":totalN,
                            "sampled_positive_region_species_n":totalK,
                            "candidate_added_host_cells_n":totalA,
                            "observed_positives_on_added_host_cells_n":totalO,
                            "continent_conditional_expected_overlap":expect,
                            "observed_minus_expected_cells":totalO-expect,
                            "observed_fraction_positive_in_added_cells":totalO/totalK if totalK else None,
                            "expected_fraction_positive_in_added_cells":expect/totalK if totalK else None,
                            "observed_minus_expected_pp":100*(totalO-expect)/totalK if totalK else None,
                            "overlap_observed_to_expected_ratio":totalO/expect if expect else None,
                            "sampling_effort_proxy_regions_with_any_other_focal_species":
                                len([v for v in effort_by_region.values() if v]),
                            "null":"Within each species×continent(/coarse-climate) conditioned on observed positive sampled regions and number of added host regions, use hypergeometric uniform placement among focal-other-butterfly effort-supported native-unsupported candidate regions. This is NOT a causal/true-absence test.",
                        }
                        if enough:
                            rng=np.random.default_rng(
                                SEED+len(valid)+int(q*100)+gate*71+
                                sum(map(ord,label))%1297)
                            B=main_n if (source,window,gate,q,stratification)==PRIMARY_ARGS else nb
                            sim=np.zeros(B,dtype=np.int32)
                            for x in valid:
                                sim+=rng.hypergeometric(x["K"],x["N"]-x["K"],
                                                       x["A"],size=B)
                            out["null_draws"]=B
                            out["null_sim_q025_q500_q975"]=[
                                float(np.quantile(sim,z)) for z in (.025,.5,.975)
                            ]
                            out["conditional_geographical_null_right_tail_not_biological_p"] = (
                                int(np.sum(sim>=totalO))+1)/(B+1)
                            bygen=defaultdict(lambda:{"obs":0.,"exp":0.,"pos":0.})
                            for x in valid:
                                gg=bygen[x["genus"]]
                                gg["obs"]+=x["O"];gg["exp"]+=x["expected"];gg["pos"]+=x["K"]
                            genera=sorted(bygen)
                            boots=[]
                            Bboot=2999 if B==main_n else 499
                            for _ in range(Bboot):
                                picks=rng.integers(0,len(genera),len(genera))
                                a=sum(bygen[genera[i]]["obs"] for i in picks)
                                e=sum(bygen[genera[i]]["exp"] for i in picks)
                                pp=sum(bygen[genera[i]]["pos"] for i in picks)
                                if pp>0:boots.append(100*(a-e)/pp)
                            out["bootstrap_genus_clusters"]=len(genera)
                            out["bootstrap_draws_valid"]=len(boots)
                            out["genus_bootstrap_excess_pp_95ci"]=[
                                float(np.quantile(boots,z)) for z in (.025,.975)
                            ] if boots else None
                        results[label]=out
    return results


def run(args):
    receipt=json.loads(args.gbif_manifest.read_text(encoding="utf-8"))
    if receipt.get("schema")!="chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01":
        raise ValueError("Invalid external GBIF sample provenance")
    if not (1<=receipt["species_requested"]<=239):
        raise ValueError("GBIF panel size not audited")
    # Fail closed on full panel under strict request, but do NOT promote
    # unread/missing butterfly species to negative/absence evidence.
    if args.require_all_species and receipt["species_requested"]!=239:
        raise ValueError("Cannot publish 239-panel validation from smoke test")
    observations=source_observation_counts(args.gbif_counts,receipt)
    model,regions,climate,prov=source_augmented_model(args.source_root)
    grouped=group_results({"model":model,"climate":climate},regions,receipt,observations,
                          nb=args.null_draws_sensitivity,main_n=args.null_draws_primary)
    main_key="all_accepted_bases|pooled|E3|q0.05|level1"
    primary=grouped[main_key]
    background_reference=None
    if args.original_gbif_manifest is not None:
        original=json.loads(args.original_gbif_manifest.read_text(encoding="utf-8"))
        original_species={r["species"] for r in original["species_receipts"]
                          if r["status"]=="COMPLETE"}
        if len(original_species)!=135:
            raise ValueError("Original GBIF background 135 source cohort has drifted")
        current_species={r["species"] for r in receipt["species_receipts"]
                         if r["status"]=="COMPLETE"}
        if not original_species<=current_species:
            raise ValueError("Recovered GBIF receipt lost original complete source species")
        source_args={"model":model,"climate":climate}
        # 1: source invariable original 135 target+background: must exactly
        #     reproduce the original 135 study, even in a merged 239 receipt.
        old=group_results(source_args,regions,receipt,observations,
            nb=args.null_draws_sensitivity,main_n=args.null_draws_primary,
            target_species=original_species,effort_reference_species=original_species)
        old_primary=old[main_key]
        if (old_primary["sampled_positive_region_species_n"]!=119 or
            old_primary["observed_positives_on_added_host_cells_n"]!=74 or
            abs(old_primary["continent_conditional_expected_overlap"]
                -69.61846219208861)>1e-7):
            raise ValueError("Original 135-butterfly GBIF ecological baseline changed")
        # 2: newly completed species are allowed to contribute target GBIF
        #     presences while keeping the background effort panel unchanged.
        kept=group_results(source_args,regions,receipt,observations,
            nb=args.null_draws_sensitivity,main_n=args.null_draws_primary,
            effort_reference_species=original_species)
        background_reference={
            "original_initial135_source_invariant_primary":old_primary,
            "expanded_targets_fixed_original135_effort_primary":kept[main_key],
            "expanded_targets_fixed_original135_effort_all_combinations":kept,
            "interpretation":"Posthoc sampling-denominator diagnostic: an expanded butterfly taxon pool also expands E_-i regional effort proxy; this separates those changes. No detection or colonisation effect is estimated."
        }
    source_diagnostics={
        "gbif_species_requested":receipt["species_requested"],
        "gbif_resolved_complete":receipt["species_taxon_resolved_and_complete"],
        "gbif_unresolved_species":sorted(
            x["species"] for x in receipt["species_receipts"]
            if x["status"]!="COMPLETE"),
        "gbif_sampled_mapped_records":receipt["outputs"]["source_records_mapped"],
        "gbif_sampled_region_basis_rows":receipt["outputs"]["aggregated_region_basis_rows"],
        "taxon_query_saturation_count":sum(
            any(w.get("total_search_more_than_100000",False)
                for w in x.get("windows",{}).values())
            for x in receipt["species_receipts"]),
        "sample_provenance_source":"GBIF occurrence REST API bounded deterministic ordinal pages, two 8-year windows, not an archived DOI download",
    }
    return {
        "schema":"chocho_239_butterflies_effort_climate_matched_gbif_presence_macro_v0.1",
        "status":primary["status"],
        "source_scope":"239 frozen butterfly host-resource panel with partial BCE accepted host catalogue (author-only WCVP/HOSTS); two bounded 2010-2025 GBIF sampling windows",
        "analysis_date":"2026-10-19",
        "original_source":prov,
        "GBIF_sample_metadata":source_diagnostics,
        "primary":primary,
        "all_combinations":grouped,
        "fixed_original_135_background_diagnostic":background_reference,
        "no_inference_gates":[
            "No true absence: GBIF sample absence is not an ecological detection/non-detection survey",
            "No adult breeding, larval use, egg allocation or butterfly population recruitment",
            "No causal plant introduction dates, changing GBIF data over time invalidates colonisation inference",
            "Unequal WGSRPD3 land area and background biological observer effort unmeasured",
            "Other 239 butterfly species' distinct sample coverage is a partial effort PROXY, not survey effort",
            "Year windows are comparisons of historical GBIF sample coverage, not independent landscape experiments",
            "Systematic sampled ordinal GBIF result pages are not a random occurrence or complete GBIF download"
        ]
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--gbif-counts",type=Path,required=True)
    p.add_argument("--gbif-manifest",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    p.add_argument("--null-draws-primary",type=int,default=9999)
    p.add_argument("--null-draws-sensitivity",type=int,default=999)
    p.add_argument("--require-all-species",action="store_true")
    p.add_argument("--original-gbif-manifest",type=Path,
                   help="Original 135-complete GBIF receipt to separate target taxa from internal effort-panel expansion")
    args=p.parse_args()
    result=run(args)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",
                                encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "gbif":result["GBIF_sample_metadata"],
        "primary":result["primary"],
        "sensitivity":{k:{
            "status":v["status"],
            "positive_n":v["sampled_positive_region_species_n"],
            "observed_minus_expected_pp":v["observed_minus_expected_pp"],
            "null_tail":v.get("conditional_geographical_null_right_tail_not_biological_p")
        } for k,v in result["all_combinations"].items()}
    },sort_keys=True))

#!/usr/bin/env python3
"""Independent 239-butterfly HOSTS vs partial BCE accepted-host macroclimate audit.

This is a source-completeness robustness check on a global potential
host-resource geography analysis. It is NOT a new published plant or
butterfly field observation, and not independent global host-catalogue
replication: BCE candidate pairs are Europe-biased.

All WCVP/HOSTS raw inputs are from the author-only SHA-verified offline
source archive; no plant×region source data are emitted. Results are
butterfly-level counts and source-qualified aggregate comparisons.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from analyze_butterfly_climate_analog_macro_20261019 import source_model, qbounds

QS=(0.05,0.10)
EXPECTED_BASE={"native":26530,"contemporary":41083,"added":14553}
EXPECTED_AUG={"native":29852,"contemporary":44819,"added":14967}


def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as handle:
        return list(csv.DictReader(handle))


def reconstruct(source_root):
    source_root=Path(source_root)
    # Guarantees exact original 239 × 365 resource indices, original 26,530 /
    # 41,083 / 14,553 source totals and SHA integrity of archive members.
    model,regions,climate,orig,source_prov=source_model(source_root)
    sys.path.insert(0,str(source_root/"code"))
    from analyze_bce_resource_novelty_source_sensitivity import (
        load_candidate_host_links, host_geography, combine_host_ranges
    )
    data=source_root/"inputs"
    candidate,metadata=load_candidate_host_links(data/"bce/accepted_host_edges.json")
    assert len(candidate)==81 and sum(map(len,candidate.values()))==1027
    original_metric={r["species"]:r for r in read_csv(
        data/"anthropogenic_species_metrics.csv")}
    assert len(original_metric)==239
    hosts=defaultdict(set)
    for row in read_csv(data/"original/insect_host_accepted.csv"):
        sp=row["insect_species"].strip()
        host=row["accepted_plant_name_id"].strip()
        if sp in original_metric and host:
            hosts[sp].add(host)
    assert len(hosts)==239 and all(hosts[s] for s in original_metric)
    assert set(candidate).issubset(hosts)
    assert not any(hosts[s]&candidate[s] for s in candidate)
    extra=set().union(*candidate.values())
    native0=host_geography(data/"original/native_extant_nondoubtful_wgsrpd3.csv")
    cont0=host_geography(data/"original/extant_nondoubtful_wgsrpd3.csv")
    native_plus=combine_host_ranges(native0,
        host_geography(data/"bce/additional_native_wgsrpd3.csv"),extra)
    cont_plus=combine_host_ranges(cont0,
        host_geography(data/"bce/additional_contemporary_wgsrpd3.csv"),extra)
    level1=json.loads((data/"level1_lookup.json").read_text())["units"]
    assert set(regions)<=set(level1)
    assert set(regions)<=set(climate)
    # The 365-region climate crossfit is the fixed *analysis universe*;
    # accepted BCE hosts can introduce geographic codes beyond that grid.
    # In this version ANT and HMD are in introduced resources for one
    # butterfly, without source climate vectors. Keep geographic source
    # totals untouched, then transparently condition ONLY the climate
    # analysis on all climate-covered added cells (not on an entire species).
    coverage=set(climate)
    assert len(coverage)==365
    assert set(regions)==coverage
    maps={}
    for name,augmented,nh,ch in (
        ("original_HOSTS",False,native0,cont0),
        ("BCE_partial_augmented",True,native_plus,cont_plus),
    ):
        rows=[]
        for sp in sorted(original_metric):
            links=hosts[sp] | (candidate.get(sp,set()) if augmented else set())
            native=set().union(*(nh.get(h,set()) for h in links))
            contemporary=set().union(*(ch.get(h,set()) for h in links))
            assert native<=contemporary
            newly=contemporary-native
            n=len(native);c=len(contemporary);a=len(newly)
            if not augmented:
                old=original_metric[sp]
                assert (n,c,a)==(int(old["native_resource_units"]),
                                  int(old["contemporary_resource_units"]),
                                  int(old["introduced_added_units"])),sp
            natcl=sorted(native & coverage)
            covered_additions=newly & coverage
            missing_additions=newly-coverage
            row={
                "species":sp,"genus":sp.split(" ")[0],
                "eligible_native_climate_regions_ge_10":len(natcl)>=10,
                "has_new_BCE_accepted_host_links":sp in candidate,
                "native_regions":n,"contemporary_regions":c,
                "newly_supported_regions":a,
                "newly_supported_climate_evaluable_regions":len(covered_additions),
                "added_regions_missing_climate":len(missing_additions),
                "missing_climate_region_codes":sorted(missing_additions),
                "climate_analog":{},
            }
            if len(natcl)>=10:
                for q in QS:
                    tlo,thi=qbounds([climate[x][0] for x in natcl],q)
                    plo,phi=qbounds([climate[x][2] for x in natcl],q)
                    analog={code for code,value in climate.items()
                            if tlo<=value[0]<=thi and plo<=value[2]<=phi}
                    O=len(covered_additions&analog)
                    E=0.0
                    for broad in set(level1.values()):
                        candidate_space={
                            code for code in regions
                            if level1[code]==broad and code not in native
                        }
                        if not candidate_space:continue
                        add_in_area=len(covered_additions&candidate_space)
                        analog_in_area=len(analog&candidate_space)
                        E+=add_in_area*analog_in_area/len(candidate_space)
                    covered_A=len(covered_additions)
                    assert 0<=O<=covered_A and 0<=E<=covered_A+1e-8
                    row["climate_analog"][str(q)]={
                        "added":covered_A,
                        "aligned":O,
                        "expected_continent_uniform":E,
                        "aligned_fraction":O/covered_A if covered_A else None,
                        "excess_fraction":(O-E)/covered_A if covered_A else None,
                    }
            rows.append(row)
        totals={
            "native":sum(r["native_regions"] for r in rows),
            "contemporary":sum(r["contemporary_regions"] for r in rows),
            "added":sum(r["newly_supported_regions"] for r in rows),
        }
        assert totals== (EXPECTED_AUG if augmented else EXPECTED_BASE),(name,totals)
        missing={
            row["species"]:row["missing_climate_region_codes"]
            for row in rows if row["added_regions_missing_climate"]>0
        }
        maps[name]={
            "rows":rows,"totals":totals,
            "climate_source_coverage":{
                "source_newly_supported_all":totals["added"],
                "newly_supported_with_climate":sum(
                    r["newly_supported_climate_evaluable_regions"] for r in rows),
                "newly_supported_without_climate":sum(
                    r["added_regions_missing_climate"] for r in rows),
                "butterfly_species_affected":len(missing),
                "missing_region_codes_by_species":missing
            }
        }
    return maps,metadata,source_prov


def aggregate(rows,*,selected=None,q=0.05):
    rows=[r for r in rows if r["eligible_native_climate_regions_ge_10"]
          and (selected is None or r["species"] in selected)]
    evaluated=[r for r in rows if r["newly_supported_climate_evaluable_regions"]>0]
    key=str(q)
    O=sum(r["climate_analog"][key]["aligned"] for r in evaluated)
    E=sum(r["climate_analog"][key]["expected_continent_uniform"] for r in evaluated)
    A=sum(r["newly_supported_climate_evaluable_regions"] for r in evaluated)
    if A<=0:raise ValueError("no evaluable newly introduced host regions")
    return {
        "eligible_species":len(rows),
        "expanded_eligible_species":len(evaluated),
        "newly_supported_region_cells":A,
        "observed_analog_cells":O,
        "continent_uniform_expected_analog_cells":E,
        "observed_analog_fraction":O/A,
        "expected_analog_fraction":E/A,
        "observed_minus_expected_fraction":(O-E)/A,
        "expanded_species_with_positive_enrichment":sum(
            r["climate_analog"][key]["excess_fraction"]>0 for r in evaluated),
        "expanded_species_with_negative_enrichment":sum(
            r["climate_analog"][key]["excess_fraction"]<0 for r in evaluated),
        "expanded_species_with_exact_zero":sum(
            r["climate_analog"][key]["excess_fraction"]==0 for r in evaluated),
    }


def run(source):
    maps,metadata,provenance=reconstruct(source)
    old={r["species"]:r for r in maps["original_HOSTS"]["rows"]}
    new={r["species"]:r for r in maps["BCE_partial_augmented"]["rows"]}
    assert len(old)==len(new)==239
    affected=set(k for k in old if new[k]["has_new_BCE_accepted_host_links"])
    assert len(affected)==81
    results={}
    for q in QS:
        tag=str(q)
        a=aggregate(list(old.values()),q=q)
        b=aggregate(list(new.values()),q=q)
        o=aggregate(list(old.values()),selected=affected,q=q)
        t=aggregate(list(new.values()),selected=affected,q=q)
        # Source sensitivity is a paired butterfly×catalogue comparison,
        # not a butterfly temporal transition or independent global sample.
        shared=[
            (old[s]["climate_analog"][tag],new[s]["climate_analog"][tag],s)
            for s in sorted(affected)
            if old[s]["eligible_native_climate_regions_ge_10"]
            and new[s]["eligible_native_climate_regions_ge_10"]
            and old[s]["newly_supported_regions"]>0
            and new[s]["newly_supported_regions"]>0
        ]
        differences=np.array([y["excess_fraction"]-x["excess_fraction"]
                             for x,y,_ in shared],dtype=float)
        result={
            "original_all_239_source_eligible":a,
            "partially_augmented_all_239_source_eligible":b,
            "original_81_affected_subset":o,
            "augmented_81_affected_subset":t,
            "paired_affected_species_positive_added_both_n":len(shared),
            "paired_affected_species_mean_change_in_excess_fraction":float(np.mean(differences)),
            "paired_affected_species_median_change_in_excess_fraction":float(np.median(differences)),
            "paired_affected_species_number_enrichment_increased":int(np.sum(differences>1e-12)),
            "paired_affected_species_number_enrichment_decreased":int(np.sum(differences<-1e-12)),
            "paired_affected_species_number_enrichment_unchanged":int(np.sum(abs(differences)<=1e-12)),
            "global_aggregate_excess_pp_difference_augmented_minus_original":
                100*(b["observed_minus_expected_fraction"]-
                     a["observed_minus_expected_fraction"]),
        }
        results[tag]=result
    orig=results["0.05"]["original_all_239_source_eligible"]
    assert orig["eligible_species"]==234
    assert orig["expanded_eligible_species"]==204
    assert orig["newly_supported_region_cells"]==14481
    assert orig["observed_analog_cells"]==9451
    assert abs(orig["continent_uniform_expected_analog_cells"]-8251.908515093264)<1e-5
    base_strict=results["0.1"]["original_all_239_source_eligible"]
    assert base_strict["observed_analog_cells"]==6948
    assert abs(base_strict["continent_uniform_expected_analog_cells"]-5834.664323099673)<1e-5
    aug=results["0.05"]["partially_augmented_all_239_source_eligible"]
    assert aug["eligible_species"]==235
    assert aug["expanded_eligible_species"]==212
    assert aug["newly_supported_region_cells"]==14961
    assert aug["observed_analog_cells"]==9949
    assert abs(aug["continent_uniform_expected_analog_cells"]-8712.829707575029)<1e-5
    second_aug=results["0.1"]["partially_augmented_all_239_source_eligible"]
    assert second_aug["newly_supported_region_cells"]==14961
    assert second_aug["observed_analog_cells"]==7370
    assert abs(second_aug["continent_uniform_expected_analog_cells"]-6119.869308104752)<1e-5
    assert maps["original_HOSTS"]["climate_source_coverage"]["newly_supported_without_climate"]==0
    aud=maps["BCE_partial_augmented"]["climate_source_coverage"]
    assert aud["newly_supported_without_climate"]==2
    assert aud["missing_region_codes_by_species"]=={
        "Coenonympha oedippus":["ANT","HMD"]
    }
    return {
        "schema":"chocho_butterfly_climate_analogue_BCE_partial_source_sensitivity_v0.1",
        "status":"SOURCE_RECONSTRUCTED_EXPLORATORY_MACRO_CLIMATE_SENSITIVITY",
        "original_source_SHA256":provenance["source_frozen_manifest_sha256"],
        "additional_hosts_source_SHA256":metadata["bce_accepted_edges_sha256"],
        "candidate_host_pairs":1027,
        "butterflies_with_partially_augmented_hosts":81,
        "unique_newly_added_accepted_plant_IDs":679,
        "whole_panel_resource_totals":{
            name:obj["totals"] for name,obj in maps.items()
        },
        "climate_source_coverage":{
            name:obj["climate_source_coverage"] for name,obj in maps.items()
        },
        "comparisons_by_native_host_climate_quantile":results,
        "claim_limitations":[
            "European-biased BCE/Clarke supplementary host links; NOT an independent globally representative second dataset.",
            "Changing reference native-host resources can classify a formerly newly introduced region as native-supported. These are host-catalogue differences, NOT temporal colonisation or plant loss.",
            "The null fixes each butterfly's newly supported region counts per broad continent and uniform allocation within Level1, NOT the introduced plant×region network.",
            "Climate estimands use exact 365 region vectors; two BCE-introduced cells (ANT, HMD) without climate are retained in geographic totals but excluded at the cell level from climate denominators. The corresponding Coenonympha oedippus has 60 climate-evaluable newly supported cells and is NOT discarded.",
            "All numbers are potential recorded larval host-resource regions rather than occupied habitats, adult recruitment, or fitness.",
            "Do not publicly redistribute author-only WCVP/HOSTS/TDWG source archive without rights approval.",
            "The ecology hypothesis was explored after source results were inspected. This is source-conditional validation, not a confirmatory independent test."
        ]
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    args=p.parse_args()
    receipt=run(args.source_root)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":receipt["status"],
        "candidate_host_pairs":receipt["candidate_host_pairs"],
        "whole_panel_resource_totals":receipt["whole_panel_resource_totals"],
        "thresholds":{
            k:{
                "original":v["original_all_239_source_eligible"],
                "augmented":v["partially_augmented_all_239_source_eligible"],
                "BCE_81_modified_original":v["original_81_affected_subset"],
                "BCE_81_modified_augmented":v["augmented_81_affected_subset"],
                "paired_species":v["paired_affected_species_positive_added_both_n"],
                "mean_paired_change":v["paired_affected_species_mean_change_in_excess_fraction"],
            }
            for k,v in receipt["comparisons_by_native_host_climate_quantile"].items()
        }
    },sort_keys=True))

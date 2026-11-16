#!/usr/bin/env python3
"""Within-butterfly taxon matched HOSTS literature-geography comparison.

Question: When the SAME named butterfly is reported with known larval hosts
on both islands and mainland, is the island report's introduced-host fraction
or host-family flora-adjusted fraction elevated?

This is a source-composition diagnostic, NOT a measured within-species diet
switch, habitat use, survival, or a wild-feeding experiment.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from analyze_HOSTS_reported_use_island_fao_ports_20261116 import build_rows
from analyze_HOSTS_vs_WCVP_family_flora_baseline_20261116 import MODES, baseline, read, select

SEED = 20261116


def species_comparison(obs):
    """Each insect identity is one unit, each source locality equally weighted."""
    by_species = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    distinct_pairs = defaultdict(lambda:defaultdict(lambda:defaultdict(set)))
    sources = set()
    for r in obs:
        group = "island" if r["is_island"] else "mainland"
        id_ = r["species"]
        loc = r["location"]
        by_species[id_][group][loc].append(r)
        distinct_pairs[id_][group][loc].add(r["plantid"])
        sources.add(loc)

    paired = []
    for species, both in sorted(by_species.items()):
        if "island" not in both or "mainland" not in both:
            continue
        measures = {}
        for group in ("island","mainland"):
            by_loc = both[group]
            locality_means = []
            for loc,arr in by_loc.items():
                # These already are distinct insect × plant × locality records.
                if len(arr)!=len(distinct_pairs[species][group][loc]):
                    raise ValueError(f"duplicate original plant pair per locality: {species} {loc}")
                raw = float(np.mean([r["observed"] for r in arr]))
                expected = float(np.mean([r["expected_family"] for r in arr]))
                locality_means.append((raw,raw-expected,expected,len(arr),loc))
            measures[group] = {
                "localities":len(locality_means),
                "distinct_host_x_locality_pairs":sum(z[3] for z in locality_means),
                "raw":float(np.mean([z[0] for z in locality_means])),
                "family_excess":float(np.mean([z[1] for z in locality_means])),
                "family_expectation":float(np.mean([z[2] for z in locality_means])),
                "source_localities":sorted(z[4] for z in locality_means),
            }
        paired.append({
            "species":species,
            "island":measures["island"],
            "mainland":measures["mainland"],
            "difference_raw":measures["island"]["raw"]-measures["mainland"]["raw"],
            "difference_family_excess":measures["island"]["family_excess"]-measures["mainland"]["family_excess"],
            "difference_family_expectation":measures["island"]["family_expectation"]-measures["mainland"]["family_expectation"],
        })
    return paired


def source_pair_evidence(obs):
    """Same butterfly × accepted host indexed in both geographic groups."""
    pairs = defaultdict(lambda:{"island":[],"mainland":[]})
    for r in obs:
        g = "island" if r["is_island"] else "mainland"
        pairs[r["species"],r["plantid"]][g].append(r)
    common=[]
    for (sp,host),d in pairs.items():
        if not d["island"] or not d["mainland"]:
            continue
        ic={r["origin"] for r in d["island"]}
        mc={r["origin"] for r in d["mainland"]}
        if len(ic)>1 or len(mc)>1:
            # A globally shared host may have different origin across sampled
            # island/continent localities. Preserve the ambiguity explicitly.
            typ="multi_locality_origin_disagreement"
        elif ic=={"introduced"} and mc=={"native"}:
            typ="island_introduced_mainland_native"
        elif ic=={"native"} and mc=={"introduced"}:
            typ="island_native_mainland_introduced"
        elif ic==mc:
            typ="same_origin_both_regions"
        else:
            typ="other_origin_pattern"
        common.append({
            "butterfly":sp,"accepted_host_id":host,
            "host_name":d["island"][0]["plant"],
            "source_pattern":typ,
            "islands":sorted({r["location"] for r in d["island"]}),
            "mainlands":sorted({r["location"] for r in d["mainland"]}),
            "rows":len(d["island"])+len(d["mainland"]),
        })
    return common


def boot(x,labels):
    if len(x)<3:
        return {"status":"TOO_FEW_SHARED_BUTTERFLIES","n_shared_species":len(x)}
    r=np.random.default_rng(SEED)
    out={}
    for key in labels:
        v=np.array([q[key] for q in x])
        boots=np.mean(r.choice(v,size=(9999,len(v)),replace=True),axis=1)
        out[key]={
          "n_shared_butterflies":len(v),
          "mean_paired_difference":float(np.mean(v)),
          "median_paired_difference":float(np.median(v)),
          "bootstrap95_species_level_mean_descriptive":[float(z) for z in np.quantile(boots,[.025,.975])],
          "positive_species":int(sum(v>0)), "negative_species":int(sum(v<0)),
          "zero_species":int(sum(v==0)),
          "subject_to_location_and_publication_dependence":True,
        }
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("hosts","host-resolution","wcvp-host-origin","locations","fao",
                 "ports","areas","butterflies239","flora-family","out-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    source=SimpleNamespace(
      hosts=a.hosts,host_resolution=a.host_resolution,wcvp=a.wcvp_host_origin,
      locations=a.locations,fao=a.fao,ports=a.ports,areas=a.areas,butterflies239=a.butterflies239
    )
    original,locs,coverage,receipts=build_rows(source)
    fam_by_id=defaultdict(set)
    for r in read(a.host_resolution):
        if r["match_status"]=="matched":
            fam_by_id[r["accepted_plant_name_id"]].add(r["family"])
    flora=baseline(a.flora_family)
    scope_reports={}
    all_cases=[]
    for cohort in ("HOSTS_Papilionoidea","fixed239_butterflies"):
        for crop in MODES:
            rr=list(select(original,cohort,crop))
            # Match previous source's min10 population, counting usable
            # source rows *after* confirming taxonomic family coverage.
            usable=[]
            per_loc=defaultdict(list)
            for r in rr:
                code=r["codes"]
                ff=fam_by_id.get(r["accepted_plant_id"],set())
                if len(ff)!=1:continue
                family=next(iter(ff))
                region=flora.get((code,crop,family),[0,0])
                if sum(region)==0:continue
                item={
                  "species":r["insect_species"],
                  "plantid":r["accepted_plant_id"],
                  "plant":r["plant"],
                  "location":r["location"],
                  "is_island":r["comparison_group"]=="island",
                  "origin":r["origin"],
                  "observed":1 if r["origin"]=="introduced" else 0,
                  "expected_family":region[1]/sum(region),
                }
                per_loc[r["location"]].append(item)
            for loc,arr in per_loc.items():
                if len(arr)>=10:
                    usable.extend(arr)
            # De-duplicate literature observations by their inferential
            # insect × accepted plant × original location combination.
            independent_pairs={}
            for r in usable:
                key=(r["species"],r["plantid"],r["location"])
                if key in independent_pairs:
                    old=independent_pairs[key]
                    if old["origin"]!=r["origin"] or old["expected_family"]!=r["expected_family"]:
                        raise ValueError("Conflicting source origin for same species/host/locality")
                else:
                    independent_pairs[key]=r
            obs=list(independent_pairs.values())
            matched=species_comparison(obs)
            hostpairs=source_pair_evidence(obs)
            stats=boot(matched,("difference_raw","difference_family_excess","difference_family_expectation"))
            paircounts=defaultdict(int)
            for z in hostpairs:
                paircounts[z["source_pattern"]]+=1
            key=cohort+"__"+crop
            summary={
              "cohort":cohort,"crop":crop,
              "source_usable_raw_rows":len(usable),
              "unique_insect_host_original_location_records":len(obs),
              "retained_source_localities":len({r["location"] for r in obs}),
              "island_localities":len({r["location"] for r in obs if r["is_island"]}),
              "mainland_localities":len({r["location"] for r in obs if not r["is_island"]}),
              "species_in_source_localities":len({r["species"] for r in obs}),
              "shared_species_island_and_mainland":len(matched),
              "same_insect_accepted_host_indexed_both":len(hostpairs),
              "shared_host_pairs_pattern_counts":dict(paircounts),
              "paired_comparison":stats,
              "matched_species_examples":sorted(matched,key=lambda x:-abs(x["difference_raw"]))[:12],
              "source_matched_species_complete":matched,
              "source_shared_host_pairs":hostpairs,
            }
            scope_reports[key]=summary
            print("HOSTS_WITHIN_BUTTERFLY_PAIR_SCOPE",json.dumps({
                k:v for k,v in summary.items()
                if k not in ("matched_species_examples","source_matched_species_complete","source_shared_host_pairs")
            },ensure_ascii=False),flush=True)
    result={
      "schema":"chocho_same_butterfly_in_island_mainland_reported_host_20261117_v1",
      "status":"SOURCE_PAIRED_LITERATURE_REPORTS_NOT_CONFIRMED_WILD_HOST_SWITCHING",
      "source_original":receipts,
      "original_source_scope":coverage,
      "cohorts":scope_reports,
      "design":"Each original butterfly species with at least one taxonomically resolved host relation in an island source locality and at least one in a mainland source locality. Unique insect × accepted host × original locality relations, averaged within species then between groups, same-species paired differences.",
      "uncertainty":"Bootstrap original shared species identities, not geographical locality or reference literature; confounding and dependence remain. Descriptive only.",
      "limits":[
        "The same species reported from both mainland and island does not imply colonization timing, ancestry, or host switching.",
        "The source HOSTS Location is a bibliographic attribution, not ground-truthed wild feeding.",
        "The source sample is selected on having enough recorded plant-host literature, and does not represent all butterflies.",
        "The same host plant record on both shores is not evidence of a live plant being shipped between them.",
        "Different geographies have different availability and detection; the family-matched species richness is not plant abundance.",
        "FAO species/genera exclusions change eligible taxa and locations; do not treat nested analyses as independent replications.",
        "WCVP introduced flags describe current regional status and not that status at the historical HOSTS record date.",
        "Species-level bootstraps are not phylogenetically independent experimental replicates.",
      ],
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    output=a.out_dir/"same_butterfly_island_mainland_host_paired.json"
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("SOURCE_PAIRED_BUTTERFLY_AUDIT_FILE",output,flush=True)


if __name__=="__main__":
    main()

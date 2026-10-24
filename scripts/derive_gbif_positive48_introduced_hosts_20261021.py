#!/usr/bin/env python3
"""Strict WCVP/HOSTS/BCE introduced-host candidates for 48 frozen GBIF-positive country cells.

The author-only source archive is validated and reconstructed in a GitHub
workflow, never redistributed. Five focal *juvenile-label* cases are fixed
before examining botanical candidates; no inferred actual feeding.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


BASE="provenance/reviewer_defenses/results/gbif_country_20261020/"
EXPECTED_COHORT=(48,20,23,26)
FIXED_CASES=frozenset({
    ("Anteos clorinde","CU"),("Colobura dirce","VE"),
    ("Papilio demoleus","CU"),("Papilio demoleus","DO"),
    ("Papilio demoleus","HT"),
})
FIXED_HEAD="gbif_juvenile_gadm1_readout_20261021.json"
SCHEMA="chocho_host_plant_candidates_5_juvenile_added_country_cases_20261021_v1"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_country_links(path,regions):
    idx={r:i for i,r in enumerate(regions)}
    bycountry=defaultdict(lambda:{"material":set(),"guard":set()})
    with Path(path).open(newline="",encoding="utf-8") as handle:
        for r in csv.DictReader(handle):
            region=r["wgsrpd3"]; country=r["country"]
            if region not in idx:continue
            if r["material_link"]=="1":
                bycountry[country]["material"].add(idx[region])
            if r["guard_link"]=="1":
                bycountry[country]["guard"].add(idx[region])
    return bycountry


def plant_names(source_root):
    source_root=Path(source_root)
    names=defaultdict(set)
    path=source_root/"inputs/original/insect_host_accepted.csv"
    with path.open(newline="",encoding="utf8") as handle:
        for r in csv.DictReader(handle):
            pid=r.get("accepted_plant_name_id","").strip()
            name=r.get("accepted_name","").strip()
            if pid and name:names[pid].add(name)
    p=source_root/"inputs/bce/accepted_host_edges.json"
    details=json.loads(p.read_text(encoding="utf8"))
    for taxon in details["per_butterfly"]:
        for row in taxon["details"]:
            pid=row.get("accepted_plant_name_id")
            name=row.get("accepted_name")
            if pid and name:names[pid].add(name.strip())
    return {pid:sorted(v) for pid,v in names.items()}


def name_status(names):
    if len(names)!=1:
        return "AMBIGUOUS_OR_UNRESOLVED_ACCEPTED_NAME"
    tokens=names[0].split()
    # Explicitly reject infraspecific/aggregate/fuzzy names. Resolution to
    # GBIF speciesKey must be separately verified in the next step.
    if len(tokens)!=2 or not re.fullmatch(r"[A-Z][A-Za-z-]+",tokens[0]) or not re.fullmatch(r"[a-z][a-z-]+",tokens[1]):
        return "NOT_EXACT_WCVP_BINOMIAL"
    return "EXACT_SOURCE_BINOMIAL_GBIF_RESOLUTION_REQUIRED"


def assemble(model,regions,crosswalk,positive_source,juvenile_result,names):
    if positive_source["status"]!="COMPLETE_48":
        raise ValueError("48 positive country source was not complete")
    observed=positive_source["observations"]
    pairs={(r["species"],r["country"]):r for r in observed}
    if len(pairs)!=48 or len({s for s,_ in pairs})!=20 or len({c for _,c in pairs})!=23:
        raise ValueError("Positive country source cohort identity changed")
    if sum(x["added_resource_country"] for x in observed)!=26:
        raise ValueError("Country added-host exposure definition changed")
    sources={(r["species"],r["country"]):r for r in juvenile_result["cells"]}
    if set(sources[r] for r in FIXED_CASES if r not in sources):
        raise ValueError("Fixed juvenile stage metadata no longer available")
    if not FIXED_CASES<=set(sources):
        raise ValueError("Cannot find all five source-locked stage country cases")
    if any(sources[k]["added_host_country"]!=1 for k in FIXED_CASES):
        raise ValueError("Stage case lost original added-host flag")
    index={r:i for i,r in enumerate(regions)}
    sp_index={s:i for i,s in enumerate(model["species"])}
    if len(sp_index)!=239:
        raise ValueError("Augmented HOSTS butterfly source size drifted")
    output=[]
    for (sp,ct),frozen in sorted(pairs.items()):
        j=sp_index[sp]
        material=crosswalk[ct]["material"]
        guard=crosswalk[ct]["guard"]
        if not material or not guard:
            raise ValueError("Missing source country material/guard regions")
        # The original 48 country positive cells were checked to be
        # outside documented native-host resource geography.
        native=model["native_species"][j]
        if any(native >> k & 1 for k in guard):
            raise ValueError(f"{sp}/{ct}: native host inside country guard; forbidden comparison")
        candidates=[]
        for host in model["species_hosts"][j]:
            introduced=model["introduced_sets"].get(host,set())
            overlap=material & introduced
            if overlap:
                candidate={"wcvp_accepted_id":host,
                    "botanical_names":names.get(host,[]),
                    "taxonomic_name_status":name_status(names.get(host,[])),
                    "country_material_introduced_wgsrpd3":sorted(regions[i] for i in overlap)}
                candidates.append(candidate)
        added=int(bool(candidates))
        if added!=frozen["added_resource_country"]:
            raise ValueError(f"Source candidate plant/exposure contradiction: {sp}/{ct} expected {frozen['added_resource_country']}, got {added}")
        key=(sp,ct)
        output.append({"species":sp,"country":ct,"added_host_country":added,
            "candidate_introduced_host_count":len(candidates),
            "country_material_regions":sorted(regions[i] for i in material),
            "country_guard_regions":sorted(regions[i] for i in guard),
            "has_previous_juvenile_stage_labels":key in sources,
            "selected_five_before_botanical_outcomes":key in FIXED_CASES,
            # Only five source-confirmed cases carry exact host IDs/names;
            # the other 43 receive counts, not an author-only relation table.
            "candidate_plants":candidates if key in FIXED_CASES else []})
    expected_flag_count=sum(x["added_host_country"] for x in output)
    selected=[x for x in output if x["selected_five_before_botanical_outcomes"]]
    if (len(output),expected_flag_count,len(selected))!=(48,26,5):
        raise ValueError("Fixed 48/26/5 host opportunity sample changed")
    return {"schema":SCHEMA,"executed_utc":datetime.now(timezone.utc).isoformat(),
       "status":"SOURCE_VERIFIED_CANDIDATE_PLANTS_NOT_LOCAL_HOST_USE",
       "source_eligible_butterflies":len(model["species"]),
       "country_positive_cells":48,"source_added_flag_cells":26,
       "focal_stage_case_count":5,
       "focal_unique_candidate_plant_ids":len({p["wcvp_accepted_id"] for x in selected for p in x["candidate_plants"]}),
       "focal_candidate_plant_links":sum(len(x["candidate_plants"]) for x in selected),
       "focal_exact_binomial_links":sum(p["taxonomic_name_status"]=="EXACT_SOURCE_BINOMIAL_GBIF_RESOLUTION_REQUIRED" for x in selected for p in x["candidate_plants"]),
       "source_host_opportunity_only_not_measured_local_plant_occurrence":True,
       "focal_five_cases":selected,
       "all_48_species_country_candidate_counts":[{k:v for k,v in x.items() if k!="candidate_plants"} for x in output]}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--crosswalk",type=Path,required=True)
    ap.add_argument("--positive48",type=Path,required=True)
    ap.add_argument("--juvenile-readout",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    from analyze_bce_augmented_plant_network_climate_null_20261019 import source_augmented_model
    model,regions,climate,provenance=source_augmented_model(a.source_root)
    links=load_country_links(a.crosswalk,regions)
    positive=json.loads(a.positive48.read_text(encoding="utf8"))
    juv=json.loads(a.juvenile_readout.read_text(encoding="utf8"))
    names=plant_names(a.source_root)
    output=assemble(model,regions,links,positive,juv,names)
    output["source_verification"]={
        "author_archive_manifest_sha256":provenance["source_provenance"]["source_frozen_manifest_sha256"] if isinstance(provenance.get("source_provenance"),dict) else provenance.get("source_frozen_manifest_sha256"),
        "crosswalk_sha256":sha(a.crosswalk),
        "positive48_sha256":sha(a.positive48),
        "juvenile_readout_sha256":sha(a.juvenile_readout),
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("HOST_CANDIDATES_FIVE="+json.dumps({k:v for k,v in output.items() if k not in ("all_48_species_country_candidate_counts",)},separators=(",",":")),flush=True)
    print("HOST_CANDIDATES_CHECKS="+json.dumps({k:v for k,v in output.items() if k not in ("focal_five_cases","all_48_species_country_candidate_counts")},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()

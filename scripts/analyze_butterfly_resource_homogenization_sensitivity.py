#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import median


def base_binomial(name: str) -> str:
    parts = str(name).strip().split()
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else ""


def genus(name: str) -> str:
    parts = str(name).strip().split()
    return parts[0] if parts else ""


def load_descriptors(path: Path):
    out = {}
    with path.open(newline="", encoding="utf-8") as h:
        for row in csv.DictReader(h):
            out[str(row["species"]).strip()] = {
                "host_family_count": float(row["host_family_count"]),
                "native_resource_units": int(row["host_wgsrpd3_unit_count"]),
            }
    return out


def load_pairs(path: Path):
    by_butterfly: dict[str, set[str]] = defaultdict(set)
    meta: dict[str, dict[str, object]] = {}
    with path.open(newline="", encoding="utf-8") as h:
        for row in csv.DictReader(h):
            sp = str(row["insect_species"]).strip()
            hid = str(row["accepted_plant_name_id"]).strip()
            if not sp or not hid:
                continue
            by_butterfly[sp].add(hid)
            accepted = str(row.get("accepted_name") or "").strip()
            input_name = str(row.get("input_host_name") or "").strip()
            fam = str(row.get("family") or "").strip()
            if hid not in meta:
                meta[hid] = {
                    "accepted_name": accepted,
                    "input_host_names": set(),
                    "family": fam,
                }
            if input_name:
                meta[hid]["input_host_names"].add(input_name)
    return by_butterfly, meta


def load_units(path: Path):
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as h:
        for row in csv.DictReader(h):
            hid = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if hid and unit:
                out[hid].add(unit)
    return out


def crop_exclusions(meta, crop, conservative: bool) -> set[str]:
    exact = set(map(str, crop["exact_binomials"]))
    wildcard_genera = set(map(str, crop["genus_wildcards"]))
    out = set()
    for hid, m in meta.items():
        names = [str(m["accepted_name"])] + sorted(m["input_host_names"])
        if any(base_binomial(n) in exact for n in names if n):
            out.add(hid)
            continue
        if conservative and any(genus(n) in wildcard_genera for n in names if n):
            out.add(hid)
    return out


def mean_region_jaccard(species, regions, units_by_species):
    idx = {sp:i for i,sp in enumerate(species)}
    masks = {r:0 for r in regions}
    for sp in species:
        bit = 1 << idx[sp]
        for r in units_by_species[sp]:
            if r in masks:
                masks[r] |= bit
    vals=[]
    rr=list(regions)
    for i in range(len(rr)):
        a=masks[rr[i]]
        for j in range(i+1,len(rr)):
            b=masks[rr[j]]
            u=(a|b).bit_count()
            if u:
                vals.append((a&b).bit_count()/u)
    return sum(vals)/len(vals), median(vals), len(vals)


def mean_species_jaccard(species, units_by_species):
    vals=[]
    for i in range(len(species)):
        a=units_by_species[species[i]]
        for j in range(i+1,len(species)):
            b=units_by_species[species[j]]
            u=len(a|b)
            if u:
                vals.append(len(a&b)/u)
    return sum(vals)/len(vals), median(vals), len(vals)


def reconstruct(species, pairs, native, contemporary, excluded):
    n_by={}
    c_by={}
    for sp in species:
        hosts=set(pairs.get(sp,set()))-excluded
        nu=set(); cu=set()
        for hid in hosts:
            nu |= set(native.get(hid,set()))
            cu |= set(contemporary.get(hid,set()))
        n_by[sp]=nu
        c_by[sp]=cu
    return n_by,c_by


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv",type=Path,required=True)
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--native-distribution-csv",type=Path,required=True)
    ap.add_argument("--contemporary-distribution-csv",type=Path,required=True)
    ap.add_argument("--crop-json",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    d=load_descriptors(args.descriptors_csv)
    pairs,meta=load_pairs(args.insect_host_csv)
    native=load_units(args.native_distribution_csv)
    contemporary=load_units(args.contemporary_distribution_csv)
    crop=json.loads(args.crop_json.read_text(encoding="utf-8"))

    species=[sp for sp,x in d.items() if x["host_family_count"]>0 and x["native_resource_units"]>0]
    if len(species)!=239:
        raise RuntimeError(f"expected 239 species, got {len(species)}")

    baseline_n,baseline_c=reconstruct(species,pairs,native,contemporary,set())
    for sp in species:
        if len(baseline_n[sp])!=d[sp]["native_resource_units"]:
            raise RuntimeError(f"baseline native drift for {sp}")
    baseline_regions=sorted(set().union(*baseline_n.values()))
    if len(baseline_regions)!=355:
        raise RuntimeError(f"expected 355 baseline native-active regions, got {len(baseline_regions)}")

    variants={
        "baseline":set(),
        "fao_exact_binomial":crop_exclusions(meta,crop,False),
        "fao_exact_plus_spp_genus":crop_exclusions(meta,crop,True),
        "exclude_all_poaceae_hosts":{hid for hid,m in meta.items() if m["family"]=="Poaceae"},
    }
    out={}
    for label,excluded in variants.items():
        n,c=reconstruct(species,pairs,native,contemporary,excluded)
        added={sp:c[sp]-n[sp] for sp in species}
        total_native=sum(len(n[sp]) for sp in species)
        total_cont=sum(len(c[sp]) for sp in species)
        total_added=sum(len(added[sp]) for sp in species)
        expanded=sum(bool(added[sp]) for sp in species)
        region_n=mean_region_jaccard(species,baseline_regions,n)
        region_c=mean_region_jaccard(species,baseline_regions,c)
        sp_n=mean_species_jaccard(species,n)
        sp_c=mean_species_jaccard(species,c)
        out[label]={
            "excluded_host_species":len(excluded),
            "species_with_any_remaining_native_resource":sum(bool(n[sp]) for sp in species),
            "expanded_species":expanded,
            "native_butterfly_x_region_units":total_native,
            "contemporary_butterfly_x_region_units":total_cont,
            "added_butterfly_x_region_units":total_added,
            "aggregate_proportional_increase":total_added/total_native if total_native else None,
            "regional_resource_assemblage":{
                "mean_jaccard_native":region_n[0],
                "mean_jaccard_contemporary":region_c[0],
                "delta":region_c[0]-region_n[0],
                "relative_change":region_c[0]/region_n[0]-1 if region_n[0] else None,
                "pairs":region_n[2],
            },
            "butterfly_resource_geography_overlap":{
                "mean_jaccard_native":sp_n[0],
                "mean_jaccard_contemporary":sp_c[0],
                "delta":sp_c[0]-sp_n[0],
                "relative_change":sp_c[0]/sp_n[0]-1 if sp_n[0] else None,
                "pairs":sp_n[2],
            }
        }

    if out["baseline"]["added_butterfly_x_region_units"]!=14553:
        raise RuntimeError("baseline added-unit drift")
    payload={
        "schema":"chocho_butterfly_resource_homogenization_sensitivity_v0.1",
        "status":"SUCCESS_CROP_AND_POACEAE_HOMOGENIZATION_SENSITIVITY",
        "comparison_regions":"The same 355 WGSRPD3 regions with at least one native resource opportunity in the baseline reconstruction are used for every regional-composition variant.",
        "variants":out,
        "claim_boundary":"These are post-hoc sensitivity reconstructions. Crop and Poaceae exclusions change the host set and therefore the native as well as introduced resource matrix; they are not causal estimates of agriculture or grasses."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))


if __name__=="__main__":
    main()

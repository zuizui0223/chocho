#!/usr/bin/env python3
"""Source-qualified alien-established butterfly × recipient host opportunities.

Zenodo 2024 original 18 butterfly species and 35 alien-positive recipient
records, against the exact original 239 butterfly x accepted-host WCVP
range source. This is NOT a transported-insect/feeding/establishment trial.
Broad zones or Okinawa→Nansei botanical proxy do not enter strict main.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from analyze_hosts_location_wcvp_origin_20261027 import (
    accepted_pairs, read_rows, sha256, wcvp_range, load_country_atoms,
)

SCHEMA="chocho_original239_zenodo35_recipient_host_opportunity_20261027_v01"
ZENODO_SHA="05483a2a9a7c926fb5105c44058ee52a4598039dc82c5bb10d1db284b3b66fe7"
ORIGINAL_ZIP_SHA="408cb206d7ce4056cee1a6efb3db4650ff7c2663a7d6d6003e323"
EXACT_ZONES={
    "Hawaii":("HAW","Hawaii"),
    "Ogasawara":("OGA","Ogasawara-shoto"),
    "Galapagos":("GAL","Galápagos"),
    "Korea":("KOR","Korea"),
}
COUNTRY_ZONES={"Japan":"JP","New Zealand":"NZ","South Africa":"ZA","Chile":"CL"}
BOTANICAL_PROXY={"Okinawa":("NNS","Nansei-shoto")}
BROAD_UNRESOLVED={"Europe","North America"}
ALL_ZONE_NAMES=set(EXACT_ZONES)|set(COUNTRY_ZONES)|set(BOTANICAL_PROXY)|BROAD_UNRESOLVED
AREA_HOST_FRACTION=.01
SOURCE_COVERAGE_THRESHOLDS=(.95,.90)
RICHNESS_BINS=("1-5","6-20","21-50","51+")


def host_bin(n):
    if n<1:
        raise ValueError("Globally known accepted host count must be positive")
    if n<=5:return "1-5"
    if n<=20:return "6-20"
    if n<=50:return "21-50"
    return "51+"


def zone_mapping(zone,regions,atoms,meta,*,coverage=.95):
    """Fail closed for country source coverage or coarser biological ranges."""
    if zone in BROAD_UNRESOLVED:
        return {"source_zone":zone,"mapping_status":"UNSCORABLE_SUPRA_COUNTRY_REGION"}
    if zone in EXACT_ZONES or zone in BOTANICAL_PROXY:
        proxy=zone in BOTANICAL_PROXY
        code,name=(BOTANICAL_PROXY if proxy else EXACT_ZONES)[zone]
        if regions.get(code)!=name:
            raise ValueError("Original Level3 botany name/code changed: "+zone)
        return {
            "source_zone":zone,
            "mapping_status":(
                "BOTANICAL_NANSEI_PROXY_NOT_EXACT_OKINAWA" if proxy else
                "EXACT_WCVP_LEVEL3_SOURCE_REGION"),
            "source_region_code":code,
            "source_region_name":name,
            "primary_geographic_evidence":not proxy,
        }
    if zone not in COUNTRY_ZONES:
        raise ValueError("Unknown independent 2024 alien receiving zone: "+zone)
    iso=COUNTRY_ZONES[zone]
    m=meta.get(iso)
    if not m or iso not in atoms or not m.get("is_geo_safe",False):
        return {"source_zone":zone,
                "mapping_status":"COUNTRY_ISO_OR_SOURCE_GEO_UNSAFE",
                "original_country_iso2":iso}
    score=float(m["guard_coverage"])
    if score+1e-12 < coverage:
        return {"source_zone":zone,
                "mapping_status":"ORIGINAL_COUNTRY_SOURCE_AREA_COVERAGE_BELOW_GATE",
                "original_country_iso2":iso,
                "source_original_mapped_area_fraction":score}
    return {
        "source_zone":zone,
        "mapping_status":"ORIGINAL_ISO_COUNTRY_ATOM_AREA_RECONSTRUCTION",
        "original_country_iso2":iso,
        "source_original_mapped_area_fraction":score,
        "mapped_source_botanical_atom_count":len(atoms[iso]),
    }


def country_plant(pid,iso,atoms,meta,native,current,*,threshold=AREA_HOST_FRACTION):
    """Multiple disjoint atoms may have the same plant native and introduced."""
    area=float(meta[iso]["country_area_km2"])
    if area<=0:
        raise ValueError("Invalid source original country land area")
    intro=native_count=anywild=0.0
    for regions,area_atom in atoms[iso]:
        # Status can be mixed if overlapping L3 botanical geometry or
        # plant native in one L3 and introduced in another.
        native_in=any(reg in native.get(pid,set()) for reg in regions)
        intro_in=any(reg in current.get(pid,set()) and
                     reg not in native.get(pid,set()) for reg in regions)
        if native_in:native_count+=area_atom
        if intro_in:intro+=area_atom
        if native_in or intro_in:anywild+=area_atom
    return native_count/area>=threshold-1e-12, intro/area>=threshold-1e-12,anywild/area


def classify_plant(pid,mapping,atoms,meta,native,current):
    code=mapping.get("source_region_code")
    if code is not None:
        n=code in native.get(pid,set())
        i=code in current.get(pid,set()) and not n
    elif mapping.get("original_country_iso2"):
        n,i,_=country_plant(pid,mapping["original_country_iso2"],
                            atoms,meta,native,current)
    else:
        raise ValueError("Attempt to use botanically unresolved source zone as absence")
    if n and i:return "BOTH_NATIVE_AND_INTRODUCED_WILD"
    if n:return "NATIVE_WILD"
    if i:return "INTRODUCED_WILD"
    return "NO_WCVP_WILD_RECORDED_OR_BELOW_AREA_GATE"


def species_source_record(sp,zone,hostIDs,mapping,atoms,meta,native,current,source_records):
    statuses=Counter(classify_plant(pid,mapping,atoms,meta,native,current)
                     for pid in hostIDs)
    n=statuses["NATIVE_WILD"]+statuses["BOTH_NATIVE_AND_INTRODUCED_WILD"]
    i=statuses["INTRODUCED_WILD"]+statuses["BOTH_NATIVE_AND_INTRODUCED_WILD"]
    if n and i: status="BOTH_NATIVE_AND_INTRODUCED_POTENTIAL_HOSTS"
    elif i: status="ONLY_INTRODUCED_POTENTIAL_HOSTS"
    elif n: status="ONLY_NATIVE_POTENTIAL_HOSTS"
    else: status="NO_WCVP_WILD_KNOWN_HOST_RECORDED"
    if sum(statuses.values())!=len(hostIDs):
        raise ValueError("Source botanical host species accounting failed")
    return {
        "species":sp,
        "independent_2024_recipient_zone":zone,
        "source_geographic_mapping_status":mapping["mapping_status"],
        "accepted_original_known_global_larval_hosts":len(hostIDs),
        "native_only_wild_source_hosts":statuses["NATIVE_WILD"],
        "introduced_only_wild_source_hosts":statuses["INTRODUCED_WILD"],
        "source_plant_native_and_introduced_in_country_atoms":statuses[
            "BOTH_NATIVE_AND_INTRODUCED_WILD"],
        "unlisted_or_below_area_gate_known_hosts":statuses[
            "NO_WCVP_WILD_RECORDED_OR_BELOW_AREA_GATE"],
        "source_any_native_WCVP_wild_host":bool(n),
        "source_any_introduced_WCVP_wild_host":bool(i),
        "source_plant_host_support_class":status,
        "original_positive_established_source_region_records":[
            {
                "historical_year":r.get("year_source_first_report"),
                "intentional_release":r.get("source_intentional_release_flag"),
                "indoors_outdoors":r.get("source_established_indoor_outdoor_flag"),
                "eradicated":r.get("source_eradication_flag"),
            } for r in source_records if r["recipient_zone"]==zone
        ],
        "source_global_HOSTS_link_not_direct_local_larval_feeding":True,
    }


def hypergeom(n,k,m):
    if not (n>0 and 0<=k<=n and 0<=m<=n):
        raise ValueError("Impossible richness-binned positive source population")
    d=math.comb(n,m)
    return {v:math.comb(k,v)*math.comb(n-k,m-v)/d
            for v in range(max(0,m+k-n),min(k,m)+1)}


def artificial_reference(background,positive):
    by_bin=defaultdict(list)
    marked=Counter(host_bin(r["accepted_original_known_global_larval_hosts"])
                   for r in positive)
    for row in background:
        by_bin[host_bin(row["accepted_original_known_global_larval_hosts"])].append(row)
    pmf={0:1.0}
    mean=0.0
    bins=[]
    for b in RICHNESS_BINS:
        rows=by_bin[b]
        n=len(rows)
        k=sum(r["source_any_introduced_WCVP_wild_host"] for r in rows)
        m=marked[b]
        if m and n==0:
            raise ValueError("Original source host richness stratum unrepresented")
        bins.append({"source_original_host_breadth_bin":b,"source_all239":n,
                     "source_with_intro_wild_host":k,
                     "independently_established_positive_species":m})
        if not n: continue
        mean+=m*k/n
        if m==0:continue
        tmp=defaultdict(float)
        for i,pi in pmf.items():
            for j,pj in hypergeom(n,k,m).items():
                tmp[i+j]+=pi*pj
        pmf=dict(tmp)
    if not math.isclose(sum(pmf.values()),1,abs_tol=1e-8):
        raise ValueError("Convolved artificial source reference probabilities drifted")
    obs=sum(v["source_any_introduced_WCVP_wild_host"] for v in positive)
    ks=sorted(pmf)
    def q(frac):
        total=0
        for v in ks:
            total+=pmf[v]
            if total>=frac-1e-12:return v
        return ks[-1]
    return {
        "source_established_butterflies_with_intro_wild_hosts_observed":obs,
        "source_accepted_breadth_matched_artificial_expected":float(mean),
        "source_reference_95pct_interval":[q(.025),q(.975)],
        "source_tail_NOT_an_ecological_p":sum(v for k,v in pmf.items() if k>=obs),
        "source_host_breadth_bin_margins":bins,
        "NO_butterfly_arrival_or_nonestablishment_data_for_this_reference":True,
    }


def run(args):
    obj=json.loads(args.zenodo_2024_positive_json.read_text())
    if obj.get("schema")!="chocho_239_butterfly_original_Zenodo_positive_alien_establishment_20261027_v01":
        raise ValueError("Wrong original 2024 independent Zenodo alien positive source")
    if obj["source_zenodo"]["record_id"]!=11519262 or \
       obj["source_zenodo"]["sources"]["occurrence_table_20230317.csv"]["sha256"]!=\
       "05483a2a9a7c926fb5105c44058ee52a4598039dc82c5bb10d1db284b3b66fe7":
        raise ValueError("Original Zenodo 2024 fixed source taxonomy SHA changed")
    source=obj["source_positive_matches"]
    original18=source["source_positive_original239_species"]
    positives=[(item,zone) for item in original18
               for zone in item["source_documented_non_native_recipient_areas"]]
    if len(original18)!=18 or len(positives)!=35 or \
       set(zone for _,zone in positives)!=ALL_ZONE_NAMES:
        raise ValueError("Original 18 butterflies × 35 established recipient records changed")
    root=args.author_source
    manifest=json.loads((root/"SOURCE_MANIFEST.json").read_text())
    if manifest.get("inputs_count")!=11:
        raise ValueError("Not the original signed author-only accepted-host archive")
    for path,digest in manifest["sha256"].items():
        if sha256(root/path)!=digest:
            raise ValueError("Original author-only WCVP input SHA256 changed "+path)
    original_list=read_rows(root/"inputs/anthropogenic_species_metrics.csv")
    focal={r["species"] for r in original_list}
    if len(original_list)!=239 or len(focal)!=239:
        raise ValueError("Incorrect original 239 host-resource taxa")
    accepted,qc,_=accepted_pairs(
        root/"inputs/original/insect_host_accepted.csv",focal)
    if qc["source_ambiguous_accepted_id_names"]:
        raise ValueError("Ambiguous original accepted HOSTS taxon")
    hosts=defaultdict(set)
    for (sp,name),pid in accepted.items():
        hosts[sp].add(pid)
    if sum(map(len,hosts.values()))!=2596:
        raise ValueError("Wrong original 2596 accepted butterfly host pairs")
    native,current=wcvp_range(
        root/"inputs/original/native_extant_nondoubtful_wgsrpd3.csv",
        root/"inputs/original/extant_nondoubtful_wgsrpd3.csv")
    levels=json.loads((root/"inputs/level3.geojson").read_text())["features"]
    regions={f["properties"]["LEVEL3_COD"]:f["properties"]["LEVEL3_NAM"]
             for f in levels}
    if len(regions)!=len(levels) or len(levels)!=369:
        raise ValueError("Original WCVP Level3 region map changed")
    atoms,meta,atom_regions,geography=load_country_atoms(
        args.country_atoms,args.country_atoms_manifest,args.country_summary)
    if any(code not in regions for code in atom_regions):
        raise ValueError("Source country disjoint atom has botanical code not in source geography")
    result={}
    for cutoff in SOURCE_COVERAGE_THRESHOLDS:
        zone_mapping_by_name={
            zone:zone_mapping(zone,regions,atoms,meta,coverage=cutoff)
            for zone in ALL_ZONE_NAMES
        }
        source_outputs=[]
        source_controls=[]
        by_region=defaultdict(list)
        for source_positive,zone in positives:
            mapping=zone_mapping_by_name[zone]
            if mapping["mapping_status"] in (
                "EXACT_WCVP_LEVEL3_SOURCE_REGION",
                "ORIGINAL_ISO_COUNTRY_ATOM_AREA_RECONSTRUCTION"):
                status="MAIN_COMPARABLE"
            elif mapping["mapping_status"]=="BOTANICAL_NANSEI_PROXY_NOT_EXACT_OKINAWA":
                status="SENSITIVITY_BOTANICAL_PROXY"
            else:
                status="UNSCORABLE_ORIGINAL_SOURCE_GEOGRAPHY"
            if status=="UNSCORABLE_ORIGINAL_SOURCE_GEOGRAPHY":
                source_outputs.append({
                    "species":source_positive["species"],
                    "recipient_zone":zone,
                    "geographic_status":status,
                    "source_mapping_status":mapping["mapping_status"],
                    "source_gap_is_not_absence_of_larval_hosts":True,
                })
                continue
            sp=source_positive["species"]
            entry=species_source_record(
                sp,zone,hosts[sp],mapping,atoms,meta,native,current,
                source_positive["original_source_region_and_year_records"])
            entry["geographic_status"]=status
            source_outputs.append(entry)
            by_region[zone].append(entry)
        for zone in sorted(by_region):
            mapping=zone_mapping_by_name[zone]
            background=[
                species_source_record(sp,zone,hosts[sp],mapping,atoms,meta,native,current,[])
                for sp in sorted(focal)
            ]
            main=by_region[zone]
            cats=Counter(r["source_plant_host_support_class"] for r in main)
            source_controls.append({
                "zone":zone,
                "mapping_scope":mapping["mapping_status"],
                "source_original_geo_coverage":mapping.get("source_original_mapped_area_fraction"),
                "positive_species":len(main),
                "positive_species_support":dict(cats),
                "source_background_all_original239_species_with_introduced_hosts":sum(
                    r["source_any_introduced_WCVP_wild_host"] for r in background),
                "source_background_all_original239_species_with_native_hosts":sum(
                    r["source_any_native_WCVP_wild_host"] for r in background),
                "source_239_taxon_BREADTH_matched_source_reference":
                    artificial_reference(background,main),
                "reference_NOT_butterfly_arrival_establishment_success":True,
            })
        if len(source_outputs)!=35:
            raise ValueError("35 source positive species×zone identities not conserved")
        counts=Counter(v["geographic_status"] for v in source_outputs)
        support=Counter(v["source_plant_host_support_class"] for v in source_outputs
                        if v["geographic_status"]=="MAIN_COMPARABLE")
        result[str(cutoff)]={
            "source_positive35_mapped_geography_counts":dict(counts),
            "source_positive35_principal_recipient_known_host_support":dict(support),
            "positive_original35_each_with_scope_and_host_counts":source_outputs,
            "per_recipient_positive_and_original239_background":source_controls,
            "all_Zenodo_receiver_zone_to_exact_WCVP_source_mapping":zone_mapping_by_name,
        }
    primary=result[str(.95)]
    hawaii=[r for r in primary["positive_original35_each_with_scope_and_host_counts"]
            if r["geographic_status"]=="MAIN_COMPARABLE" and
            r["independent_2024_recipient_zone"]=="Hawaii"]
    if len(hawaii)!=10 or sum(x["introduced_only_wild_source_hosts"] for x in hawaii)!=82 or \
       sum(x["native_only_wild_source_hosts"] for x in hawaii)!=5 or \
       sum(x["unlisted_or_below_area_gate_known_hosts"] for x in hawaii)!=251:
        raise ValueError("Independent earlier ten Hawaii WCVP original accepted host result drift")
    if primary["source_positive35_mapped_geography_counts"]!={
        "MAIN_COMPARABLE":24,
        "SENSITIVITY_BOTANICAL_PROXY":3,
        "UNSCORABLE_ORIGINAL_SOURCE_GEOGRAPHY":8,
    }:
        raise ValueError("Original source 35 geographic mapping categories changed")
    if result[str(.90)]["source_positive35_mapped_geography_counts"]!={
        "MAIN_COMPARABLE":25,
        "SENSITIVITY_BOTANICAL_PROXY":3,
        "UNSCORABLE_ORIGINAL_SOURCE_GEOGRAPHY":7,
    }:
        raise ValueError("Original source Chile <95 country coverage sensitivity changed")
    output={
        "schema":SCHEMA,
        "status":"SOURCE_35_ALIEN_ESTABLISHMENT_POSITIVES_RECIPIENT_WCVP_HOST_GEOGRAPHY_NO_ARRIVAL_ESTIMATION",
        "source":{
            "original_2024_zenodo_positive_json_sha256":sha256(args.zenodo_2024_positive_json),
            "original_2024_matched_aliens":18,
            "original_2024_positive_species_recipient_records":35,
            "original_239_butterflies":239,
            "original_accepted_butterfly_host_pairs":2596,
            "original_author_zip_SHA256":ORIGINAL_ZIP_SHA,
            "original_county_atoms":geography,
            "original_WCVP_level3_regions":369,
            "country_species_host_min_land_fraction":AREA_HOST_FRACTION,
        },
        "geographic_source_guard95_and_sensitivity90":result,
        "scientific_limits":[
            "These are independent established alien butterfly-positive source records, not all transported taxa or arrival attempts.",
            "No candidate without 2024 positive source is treated as a failed/nonestablished butterfly.",
            "Original global HOSTS documented larval-host edges are potential source food only; no local wild larval feeding confirmation.",
            "Original WCVP wild accepted plant distributions omit cultivated-only nurseries, gardens, food crops and transient botanical plants.",
            "Broad Europe/North America and botanically broader Okinawa→Nansei regions cannot be treated as exact plant locality tests.",
            "Countries with source geographic coverage below a predeclared cutoff remain unscored, not evidence of host absence.",
            "For country-wide source positive entries the actual butterfly establishment point may not overlap the mapped source host polygon.",
            "A source butterfly may be listed as historic/eradicated/indoors or with a suspect first year; preserve flag, no current persistence asserted.",
            "Per-region hypergeometric references are artificial taxonomic documentation exchangeability, NOT ecological p values or invasion probabilities.",
            "Recip zone names such as Japan and Okinawa/Ogasawara overlap, individual butterflies recur across regions, no cross-row independence.",
        ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "status":output["status"],
        "primary35_mapping":primary["source_positive35_mapped_geography_counts"],
        "primary_24_host_status":primary["source_positive35_principal_recipient_known_host_support"],
        "receiving_regions":[
            {k:r[k] for k in ("zone","positive_species",
                              "source_background_all_original239_species_with_introduced_hosts",
                              "source_background_all_original239_species_with_native_hosts",
                              "positive_species_support",
                              "source_239_taxon_BREADTH_matched_source_reference")}
            for r in primary["per_recipient_positive_and_original239_background"]
        ],
        "nonintroduction_WCVP_support_exceptions":[
            r for r in primary["positive_original35_each_with_scope_and_host_counts"]
            if r.get("geographic_status")=="MAIN_COMPARABLE" and
               not r.get("source_any_introduced_WCVP_wild_host")
        ],
        "area90_mapping":result[str(.90)]["source_positive35_mapped_geography_counts"]
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--zenodo-2024-positive-json",type=Path,required=True)
    p.add_argument("--country-atoms",type=Path,required=True)
    p.add_argument("--country-atoms-manifest",type=Path,required=True)
    p.add_argument("--country-summary",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args)

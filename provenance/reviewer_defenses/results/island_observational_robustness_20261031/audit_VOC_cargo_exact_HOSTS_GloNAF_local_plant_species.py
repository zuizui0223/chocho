#!/usr/bin/env python3
"""Historical VOC commodity genus -> EXACT modern island plant + larval host.

The historical cargo commodity is known only at genus/category, not the
transported biological plant species. This is a strict contemporary
*catalogue triad*: a butterfly independently indexed in that island in
2011–2017, the exact HOSTS larval plant binomial recorded for that butterfly,
and the same plant binomial reported as a naturalized alien in the island's
GloNAF region. The genus of that plant must also match a predeclared
historical cargo commodity category on an original 1595–1698 VOC route.

This does NOT identify live historic plant cargo, colonial introduction
dates, port-local plant distribution or local butterfly larval feeding.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys
from zipfile import ZipFile

BASE=Path(__file__).resolve().parent
CONTINUITY=BASE.parent/"island_host_continuity_20261030/island_host_continuity_events.csv"
ISLAND_REGION=BASE/"island_port_colonial_10km.islands.csv"
CARGO=BASE/"voc_17c_cargo_on_verified_island_routes_v1.json"
MAPLOCK=BASE/"voc_17c_commodity_to_HOSTS_genus_lock_20261101.json"
LOCK=BASE/"voc_cargo_host_species_naturalized_islands_lock_20261102.json"
PANEL=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
HOSTS=Path("_external/HOSTS/resource.csv")
GLO_NAF=Path("_temp_historical_sources/lenzner_6640264.zip")
TDWG=Path("_external/tdwg/geojson/level3.geojson")
OUT=BASE/"voc_cargo_exact_HOSTS_GloNAF_local_plant_species_v1.json"
PAIRS=BASE/"voc_cargo_exact_HOSTS_GloNAF_local_plant_species_pairs_v1.csv"
SIX={"12":"SUM","17":"SUL","19":"JAW","31":"SRL","43":"JAP","44":"TAI"}
CUTOFFS=(.01,.05,.1)
HISTORICAL_MD5="16f5585d2500676951a4b664bbf0ffde"


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def normalized_host_binomial(genus,epithet):
    g=str(genus or "").strip()
    e=str(epithet or "").strip()
    if not g or not e or "\ufffd" in (g+e):
        return None
    s=e if e.casefold().startswith(g.casefold()+" ") else g+" "+e
    if len(s)<5 or " sp." in s.casefold():
        return None
    return s.casefold()


def record_first_butterflies(path=CONTINUITY):
    source=read_rows(path)
    if len(source)!=5665:
        raise RuntimeError("Original GBIF 2011–17 island species source changed")
    events=defaultdict(set)
    for r in source:
        yr=int(r["first_record_year"])
        if not (2011<=yr<=2017):
            raise RuntimeError("Archived source butterfly first-record years outside window")
        island=str(r["landmass_id"])
        if island in SIX:
            butterfly=str(r["species"]).strip().casefold()
            if butterfly in events[island]:
                raise RuntimeError("Duplicate original butterfly first event in "+island)
            events[island].add(butterfly)
    if sum(len(events[k]) for k in SIX)!=49:
        raise RuntimeError("Original six-island first indexed species count changed")
    return dict(events)


def load_sources():
    with LOCK.open(encoding="utf8") as f:
        lock=json.load(f)
    if lock["schema"]!="chocho_VOC_cargo_exact_HOSTS_host_species_GloNAF_naturalized_region_v1":
        raise RuntimeError("Exact species triad source lock missing")
    root=read_rows(ISLAND_REGION)
    region_of={r["landmass_id"]:r for r in root}
    for island,region in SIX.items():
        r=region_of.get(island)
        if r is None or r["region_codes_mapped"]!=region or r["eligible_for_history_test"]!="1":
            raise RuntimeError("Original 51-island WPI/TDWG3 source matched regions changed")
    with CARGO.open(encoding="utf8") as f:
        cargo=json.load(f)
    if cargo["verified_historic_island_edges"]!=9 or cargo["historical_source_original_cargo_rows"]!=7737:
        raise RuntimeError("Actual 9 original historic VOC island route source changed")
    pairs=cargo["pair_specific_aggregates"]
    if len(pairs)!=9 or sum(r["source_cargo_rows_joined_by_unrepeated_REISN"]>0 for r in pairs)!=7:
        raise RuntimeError("Original 7 cargo-supported routes changed")
    with MAPLOCK.open(encoding="utf8") as f:
        name_map=json.load(f)
    if name_map.get("schema")!="chocho_17c_VOC_commodity_name_to_focal239_HOSTS_genus_screen_v1":
        raise RuntimeError("Original commodity genus lock drift")
    category_genera={
        k:{v.casefold() for v in values}
        for k,values in name_map["genus_mapping"].items()
    }
    butterflies=record_first_butterflies()
    return pairs,category_genera,butterflies


def panel_host_species(hosts=HOSTS,panel=PANEL):
    focal={r["species"].casefold() for r in read_rows(panel)}
    if len(focal)!=239:
        raise RuntimeError("Original 239 butterfly resource panel changed")
    native=defaultdict(set)
    raw=read_rows(hosts)
    for r in raw:
        b=(r["Insect Genus"].strip()+" "+r["Insect Species"].strip()).casefold()
        if b not in focal:
            continue
        plant=normalized_host_binomial(r["Hostplant Genus"],r["Hostplant Species"])
        if plant:
            native[b].add(plant)
    if len(raw)!=140485 or len(native)!=239:
        raise RuntimeError("Source HOSTS 140485 records / 239 butterfly taxonomy mismatch")
    return dict(native),{"original_HOSTS_rows":len(raw),
                         "focal_butterfly_species":len(focal),
                         "focal_butterflies_with_any_binomial":len(native)}


def source_region_masks(glonaf=GLO_NAF,tdwg=TDWG,thresholds=CUTOFFS):
    """Rebuild ONLY original six TDWG regions from source GloNAF polygons.

    Numeric bit-position genome is the exact same original GloNAF source
    string order used by the public matched-284-region analysis. Never
    synthesize species records for an overlapping region with no source.
    """
    from shapely.geometry import shape
    from shapely.ops import transform as geo_transform
    from shapely.strtree import STRtree
    from shapely.validation import make_valid
    from shapely.errors import GEOSException
    from _colonial_source_imports import load_glonaf_helpers
    helpers=load_glonaf_helpers()
    if helpers.check_source(glonaf) is not None:
        pass
    source_bits,source_audit=helpers.read_glonaf_taxa(glonaf)
    _,geom,project,invalid=helpers.read_geometries(glonaf)
    # read_geometries returns projected polygons and original source
    # region IDs, not botanical taxon bit positions.
    original_ids,geometry,project,invalid=helpers.read_geometries(glonaf)
    idx=STRtree(geometry)
    archived=json.loads(tdwg.read_text(encoding="utf8"))
    target={}
    for feature in archived["features"]:
        code=str(feature["properties"].get("LEVEL3_COD") or "")
        if code not in SIX.values():
            continue
        pg=shape(feature["geometry"])
        pg=make_valid(pg) if not pg.is_valid else pg
        pg=geo_transform(project.transform,pg)
        if not pg.is_valid:
            pg=make_valid(pg)
        if pg.is_empty or pg.area<=0:
            raise RuntimeError("Invalid mapped botanical region "+code)
        target[code]=pg
    if len(target)!=6:
        raise RuntimeError("Not all six predeclared TDWG region codes present")
    covered={q:{} for q in thresholds}
    nsrc={q:{} for q in thresholds}
    for code,g in target.items():
        overlaps=[]
        for k in idx.query(g):
            k=int(k)
            try:
                intersect=g.intersection(geometry[k])
            except GEOSException:
                intersect=make_valid(g).intersection(make_valid(geometry[k]))
            if not intersect.is_empty and intersect.area>0:
                overlaps.append((k,float(intersect.area)))
        for q in thresholds:
            bits=0
            included=0
            for k,area in overlaps:
                if helpers.classify_source_contribution(g.area,geometry[k].area,area,q):
                    bits|=source_bits.get(original_ids[k],0)
                    included+=1
            covered[q][code]=bits
            nsrc[q][code]=included
    return covered,nsrc,source_audit,invalid


def source_canonical_glonaf_aliases(glonaf=GLO_NAF):
    from _colonial_source_imports import load_glonaf_index
    aliases,audit=load_glonaf_index(glonaf)
    return aliases,audit


def local_triads(community,hosts,categories,aliases,region_bits):
    """All butterfly × exact HOSTS plant binomial × category matches.

    Return ONLY source supported *modern naturalized* records; genus-only
    matches do not qualify. The commodity category must have been actually
    recorded on the fixed historical route.
    """
    out=set()
    for butterfly in community:
        for plant in hosts.get(butterfly,set()):
            genus=plant.split(" ",1)[0]
            source_idxs=aliases.get(plant,())
            if not source_idxs:
                continue
            if not any(region_bits & (1<<i) for i in source_idxs):
                continue
            for category in categories:
                if genus in categories[category]:
                    out.add((butterfly,plant,category))
    return out


def aggregate_pair(a,b,cats,first,hosts,category_genera,aliases,masked):
    left=local_triads(first[a],hosts,
                      {k:category_genera[k] for k in cats},
                      aliases,masked[SIX[a]])
    right=local_triads(first[b],hosts,
                      {k:category_genera[k] for k in cats},
                      aliases,masked[SIX[b]])
    same=left&right
    return {
        "n_left":len(left),"n_right":len(right),"n_both_exact_triplets":len(same),
        "left":sorted(left),"right":sorted(right),"identical_butterfly_host_species_category_both":sorted(same),
        "n_left_unique_plants":len({p for _,p,_ in left}),
        "n_right_unique_plants":len({p for _,p,_ in right}),
        "n_both_unique_plants":len({p for _,p,_ in same}),
    }


def analyze():
    pairs,categories,events=load_sources()
    hosts,host_audit=panel_host_species()
    # Import crosswalk from pinned original source branch without copying
    # its data files, so the same GloNAF species-bit indexing is used.
    aliases,alias_audit=source_canonical_glonaf_aliases()
    masked,source_regions,glonaf_audit,invalid=source_region_masks()
    if alias_audit["distinct_glonaf_names"]!=glonaf_audit["distinct_species"]:
        raise RuntimeError("Original GloNAF bit-position source mismatch")
    comparisons={}
    for q in CUTOFFS:
        entries=[]
        for r in pairs:
            a,b=str(r["island_a"]),str(r["island_b"])
            cargo_rows=int(r["source_cargo_rows_joined_by_unrepeated_REISN"])
            if a not in SIX or b not in SIX:
                raise RuntimeError("Original historical pair island ID changed")
            cc={k for k,num in r["source_cargo_category_counts"].items()
                if num>0 and categories.get(k)}
            matched=aggregate_pair(a,b,cc,events,hosts,categories,aliases,masked[q])
            entries.append({
                "island_a":a,"island_b":b,
                "historical_route_1595_1698_leg_rows":r["observed_original_voyage_leg_rows"],
                "unambiguously_linked_original_cargo_rows":cargo_rows,
                "historically_recorded_commodity_categories":";".join(sorted(cc)),
                "eligible_historical_cargo_data":int(cargo_rows>0),
                "newly_indexed_butterfly_species_a":len(events[a]),
                "newly_indexed_butterfly_species_b":len(events[b]),
                "same_butterfly_exact_naturalized_HOSTS_species_both":int(
                    len(matched["identical_butterfly_host_species_category_both"])>0),
                "any_island_exact_naturalized_HOSTS_species":int(bool(matched["left"] or matched["right"])),
                "each_island_exact_naturalized_HOSTS_species":int(bool(matched["left"]) and bool(matched["right"])),
                "exact_source_taxon_category_match_a":matched["n_left"],
                "exact_source_taxon_category_match_b":matched["n_right"],
                "exact_source_taxon_category_match_shared_butterfly_and_plant":matched["n_both_exact_triplets"],
                "exact_source_triads_a":";".join("|".join(x) for x in matched["left"]),
                "exact_source_triads_b":";".join("|".join(x) for x in matched["right"]),
                "exact_source_triads_both":";".join("|".join(x) for x in matched["identical_butterfly_host_species_category_both"]),
                "source_plant_observations":("MODERN_REGION_NOT_PORT_SPECIFIC"
                    if (a=="43" or b=="43") else "MODERN_ISLAND_REGION_APPROXIMATION"),
            })
        relevant=[x for x in entries if x["eligible_historical_cargo_data"]]
        comparisons[str(q)]={
            "threshold":q,
            "source_target_WGSRPD3_region_GLoNAF_polygon_counts":source_regions[q],
            "original_historic_positive_island_pairs":len(entries),
            "history_routes_with_actual_cargo_ID_rows":len(relevant),
            "historical_cargo_links_with_exact_species_naturalized_host_match_one_endpoint":sum(
                x["any_island_exact_naturalized_HOSTS_species"] for x in relevant),
            "historical_cargo_links_with_exact_species_naturalized_host_match_both_endpoints":sum(
                x["each_island_exact_naturalized_HOSTS_species"] for x in relevant),
            "historical_cargo_links_with_same_butterfly_same_naturalized_host_both":sum(
                x["same_butterfly_exact_naturalized_HOSTS_species_both"] for x in relevant),
            "pair_rows":entries,
        }
    return {
        "schema":"chocho_historical_VOC_cargo_modern_GloNAF_exact_larval_host_naturalized_local_v1",
        "status":"HISTORICAL_COMMODITY_GENUS_X_MODERN_EXACT_HOST_PLANT_NATURALIZED_REGION_SCREEN_NOT_SPECIES_TRANSPORT",
        "source_sha256":{
            "historical_route_cargo_receipt":sha256(CARGO),
            "precommitted_commodity_to_genus":sha256(MAPLOCK),
            "precommitted_science_question":sha256(LOCK),
            "local_source_first_butterflies":sha256(CONTINUITY),
            "pinned_raw_HOSTS":sha256(HOSTS),
            "pinned_GloNAF_Zenodo_archive":sha256(GLO_NAF),
            "frozen_239_butterfly_panel":sha256(PANEL),
            "official_WGSRPD3_polygon_geojson":sha256(TDWG),
        },
        "original_Hosts_source_audit":host_audit,
        "original_GloNAF_taxon_aliases_audit":alias_audit,
        "original_GloNAF_region_source_audit":glonaf_audit,
        "original_GloNAF_geometry_issues":invalid,
        "source_island_to_region":SIX,
        "predeclared_area_overlap_thresholds":CUTOFFS,
        "primary_threshold":"0.05",
        "comparisons":comparisons,
        "interpretation_limits":[
            "VOC 17th century commodity categories cannot resolve exact historic shipped botanical species, live planting material, or specific island leg unloading.",
            "The local GloNAF naturalized alien plant binomial is a modern regional catalogue record, not a dated 17th-century or port-local record.",
            "The same modern HOSTS butterfly/plant species link is a global documented host relation, not observed larval feeding on those islands.",
            "JAP maps all of Japan rather than just Kyushu; it is not an island-scale alien species dataset.",
            "GloNAF and TDWG polygon intersection is a coarser spatial proxy with source coverage dependence, not a unique historic herbarium plant location.",
            "Data absence or no exact species-binomial alias is unresolved taxonomy/recording, not confirmed no transported plant or ecological host.",
            "GBIF 2011–2017 first indexed butterflies do not establish historical arrival, reproduction, establishment or site-specific host use.",
            "A matching triangle in VOC commodity GENUS, modern host SPECIES and naturalized flora SPECIES is not proof of historically shipped host propagule or butterfly.",
            "No p-values or enrichment tests against unobserved historical shipping routes are performed.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--pairs",type=Path,default=PAIRS)
    args=p.parse_args()
    result=analyze()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf8")
    primary=result["comparisons"]["0.05"]["pair_rows"]
    with args.pairs.open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(primary[0]))
        writer.writeheader();writer.writerows(primary)
    print(json.dumps({
        "source":result["status"],
        "primary":{k:v for k,v in result["comparisons"]["0.05"].items()
                   if k!="pair_rows"},
        "sensitivity":{k:{j:v for j,v in x.items() if j!="pair_rows"}
                       for k,x in result["comparisons"].items()},
    },indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Point-to-TDWG3 host-resource climatic-region mapping, all PR72 cases.

All 73 original positive source species×country cases. This only tests
observational overlap with *coarse botanical resource regions*, NOT verified
local larval host use, actual point climate or realized climatic niche entry.
Raw source coordinates/TDWG geometry never leave this process.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from shapely import make_valid
from shapely.geometry import Point, shape
from shapely.ops import unary_union
from shapely.strtree import STRtree

from audit_butterfly_73_point_source_20261027 import (
    MAX_PANEL_ROWS, MAX_QUERY_ROWS, PAGE, Reader, query_url, source_roster,
)
from analyze_butterfly_host_point_proximity_20261026 import valid_point, BASIS
from analyze_bce_augmented_plant_network_climate_null_20261019 import (
    source_augmented_model,
)
from analyze_butterfly_239_GBIF_occurrence_host_opportunity_20261019 import (
    climate_and_plant_candidates,
)

SHA = "408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
SCHEMA = "chocho_73_butterfly_point_tdwg_host_climate_20261027_v01"
CATEGORIES = (
    "introduced_analog_resource_region", "introduced_novel_resource_region",
    "no_documented_added_resource_region", "native_resource_region",
    "unmapped_point", "ambiguous_tdwg_region", "climate_not_modelled_region",
)


def polygonal(geo):
    geom = make_valid(shape(geo))
    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return geom
    if geom.geom_type == "GeometryCollection":
        polys = [x for x in geom.geoms if x.geom_type in ("Polygon", "MultiPolygon")]
        if polys:
            return unary_union(polys)
    raise ValueError("Nonpolygonal or invalid original TDWG3 region geometry")


class RegionLocator:
    def __init__(self, features, regions):
        grouped = defaultdict(list)
        accepted = set(regions)
        for f in features:
            name = f["properties"]["LEVEL3_COD"]
            if name not in accepted:
                continue
            g = polygonal(f["geometry"])
            if not g.is_valid or g.is_empty or g.area <= 0:
                raise ValueError("Invalid original TDWG feature " + name)
            grouped[name].append(g)
        if set(grouped) != accepted:
            raise ValueError("Original accepted TDWG3 region geometry incomplete")
        self.codes = list(regions)
        self.polygons = [unary_union(grouped[c]) for c in regions]
        self.tree = STRtree(self.polygons)
        if any(g.is_empty or not g.is_valid for g in self.polygons):
            raise ValueError("TDWG3 geometry union invalid")

    def lookup(self, lat, lon):
        pt = Point(float(lon), float(lat))
        matches = list(self.tree.query(pt, predicate="intersects"))
        if len(matches) == 0:
            return None, "unmapped_point"
        if len(matches) > 1:
            return None, "ambiguous_tdwg_region"
        return int(matches[0]), "matched"


def resource_class(index, species, bysp, *, q):
    data = bysp[species]
    if int(data["native"]) >> index & 1:
        return "native_resource_region"
    if not int(data["added"]) >> index & 1:
        return "no_documented_added_resource_region"
    mask = data["masks"].get(q)
    if mask is None:
        return "climate_not_modelled_region"
    return ("introduced_analog_resource_region" if int(mask) >> index & 1
            else "introduced_novel_resource_region")


def classify_point(lat, lon, species, locator, bysp, *, q=.05):
    index, status = locator.lookup(lat, lon)
    if status != "matched":
        return status
    return resource_class(index, species, bysp, q=q)


def fetch_and_measure(case, reader, remaining, locator, bysp):
    audit = Counter()
    qcounts = {str(q): Counter() for q in (.05, .10)}
    strict = {str(q): Counter() for q in (.05, .10)}
    source_hashes = []
    ids = set()
    try:
        data0, sha0 = reader.get(query_url(case, 0))
        n = int(data0["count"])
        if not 0 <= n <= 100000:
            raise ValueError("Impossible full-index count")
        if n > MAX_QUERY_ROWS or n > remaining:
            return {
                "full_source_status": "PREDECLARED_INDEX_VOLUME_CAP_NOT_EVALUABLE",
                "new_live_indexed_count": n,
                "butterfly": case["species"], "country": case["country"],
                "source_country_resource_class": case["resource_class"],
            }
        for offset in list(range(0, n, PAGE)) or [0]:
            data, digest = (data0, sha0) if offset == 0 else \
                reader.get(query_url(case, offset))
            if int(data["count"]) != n:
                raise ValueError("Source count drift within paginated query")
            records = data["results"]
            if not isinstance(records, list) or \
               len(records) != min(PAGE, max(0, n-offset)):
                raise ValueError("Source page incomplete")
            source_hashes.append({
                "offset": offset, "rows": len(records),
                "raw_body_sha256": digest,
            })
            for original in records:
                recid = original.get("key")
                if recid is None or recid in ids:
                    raise ValueError("Source GBIF key missing/duplicated")
                ids.add(recid)
                point, status = valid_point(original,
                    case["gbif_accepted_species_key"], case["country"])
                if point is not None and original.get("basisOfRecord") not in BASIS:
                    status, point = "unaccepted_basis", None
                audit[status] += 1
                if point is None:
                    continue
                for q in (.05, .10):
                    result = classify_point(point["lat"],point["lon"],
                                            case["species"],locator,bysp,q=q)
                    qcounts[str(q)][result] += 1
                    if point["strict"]:
                        strict[str(q)][result] += 1
        if len(ids) != n:
            raise ValueError("Original GBIF indexed key volume changed")
        eligible = audit["valid"]
        if any(sum(qcounts[str(q)].values()) != eligible for q in (.05, .10)):
            raise ValueError("Point map classification failed conservation")
        if any(k not in CATEGORIES for x in qcounts.values() for k in x):
            raise ValueError("Unregistered point resource category")
        return {
            "full_source_status": "COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING",
            "butterfly": case["species"], "country": case["country"],
            "source_country_resource_class": case["resource_class"],
            "gbif_species_key": case["gbif_accepted_species_key"],
            "original_GBIF_indexed_count": case["original_gbif_record_count"],
            "new_live_indexed_count": n,
            "source_count_delta": n - case["original_gbif_record_count"],
            "coordinate_valid_index_records": eligible,
            "coordinate_strict_le1000m_records": sum(strict["0.05"].values()),
            "original_country_area_added_resource_fraction":
                case["added_host_country_area_fraction"],
            "original_country_area_analog_resource_fraction":
                case["analog_host_country_area_fraction"],
            "original_country_area_novel_resource_fraction":
                case["novel_host_country_area_fraction"],
            "source_quality_disposition": dict(audit),
            "regional_climate_point_classes": {
                key: {c: int(value) for c,value in value.items()}
                for key, value in qcounts.items()
            },
            "strict_precision_point_classes": {
                key: {c: int(value) for c,value in value.items()}
                for key, value in strict.items()
            },
            "source_body_sha256_pages": source_hashes,
        }
    except Exception as e:
        return {
            "full_source_status": "INCOMPLETE_SOURCE_POINT_MAPPING",
            "butterfly": case["species"], "country": case["country"],
            "source_country_resource_class": case["resource_class"],
            "reason": f"{type(e).__name__}: {str(e)}"[:250],
            "source_body_sha256_pages": source_hashes,
        }


def group_report(cases, *, q="0.05"):
    merged = [r for r in cases if
              r["full_source_status"]=="COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING"]
    combined = Counter()
    nvalid = 0
    frac_cases = []
    for row in merged:
        observed = Counter(row["regional_climate_point_classes"][q])
        combined.update(observed)
        n = row["coordinate_valid_index_records"]
        nvalid += n
        if n > 0:
            overlap = observed["introduced_analog_resource_region"] + \
                      observed["introduced_novel_resource_region"]
            frac_cases.append(overlap/n)
    matched = (combined["introduced_analog_resource_region"]+
               combined["introduced_novel_resource_region"])
    return {
        "cases_source_complete": len(merged),
        "coordinate_valid_index_records": nvalid,
        "point_record_weighted_categories": dict(combined),
        "point_record_weighted_within_any_introduced_host_region": (
            matched/nvalid if nvalid else None),
        "species_country_equal_case_mean_introduced_resource_overlap": (
            sum(frac_cases)/len(frac_cases) if frac_cases else None),
        "cases_with_any_point_in_introduced_resource_region": sum(
            r["regional_climate_point_classes"][q].get(
                "introduced_analog_resource_region",0) +
            r["regional_climate_point_classes"][q].get(
                "introduced_novel_resource_region",0)>0 for r in merged),
        "cases_eligible_for_case_mean":len(frac_cases),
        "unmapped_or_ambiguous_indices":combined["unmapped_point"]+
                                       combined["ambiguous_tdwg_region"],
    }


def run(args):
    if not (args.source_root/"SOURCE_MANIFEST.json").is_file():
        raise ValueError("Original author-only source archive absent; fail closed")
    model,regions,climate,prov = source_augmented_model(args.source_root)
    bysp,classes,index,meta = climate_and_plant_candidates(
        model, regions, climate
    )
    if len(regions)!=365 or len(bysp)!=239 or \
       sum(x["masks"].get(.05) is not None for x in bysp.values()) != 235:
        raise ValueError("Source climate/geography model changed")
    source_geometry=(args.source_root/"inputs/level3.geojson")
    raw=source_geometry.read_bytes()
    geom_sha=hashlib.sha256(raw).hexdigest()
    features=json.loads(raw)["features"]
    loc=RegionLocator(features,regions)
    source=source_roster(args.roster)
    if not set(r["species"] for r in source["roster"]) <=set(bysp):
        raise ValueError("Butterfly names do not match source accepted host graph")
    reader=Reader()
    current=MAX_PANEL_ROWS
    rows=[]
    for case in sorted(source["roster"],key=lambda r:(r["species"],r["country"])):
        result=fetch_and_measure(case,reader,current,loc,bysp)
        rows.append(result)
        if result["full_source_status"]=="COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING":
            current-=result["new_live_indexed_count"]
    complete=len([r for r in rows if r["full_source_status"]==
                  "COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING"])
    group={
        key:{
            q:group_report([r for r in rows if
                            r["source_country_resource_class"]==key],q=q)
            for q in ("0.05","0.1")
        } for key in ("novel","analogue")
    }
    # Sensitivity source case exclusion is based on the known PR72 original GBIF
    # index counts (not selected after the observed co-location outcome).
    original_top2=sorted(source["roster"],
                         key=lambda r:-r["original_gbif_record_count"])[:2]
    drop={(r["species"],r["country"]) for r in original_top2}
    drop_view={
        q:group_report([r for r in rows if (r["butterfly"],r["country"])
                        not in drop],q=q)
        for q in ("0.05","0.1")
    }
    expected_counts=sum(r["new_live_indexed_count"] for r in rows
        if r["full_source_status"]=="COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING")
    report={
        "schema": SCHEMA,
        "status":("COMPLETE_73_POINT_TO_REGIONAL_RESOURCE_MAPPING"
                  if complete==73 else "INCOMPLETE_SOURCE_OR_REGION_MAPPING"),
        "source_original_author_only_zip_sha256":SHA,
        "source_tdwg_geojson_sha256":geom_sha,
        "original_source_geo_model": {
            "species": len(bysp),
            "regions":len(regions),
            "climate_evaluable_butterflies":235,
            "source_accepted_host_additions":1027,
            "source_CLIMATIC_resolution":"CHELSA 1981–2010 BIO1/BIO12 WGSRPD3 regional vectors",
        },
        "source_original_country_butterfly_cases":73,
        "original_indexed_butterfly_points":6076,
        "new_live_source_complete_cases":complete,
        "new_live_index_count_complete_cases":expected_counts,
        "source_original_accepted_count_matches":sum(
            r.get("source_count_delta")==0 for r in rows if
            r["full_source_status"]=="COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING"),
        "network_requests":reader.calls,
        "network_retries":reader.retries,
        "source_country_class_contrasts":group,
        "original_index_record_dominant_top_two_exclusion": {
            "source_selected_cases":original_top2,
            "remaining":drop_view,
        },
        "cases": rows,
        "geometry_inference_limit": (
            "An indexed butterfly point in a coarse WCVP introduced host "
            "resource region is NOT confirmed co-site host accessibility, "
            "plant use, or a measured climatic niche shift. Temperature/"
            "precipitation remain regional CHELSA proxies, not site climate."
        ),
        "no_raw_coordinates_or_raw_geometries_exported":True,
        "executed_utc":datetime.now(timezone.utc).isoformat(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n",
                           encoding="utf-8")
    print(json.dumps({
        "status":report["status"],
        "source_complete_cases":complete,
        "new_live_point_indices":expected_counts,
        "source_country_class_q05":{k:group[k]["0.05"] for k in group},
        "top_two_removed_q05":drop_view["0.05"],
        "requests":reader.calls
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--roster",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())

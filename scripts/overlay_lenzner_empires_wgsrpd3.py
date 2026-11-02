#!/usr/bin/env python3
"""Geographic overlay of pinned Lenzner et al. (2022) empire polygons on WGSRPD3 islands.

Classification uses measured GIS geometry rather than contemporary country
ownership, modern ports, butterfly occurrences or outcomes. Every region is
retained, including gaps and ambiguous historical exposure.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import shapefile
from pyproj import Transformer
from shapely import make_valid
from shapely.geometry import shape as geo_shape
from shapely.ops import transform, unary_union
from shapely.strtree import STRtree

EMPIRES=("Great_Britain","Spain","Portugal","Netherlands")


def readcsv(path):
    with path.open(encoding="latin-1" if path.name == "Empires.csv" else "utf-8", newline="") as h:
        yield from csv.DictReader(h)


def csvwrite(path, rows, fields):
    with path.open("w",encoding="utf-8",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    p=argparse.ArgumentParser()
    for name in ("wgsrpd3-geojson","q1-island-lookup","lenzner-empires","lenzner-shapefile","output-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    lookup = {r["area_code_l3"]: r for r in readcsv(a.q1_island_lookup) if r["class"]=="island"}
    if len(lookup)!=125: raise RuntimeError(f"Expected frozen 125 WGSRPD island areas; got {len(lookup)}")
    orig=json.loads(a.wgsrpd3_geojson.read_text(encoding="utf-8"))
    trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform

    island_shapes=defaultdict(list)
    for feature in orig["features"]:
        code=str((feature.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        if code in lookup and feature.get("geometry"):
            geom=make_valid(geo_shape(feature["geometry"]))
            if not geom.is_empty: island_shapes[code].append(transform(trans,geom))
    source_history=defaultdict(lambda:defaultdict(list))
    for row in readcsv(a.lenzner_empires):
        empire=row["colonizer"]
        if empire not in EMPIRES:continue
        key=str(row["OBJIDsic"]).strip()
        try: y1=int(row["from"]); y2=int(row["to"])
        except (TypeError,ValueError): continue
        if y2 < y1:raise RuntimeError(f"Invalid date interval in source: {row}")
        source_history[key][empire].append((y1,y2))

    reader=shapefile.Reader(str(a.lenzner_shapefile), encoding='latin-1')
    fieldnames=[f[0] for f in reader.fields[1:]]
    if "OBJIDsic" not in fieldnames:
        raise RuntimeError(f"Source shapefile lacks expected OBJIDsic: {fieldnames}")
    idx=fieldnames.index("OBJIDsic")
    geoms=[]
    ids=[]
    skipped=0
    for rec, shp in zip(reader.iterRecords(),reader.iterShapes()):
        key=str(rec[idx]).strip()
        if not key: continue
        try:
            geom=make_valid(geo_shape(shp.__geo_interface__))
            projected=transform(trans,geom)
            if projected.is_empty or projected.area<=0:continue
        except Exception:
            skipped+=1
            continue
        ids.append(key)
        geoms.append(projected)
    if len(geoms)<1000:raise RuntimeError(f"GloNAF shapefile read too few valid regions: {len(geoms)}")

    tree=STRtree(geoms)
    result=[]
    unmapped=[]
    for code,row in sorted(lookup.items()):
        shapes=island_shapes.get(code,[])
        if not shapes:
            unmapped.append(code)
            continue
        target=unary_union(shapes)
        region_area=target.area
        parts=[]
        for k in tree.query(target):
            j=int(k)
            g=geoms[j]
            if not g.intersects(target):
                continue
            intersection=target.intersection(g)
            if intersection.is_empty or intersection.area<=0:continue
            parts.append((ids[j],intersection))
        union=unary_union([x[1] for x in parts]) if parts else None
        covered=0 if union is None else min(region_area,union.area)
        coverage=0 if region_area<=0 else covered/region_area
        areas_by_empire=defaultdict(float)
        year_sum_by_empire=defaultdict(float)
        matched_ids=[]
        for key, geom in parts:
            matched_ids.append(key)
            for empire, times in source_history.get(key,{}).items():
                areas_by_empire[empire] += geom.area
                year_sum_by_empire[empire] += geom.area*sum(end-start+1 for start,end in times)
        fractions={k:(areas_by_empire[k]/covered if covered>0 else 0) for k in EMPIRES}
        historic_weighted_years={k:(year_sum_by_empire[k]/covered if covered>0 else 0) for k in EMPIRES}
        matches=sorted([k for k in EMPIRES if fractions[k]>=0.75])
        strict=(coverage>=0.80 and len(matches)==1 and
            max([fractions[k] for k in EMPIRES if k != matches[0]],default=0)<=0.2)
        assigned=matches[0] if strict else ("multi_empire" if len(matches)>1 else "ambiguous_or_uncovered")
        result.append({
            "area_code_l3":code,
            "region_name":row["region_name"],
            "level1":row["level1"],
            "island_geographic_area_km2":region_area/1e6,
            "lenzner_polygon_coverage":coverage,
            "lenzner_region_matches":len(set(matched_ids)),
            "lenzner_region_ids":";".join(sorted(set(matched_ids))),
            "classified_historical_empire":assigned,
            "strict_historical_empire":int(strict),
            **{empire+"_area_fraction":fractions[empire] for empire in EMPIRES},
            **{empire+"_weighted_exposure_years":historic_weighted_years[empire] for empire in EMPIRES}
        })
    if len(result)+len(unmapped)!=125:
        raise RuntimeError("Not all pinned Q1 islands represented")
    csvwrite(a.output_dir/"island_lenzner_empire_crosswalk.csv",result,list(result[0]))
    counts=defaultdict(int)
    for r in result:
        counts[r["classified_historical_empire"]]+=1
    payload={
      "schema":"chocho_lenzner2022_wgsrpd3_island_overlay_v0.1",
      "source":{"study_doi":"10.1038/s41559-022-01865-1",
                "lenzner_git_commit":"7aae72aadd39ae3827ece6935c6740dc0aace030",
                "wgsrpd_git_commit":"52da7828aba9d461dd133c27b3bd7a4407161f54"},
      "source_shape_regions_valid":len(geoms),
      "source_shape_skipped":skipped,
      "frozen_island_regions":125,
      "mapped_island_regions":len(result),
      "no_geometry_islands":unmapped,
      "high_confidence_classified_islands":sum(r["strict_historical_empire"] for r in result),
      "counts_by_assignment":dict(counts),
      "assumptions":[
        "Area overlaps are calculated in EPSG:6933 from the published GloNAF polygons, not geographic names.",
        "Historic exposure is multi-valued; islands with competing empires or poor geometric coverage are left ambiguous.",
        "Empire rule requires at least 80% area coverage, at least 75% of covered area in one empire and at most 20% in each other.",
        "Empire membership is a historic region-level association, never direct evidence of port arrivals or plant movements.",
        "Source GloNAF and WGSRPD mapping differences can be pronounced for small archipelagos, oceanic islands and mixed landmasses."
      ],
    }
    (a.output_dir/"empire_geographic_overlay_audit.json").write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2,ensure_ascii=False))
    return 0
if __name__=="__main__":
    raise SystemExit(main())

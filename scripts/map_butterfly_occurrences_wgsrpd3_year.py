#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from shapely.geometry import Point, shape
from shapely.strtree import STRtree


def load_level3(path: Path):
    payload=json.loads(path.read_text(encoding="utf-8"))
    geoms=[]; codes=[]
    for feature in payload.get("features",[]):
        props=feature.get("properties") or {}
        code=str(props.get("LEVEL3_COD") or "").strip()
        raw=feature.get("geometry")
        if code and raw:
            geoms.append(shape(raw)); codes.append(code)
    if len(geoms)<300 or len(set(codes))!=len(codes):
        raise RuntimeError("WGSRPD level-3 geometry/code drift")
    return geoms, tuple(codes)


def candidate_indices(tree: STRtree, geoms, point: Point):
    index_by_id={id(g):i for i,g in enumerate(geoms)}
    out=[]
    for candidate in tree.query(point):
        out.append(index_by_id[id(candidate)] if hasattr(candidate,"geom_type") else int(candidate))
    return out


def map_point(tree,geoms,codes,lon,lat):
    point=Point(float(lon),float(lat)); matches=[]
    for idx in candidate_indices(tree,geoms,point):
        if geoms[idx].covers(point):
            matches.append(codes[idx])
    if not matches:return ""
    return sorted(set(matches))[0]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--occurrences-csv",type=Path,required=True)
    ap.add_argument("--level3-geojson",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    args=ap.parse_args()

    geoms,codes=load_level3(args.level3_geojson); tree=STRtree(geoms)
    rows=[]; unmapped=0
    with args.occurrences_csv.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        required={"species","gbif_key","latitude","longitude","year"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("occurrence CSV schema drift")
        for row in reader:
            try:
                lon=float(row["longitude"]); lat=float(row["latitude"])
            except Exception:
                unmapped+=1; continue
            code=map_point(tree,geoms,codes,lon,lat)
            if not code:
                unmapped+=1; continue
            rows.append({
                "species":str(row["species"]).strip(),
                "gbif_key":str(row["gbif_key"]).strip(),
                "latitude":lat,
                "longitude":lon,
                "year":str(row["year"]).strip(),
                "basis_of_record":str(row.get("basis_of_record") or "").strip(),
                "dataset_key":str(row.get("dataset_key") or "").strip(),
                "wgsrpd3_code":code,
            })
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0]) if rows else ["species","gbif_key","latitude","longitude","year","basis_of_record","dataset_key","wgsrpd3_code"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"mapped_records":len(rows),"unmapped_records":unmapped},sort_keys=True))


if __name__=="__main__":
    main()

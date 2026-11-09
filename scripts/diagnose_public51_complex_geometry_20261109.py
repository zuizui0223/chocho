#!/usr/bin/env python3
"""Diagnose two large GSHHG sibling polygons without querying GBIF or altering source."""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile, math
from pathlib import Path
import geopandas as gpd
import shapely
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection
from shapely.ops import transform
from shapely.geometry.polygon import orient
from shapely import wkt
from pyproj import Transformer

SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"

def polygons(g):
    if g.is_empty:return []
    if isinstance(g,Polygon):return [g]
    if isinstance(g,(MultiPolygon,GeometryCollection)):
        return [p for z in g.geoms for p in polygons(z)]
    return []

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gshhg-zip",type=Path,required=True)
    ap.add_argument("--work",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    ap.add_argument("--islands",default="15,21")
    opt=ap.parse_args()
    b=opt.gshhg_zip.read_bytes()
    sha=hashlib.sha256(b).hexdigest()
    if sha!=SHA:raise ValueError("Historical source zip changed")
    opt.work.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        members=[p for p in z.namelist() if p.endswith((".shp",".shx",".dbf",".prj")) and
                 "/h/GSHHS_h_L1." in "/"+p]
        if len(members)<3:raise ValueError("High-resolution level 1 source missing")
        for p in members:z.extract(p,opt.work)
    df=gpd.read_file(next(opt.work.rglob("GSHHS_h_L1.shp")))
    df["landmass_id"]=df.sibling_id.astype(int).astype(str)
    tr=Transformer.from_crs(4326,6933,always_xy=True)
    area=lambda g:transform(tr.transform,g).area
    report={"schema":"original_GSHHG_complex_island_geometry_diagnosis_20261109",
            "shape_zip_SHA256":sha,"results":[],"no_GBIF_requests":True}
    for island in opt.islands.split(","):
        gs=df[df.landmass_id==island].geometry.to_numpy()
        if not len(gs):raise ValueError(f"Original island id {island} not found")
        combo=shapely.union_all(gs)
        parts=polygons(combo)
        summaries=[]
        for rank,g in enumerate(sorted(parts,key=area,reverse=True)[:10],1):
            original_area=area(g)
            if original_area<=0:raise ValueError("Nonpositive polygon area")
            diag={"rank":rank,"area_km2":original_area/1e6,
                  "bounds":[float(x) for x in g.bounds],
                  "original_WKT5_characters":len(wkt.dumps(g,rounding_precision=5,trim=True)),
                  "holes":len(g.interiors),
                  "exterior_nodes":len(g.exterior.coords),
                  "hole_nodes":sum(len(r.coords) for r in g.interiors),
                  "simplify":[]}
            if len(g.interiors):
                filled=Polygon(g.exterior)
                diag["fill_all_holes_area_ratio"]=area(filled)/original_area
            for tol in [0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.4]:
                q=g if tol==0 else g.simplify(tol,preserve_topology=True)
                if not q.is_valid or q.is_empty or not isinstance(q,Polygon):continue
                q=orient(q,sign=1.0)
                q5=wkt.dumps(q,rounding_precision=5,trim=True)
                a=area(wkt.loads(q5))/original_area
                diag["simplify"].append({"tol_degrees":tol,"characters":len(q5),
                  "area_retained":a,"within_3pct_area":.97<=a<=1.03,
                  "within_url_2300":len(q5)<=2300,"within_url_4300":len(q5)<=4300})
            summaries.append(diag)
        report["results"].append({
          "island_id":island,"original_rows":len(gs),"polygon_components":len(parts),
          "area_km2":area(combo)/1e6,"largest_10_components":summaries
        })
        print("ISLAND_GSHHG_GEOMETRY",island,json.dumps(report["results"][-1],ensure_ascii=False))
    opt.out_dir.mkdir(parents=True,exist_ok=True)
    (opt.out_dir/"original_complex_island_geometry_diagnostics.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Build a 126-WGSRPD3 island sea-GEOMETRY-only baseline, never ocean current.

Reuse original GSHHG 0.0625° land/water topology and 0.5° component-aware
marine graph from the 51-island barrier analysis. Attach independent fixed
official TDWG3 polygon boundary anchors to the sea component, NOT port
locations, then compute weighted shortest WATER PATHS (not exact Circuitscape
resistance, particle trajectories or seed arrival). Do not force unresolved
raster islands onto arbitrary water nodes.

Each TDWG3 archipelago region is represented by up to three deterministically
chosen coastline samples; source choice never consults WCVP similarities.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
from pyproj import Geod
from shapely.geometry import Point,shape
from shapely import points, covers
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

ORIGINAL_GRAPH_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
TDWG_SHA="52da7828aba9d461dd133c27b3bd7a4407161f54"
G=Geod(ellps="WGS84")
SNAP_MAX_KM=50.
MIN_SEPARATION_ANCHOR_KM=25.
MAX_ANCHORS=3

def sha(p):
    d=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):d.update(chunk)
    return d.hexdigest()

def rows(path):
    with path.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def polygons(g):
    if g.geom_type=="Polygon":return [g]
    if g.geom_type=="MultiPolygon":return sorted(list(g.geoms),key=lambda x:-x.area)
    return []

def nearest_wet_node(lon,lat,mask,fine_lats,fine_lons,component,size=16):
    if not (-85<lat<85 and math.isfinite(lon+lat)):return None
    nlat,nlon=mask.shape
    lo=((lon+180)%360)-180
    r=int(round((85-lat)/.0625-.5))
    c=int(round((lo+180)/.0625-.5))%nlon
    rr=np.arange(max(0,r-size),min(nlat,r+size+1),dtype=int)
    cc=(np.arange(c-size,c+size+1,dtype=int)%nlon).astype(int)
    grid=mask[np.ix_(rr,cc)]
    rpos,cpos=np.nonzero(grid>=0)
    if len(rpos)==0:return None
    rractual=rr[rpos];ccactual=cc[cpos]
    nodes=grid[rpos,cpos].astype(int)
    keep=(component[nodes]==0)
    # Source graph giant marine component is not assumed to be ID 0; take
    # the actual largest source node-component ID via code calling us.
    if not keep.any():return None
    rractual=rractual[keep];ccactual=ccactual[keep];nodes=nodes[keep]
    lonarr=fine_lons[ccactual];latarr=fine_lats[rractual]
    _,_,meter=G.inv(np.full(len(nodes),lo),np.full(len(nodes),lat),lonarr,latarr)
    ix=int(np.argmin(meter))
    if meter[ix]>SNAP_MAX_KM*1000:return None
    return {"node":int(nodes[ix]),"source_coast_candidate_lon":float(lo),
            "source_coast_candidate_lat":float(lat),
            "snapped_wet_longitude":float(lonarr[ix]),
            "snapped_wet_latitude":float(latarr[ix]),
            "snapped_to_wet_center_km":round(float(meter[ix]/1000),4),
            "fine_row":int(rractual[ix]),"fine_col":int(ccactual[ix])}

def assign_sources(code,polygon,mask,lat,lon,component,component_id):
    # Require source polygon itself to have rasterized GSHHG land centers
    # within a 1° local vicinity of sampled offshore-coast candidates.
    parts=polygons(polygon)
    origins=[]
    for g in parts[:5]:
        boundary=g.exterior
        for k in (0.05,0.20,0.40,0.60,0.80,0.95):
            sample=boundary.interpolate(k,normalized=True)
            if sample.is_empty:continue
            node=nearest_wet_node(sample.x,sample.y,mask,lat,lon,
                    np.where(component==component_id,0,1),size=16)
            if node is None:continue
            # Strong source identity gate: at least one GSHHG raster LAND
            # cell must be an exact point inside the TDWG source polygon
            # near the selected coast node; no inferred tiny islands.
            r,c=node["fine_row"],node["fine_col"]
            rr=np.arange(max(0,r-12),min(mask.shape[0],r+13))
            cc=np.arange(c-12,c+13)%mask.shape[1]
            cut=mask[np.ix_(rr,cc)]
            ir,ic=np.nonzero(cut<0)
            if not len(ir):continue
            pts=points(lon[cc[ic]],lat[rr[ir]])
            actual_land_inside=int(np.count_nonzero(covers(g,pts)))
            if actual_land_inside==0:continue
            node.update({"TDWG_code":code,"official_polygon_coast_verified":True,
              "nearby_GSHHG_raster_land_centers_in_same_TDWG_component":actual_land_inside,
              "original_TDWG_polygon_component_area_in_degree2":float(g.area)})
            if all(abs(G.inv(node["snapped_wet_longitude"],node["snapped_wet_latitude"],
                            v["snapped_wet_longitude"],v["snapped_wet_latitude"])[2])/1000
                   >=MIN_SEPARATION_ANCHOR_KM for v in origins):
                origins.append(node)
            if len(origins)>=MAX_ANCHORS:break
        if len(origins)>=MAX_ANCHORS:break
    return origins

def main(a):
    if sha(a.graph)!=ORIGINAL_GRAPH_SHA:
        raise RuntimeError("Original frozen GSHHG 0.0625° marine source SHA changed")
    pairs=rows(a.pair_input)
    expected_codes=sorted(set([r["source_region"] for r in pairs]+
                              [r["destination_region"] for r in pairs]))
    if len(expected_codes)!=126 or len(pairs)!=7875:
        raise RuntimeError("Source 126 accepted WCVP island cohort mismatch")
    raw=json.loads(a.tdwg3.read_text())
    t={}
    for f in raw["features"]:
        code=str((f.get("properties") or {}).get("LEVEL3_COD",""))
        if code in expected_codes:
            if code in t:raise RuntimeError("Original TDWG3 region duplicated")
            t[code]=shape(f["geometry"])
    if set(t)!=set(expected_codes):raise RuntimeError("TDWG original source geometry missing islands")
    with np.load(a.graph) as d:
        mask=d["fine_cell_to_node"];fine_lat=d["fine_latitude"];fine_lon=d["fine_longitude"]
        components=d["node_component"]
        row=d["node_block_row"];col=d["node_block_col"]
        ei=d["edge_i"];ej=d["edge_j"]
        block_lat=d["block_latitude"];block_lon=d["block_longitude"]
        source_comp_counts=np.bincount(components)
        giant=int(np.argmax(source_comp_counts))
        if mask.shape!=(2720,5760) or len(row)!=193516 or len(ei)!=372096:
            raise RuntimeError("Original GSHHG topology changed")
        if len(components)!=len(row) or len(col)!=len(row):
            raise RuntimeError("Marine graph node metadata missing")
        la1=block_lat[row[ei]];la2=block_lat[row[ej]]
        lo1=block_lon[col[ei]];lo2=block_lon[col[ej]]
        _,_,meter=G.inv(lo1,la1,lo2,la2)
        km=np.abs(np.asarray(meter,dtype=float))/1000.
        if np.any(km<=0) or np.any(~np.isfinite(km)):
            raise RuntimeError("Marine edge geometry invalid")
        # Same topology as the frozen circuit graph, DIFFERENT statistic
        # (edge geodesic weighted shortest path, not circuit resistance).
        graph=coo_matrix((np.concatenate([km,km]),
          (np.concatenate([ei,ej]),np.concatenate([ej,ei]))),
          shape=(len(row),len(row))).tocsr()
        del km
        attachments={}
        for code in expected_codes:
            attachments[code]=assign_sources(code,t[code],mask,fine_lat,
                                             fine_lon,components,giant)
        grouped={}
        for src in expected_codes:
            anchors=attachments[src]
            if len(anchors)==0:continue
            indices=sorted(set(x["node"] for x in anchors))
            d_all=dijkstra(graph,indices=indices,directed=False,min_only=True)
            grouped[src]={}
            for dst in expected_codes:
                if dst==src:continue
                eligible=[float(d_all[x["node"]]) for x in attachments[dst]
                          if math.isfinite(d_all[x["node"]])]
                grouped[src][dst]=min(eligible) if eligible else None
        geom=[]
        for r in pairs:
            ca,cb=r["source_region"],r["destination_region"]
            d1=grouped.get(ca,{}).get(cb)
            d2=grouped.get(cb,{}).get(ca)
            if d1 is not None and d2 is not None and
                 abs(d1-d2)>max(1.,1e-5*min(d1,d2)):
                raise RuntimeError("Unweighted marine graph is not symmetric")
            dsrc=d1 if d1 is not None else d2
            gdist=float(r["centroid_geodesic_km"])
            if dsrc is not None and (dsrc<=0 or not math.isfinite(dsrc)):
                raise RuntimeError("Unexpected zero or negative sea distance")
            rowpair={
              "source_region":ca,"destination_region":cb,
              "centroid_geodesic_km":gdist,
              "GSHHG_original_0p0625_sea_shortest_graph_km":round(dsrc,5) if dsrc is not None else "",
              "sea_shortest_over_centroid_geodesic_ratio":round(dsrc/gdist,5) if dsrc is not None and gdist>0 else "",
              "source_fixed_verified_offshore_anchors":len(attachments[ca]),
              "dest_fixed_verified_offshore_anchors":len(attachments[cb]),
              "source_sea_geometry_coverage_status":"BOTH_ANCHORED" if dsrc is not None else "NO_COMMON_VERIFIED_SEA_ANCHOR",
              "is_drift_or_lagrangian_particle_transport":False,
              "source_graph_GSHHG_sha256":ORIGINAL_GRAPH_SHA,
            }
            geom.append(rowpair)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"WCVP126_GSHHG_source_coastal_anchors_and_sea_geometry_pairs.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(geom[0]));w.writeheader();w.writerows(geom)
    nusable=sum(v["source_sea_geometry_coverage_status"]=="BOTH_ANCHORED" for v in geom)
    audit={
      "schema":"chocho_WCVP126_accepted_flora_source_locked_GSHHG_0p0625_geometry_shortest_water_path_v1",
      "status":"SOURCE_LOCKED_MARINE_SHORTEST_PATH_GEOMETRY_NOT_CIRCUIT_RESISTANCE_OR_OCEAN_CURRENTS",
      "source_original_GSHHG_graph_SHA256":ORIGINAL_GRAPH_SHA,
      "TDWG3_original_git_commit":TDWG_SHA,
      "graph_node_count":193516,
      "graph_undirected_edge_count":372096,
      "giant_source_marine_component_id":giant,
      "giant_source_graph_nodes":int(source_comp_counts[giant]),
      "source_original_126_regions":len(expected_codes),
      "source_126_regions_with_at_least_one_verified_coastal_water_anchor":
          sum(len(v)>0 for v in attachments.values()),
      "source_region_missing_verified_coastal_anchors":
          [code for code in expected_codes if not attachments[code]],
      "source_valid_original_island_pair_geometric_shortest_sea":nusable,
      "source_original_island_pairs_total":len(geom),
      "max_fixed_geom_coastal_anchors_per_region":MAX_ANCHORS,
      "anchor_minimum_separation_km":MIN_SEPARATION_ANCHOR_KM,
      "max_snap_from_official_TDWG_polygon_coast_km":SNAP_MAX_KM,
      "source_original_WPI_ports_not_used_as_natural_geographic_terminals":True,
      "coastal_anchors_by_TDWG_region":attachments,
      "limits":[
       "This reuses original 51-island GSHHG ocean graph topology but updates the calculation to 126 independently placed official TDWG3 coastal terminals.",
       "The scalar is the shortest sea path in 0.5-degree block edge geodesic lengths, not exactly 0.0625-degree Euclidean shortest path and NOT Circuitscape effective resistance.",
       "For archipelagos, at most three deterministic shoreline points cannot represent all island coasts. Raster land presence and nearest marine node gate may exclude tiny islands.",
       "The marine graph represents water connectivity only, without ocean current direction, speed, season, wind drift, ships or seeds.",
       "A coast-attachment failure is unknown model coverage, NOT a non-existent island or biologically impossible dispersal.",
       "Future ocean current transport must separately use actual directional NOAA velocities at the same region-source terminals and lag horizons.",
      ]
    }
    (a.out_dir/"GSHHG126_sea_geometry_source_coverage_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n")
    print("GSHHG126_SEA_GEOMETRY_BASELINE",json.dumps({
        "source_original_regions":audit["source_original_126_regions"],
        "valid_coast_region_anchors":audit["source_126_regions_with_at_least_one_verified_coastal_water_anchor"],
        "valid_pair_sea_distances":nusable,"original_pair_count":len(geom),
        "excluded_missing_coast":audit["source_region_missing_verified_coastal_anchors"],
        "status":audit["status"]
    }),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    for name in ("graph","pair-input","tdwg3","out-dir"):
        ap.add_argument("--"+name,type=Path,required=True)
    main(ap.parse_args())

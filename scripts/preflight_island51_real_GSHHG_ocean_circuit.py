#!/usr/bin/env python3
"""Global 1 degree OCEAN grid topology and 51 original landmass port-centroid
snap audit. Uses checksum-verified original GSHHG v2.3.7 coarse level 1 polygons
for grid-center water mask and exact source WPI2019 archipelago geocentroids.

This is a geometry-only QA. Do NOT use synthetic inferred island positions,
treat land cells as navigable water, or report any AIS/wind pathway before
independent geolocation and sea connectivity gates pass.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path

import geopandas as gpd
import numpy as np
import scipy.sparse as ss
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import shapely
from shapely.strtree import STRtree

SOURCE_GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
N_ISLANDS=51
LAT=np.arange(84.5,-85.0,-1.0)
LON=np.arange(-179.5,180.0,1.0)
R=6371.0088


def digest(path):
    v=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024*4),b""):v.update(chunk)
    return v.hexdigest()


def to_xyz(lon,lat):
    x=np.deg2rad(lon);y=np.deg2rad(lat)
    return np.stack([np.cos(y)*np.cos(x),np.cos(y)*np.sin(x),np.sin(y)],axis=-1)


def land_mask_from_GSHHG(original_zip):
    checksum=digest(original_zip)
    if checksum!=SOURCE_GSHHG_SHA:
        raise ValueError("Pinned GSHHG exact original source SHA256 mismatch")
    uri=f"zip://{original_zip.resolve()}!GSHHS_shp/c/GSHHS_c_L1.shp"
    geodata=gpd.read_file(uri)
    if not geodata.crs:raise ValueError("No source coastline CRS")
    geodata=geodata.to_crs(4326)
    polygons=geodata.geometry[geodata.geometry.notna()].to_numpy()
    if not len(polygons):raise ValueError("No source GSHHG coarse coastline land polygons")
    landtree=STRtree(polygons)
    llon,llat=np.meshgrid(LON,LAT)
    pts=shapely.points(llon.ravel(),llat.ravel())
    q=landtree.query(pts,predicate="intersects")
    land=np.zeros(len(pts),dtype=bool)
    if q.shape[1]>0:land[q[0]]=True
    return land.reshape(len(LAT),len(LON)),{
         "GSHHG_source_SHA256":checksum,
         "source_1deg_water_center_mask_method":"original GSHHS_c_L1 point intersects polygon at each 1deg midpoint",
         "coarse_polygon_features":len(polygons),
         "water_centers":int(np.sum(~land)),
         "land_centers":int(np.sum(land)),
         "coast_resolution_limit":"1 degree sea-center mask misses channels narrower than approximately 1 degree"
    }


def connect_ocean(sea,neighbors=8):
    h,w=sea.shape
    if (h,w)!=(170,360):raise ValueError("Not a fixed source 170×360 latlon grid")
    if neighbors not in (4,8):raise ValueError("Supported topology 4/8 only")
    ids=np.full(sea.shape,-1,dtype=np.int32)
    ids[sea]=np.arange(int(sea.sum()))
    rows=[];cols=[]
    dirs=[(1,0),(0,1)]
    if neighbors==8:dirs.extend([(1,1),(1,-1)])
    for dy,dx in dirs:
        for y in range(h):
            yy=y+dy
            if yy>=h:continue
            for x in range(w):
                if not sea[y,x]:continue
                xx=(x+dx)%w
                if not sea[yy,xx]:continue
                # A diagonal must not pass through the corner of a land
                # cell. The side-cell check uses the same longitude wrap.
                if dy and dx and not (sea[y,xx] and sea[yy,x]):continue
                i=int(ids[y,x]);j=int(ids[yy,xx])
                rows.extend((i,j));cols.extend((j,i))
    graph=ss.csr_matrix((np.ones(len(rows),dtype=np.uint8),(rows,cols)),
                        shape=(int(sea.sum()),)*2)
    n,components=connected_components(graph,directed=False)
    counts=np.bincount(components)
    return ids,graph,components,{
        "neighbor_definition":neighbors,
        "diagonal_corner_rule":"both orthogonal side cells must be sea" if neighbors==8 else "no diagonal edges",
        "marine_cells":int(sea.sum()),
        "undirected_edges":len(rows)//2,
        "number_of_sea_components":int(n),
        "largest_component_cells":int(np.max(counts)),
        "largest_component_fraction_of_ocean_cells":float(np.max(counts)/len(components)),
    }


def snap_island_centroids(source,sea,ids,components,*,distance_limit_km=220):
    if source.get("status")!="ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES":
        raise ValueError("Original 51 island coordinates not verified")
    points=source.get("verified_source_51_mean_port_points",{})
    if len(points)!=N_ISLANDS:raise ValueError("Exact original 51-source island target lost")
    llon,llat=np.meshgrid(LON,LAT)
    ocean_xyz=to_xyz(llon[sea],llat[sea])
    ocean_yx=np.argwhere(sea)
    tree=cKDTree(ocean_xyz)
    obs={}
    for key,p in sorted(points.items(),key=lambda x:int(x[0])):
        xyz=to_xyz(float(p["longitude"]),float(p["latitude"]))
        chord,k=tree.query(xyz,k=1)
        y,x=map(int,ocean_yx[int(k)])
        node=int(ids[y,x])
        dist=float(2*R*np.arcsin(np.clip(chord/2,0,1)))
        if dist>distance_limit_km:
            obs[key]={"status":"FAILED_NO_MARINE_CENTER_WITHIN_DISTANCE_LIMIT",
                 "longitude":float(p["longitude"]),"latitude":float(p["latitude"]),
                 "nearest_ocean_grid_node":node,"nearest_displacement_km":dist}
        else:
            obs[key]={
              "status":"SOURCE_SITE_SNAPPED_TO_NEAREST_SEA_CELL",
              "original_longitude":float(p["longitude"]),
              "original_latitude":float(p["latitude"]),
              "ocean_grid_longitude":float(LON[x]),"ocean_grid_latitude":float(LAT[y]),
              "ocean_grid_node":node,"raster_row":y,"raster_col":x,
              "displacement_km":dist,"sea_component_id":int(components[node])
            }
    node_islands={}
    for key,item in obs.items():
        if item["status"].startswith("SOURCE_SITE_"):
            node_islands.setdefault(str(item["ocean_grid_node"]),[]).append(key)
    duplicate_nodes={node:keys for node,keys in node_islands.items() if len(keys)>1}
    for keys in duplicate_nodes.values():
        for key in keys:
            obs[key]["shared_ocean_grid_node_with"]=[other for other in keys if other!=key]
    ok=[x for x in obs.values() if x["status"].startswith("SOURCE_SITE_")]
    groups={}
    for x in ok:
        key=str(x["sea_component_id"])
        groups[key]=groups.get(key,0)+1
    return obs,{
       "islands_original":len(points),
       "islands_snapped":len(ok),
       "islands_failed":len(points)-len(ok),
       "snapping_method":"independent nearest sea cell; no displacement to force distinct island nodes",
       "unique_marine_nodes":len(node_islands),
       "duplicate_marine_node_groups":duplicate_nodes,
       "all_51_distinct_marine_nodes":len(ok)==51 and len(node_islands)==51,
       "snapped_components_counts":groups,
       "all_51_in_one_marine_component":len(ok)==51 and len(groups)==1,
       "max_snap_distance_km":max((x["displacement_km"] for x in ok),default=None),
       "median_snap_distance_km":float(np.median([x["displacement_km"] for x in ok])) if ok else None,
       "distance_limit_km":distance_limit_km,
       "land_mask_not_based_on_AIS_butterfly_or_plant_outcomes":True
    }


def run(original_zip,source,*,neighbors):
    land,info=land_mask_from_GSHHG(original_zip)
    ids,graph,comp,meta=connect_ocean(~land,neighbors)
    nodes,nodesmeta=snap_island_centroids(source,~land,ids,comp)
    return {"schema":"chocho_GSHHG237_51_WPI_source_1deg_ocean_topology_and_snapping_v1",
        "status":"VALID_1DEG_OCEAN_DISCRETE_CIRCUIT_SOURCE_TOPOLOGY" if nodesmeta["all_51_in_one_marine_component"] and nodesmeta["all_51_distinct_marine_nodes"]
            else "INCOMPLETE_1DEG_OCEAN_CONNECTIVITY_DO_NOT_FIT_ALL_51",
        "GSHHG_original":info,"graph":meta,
        "WPI_centroid_marine_snapping":nodesmeta,
        "island_points":nodes,
        "notes":[
          "1 degree seacell-centroid connectivity is extremely coarse; land and strait geometry can disconnect true oceanic routes.",
          "The original 51 islands are defined by existing biological data; no point excluded silently.",
          "Island centroids are independently assigned to the nearest sea cell; shared coarse cells are retained and reported, never separated by moving an island.",
          "A 1-degree sea grid cannot resolve the shipping network itself and cannot assess current speed or wind direction.",
          "If not all islands are connected, do not fill missing marine routes with fabricated resistances or silently change sample."
        ]}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("gshhg-zip","verified-islands","out"):p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--neighbors",type=int,choices=(4,8),default=8)
    a=p.parse_args()
    source=json.loads(a.verified_islands.read_text())
    result=run(a.gshhg_zip,source,neighbors=a.neighbors)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("ISLAND51_ORIGINAL_GSHHG_1DEG_OCEAN_TOPOLOGY",json.dumps({
        "status":result["status"],"mask":result["GSHHG_original"],
        "graph":result["graph"],"anchor":result["WPI_centroid_marine_snapping"]
    },sort_keys=True),flush=True)
if __name__=="__main__":main()

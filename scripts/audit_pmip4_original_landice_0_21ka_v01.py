#!/usr/bin/env python3
"""Verified original PMIP4 ICE-6G_C global physical land & ice masks, 0 and 21 ka.

Scientific scope: passable terrestrial *substrate* only, not historical plant
occupancy, climatic habitat corridors or genetic barriers.
"""
import argparse,hashlib,io,json
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file
from scipy.ndimage import label

KNOWN={0:("5ba8c0df48a12132238ffcbf675578e0","518bdfa9c918cb17090a4380c494572ffd011a2a5b547dab57327faf7ade415a"),
       21:("9443668583800099ce9e865be90e4ad6","3decf87e11a2d2afe65a7358d904214aea2e5c59ec92ab8227f36070beb98e94")}
RADIUS_KM=6371.0088

def get_source(path,year):
    data=Path(path).read_bytes()
    if (hashlib.md5(data).hexdigest(),hashlib.sha256(data).hexdigest())!=KNOWN[year]:
        raise ValueError("published PMIP4 0/21ka source hash mismatch")
    with netcdf_file(io.BytesIO(data),mmap=False) as nc:
        lat=np.asarray(nc.variables["lat"][:],float).copy()
        lon=np.asarray(nc.variables["lon"][:],float).copy()
        land=np.asarray(nc.variables["sftlf"][:],float).copy()
        ice=np.asarray(nc.variables["sftgif"][:],float).copy()
    if (lat.shape!=(1080,) or lon.shape!=(2160,) or
        land.shape!=ice.shape!=(1080,2160) or
        not np.allclose(np.diff(lat),1/6,atol=1e-4) or
        not np.allclose(np.diff(lon),1/6,atol=1e-4) or
        abs(float(lat[0])+90-1/12)>1e-4 or abs(float(lon[0]))>1e-5):
        raise ValueError("PMIP4 geographic grid and orientation drift")
    if (not np.isfinite(land).all() or not np.isfinite(ice).all() or
        not np.isin(land,[0.,100.]).all() or not np.isin(ice,[0.,100.]).all()):
        raise ValueError("nonbinary or missing published PMIP4 land/ice mask")
    return lat,lon,land==100,ice==100

def icefree(land,ice):
    if land.shape!=ice.shape:raise ValueError("mask shape mismatch")
    return land.astype(bool)&~ice.astype(bool)

def cell_area(lat):
    lat=np.deg2rad(np.asarray(lat,float))
    half=np.deg2rad(1/12)
    return RADIUS_KM**2*np.deg2rad(1/6)*(np.sin(lat+half)-np.sin(lat-half))

def area(mask,area_per_row):
    return float(np.sum(np.count_nonzero(mask,axis=1)*area_per_row))

def component_ids(mask):
    """Four-neighbor physical geography; join cells across longitude 0/360."""
    mask=np.asarray(mask,bool)
    if mask.ndim!=2:raise ValueError("2D mask required")
    lab,count=label(mask,np.array([[0,1,0],[1,1,1],[0,1,0]],int))
    parent=np.arange(count+1,dtype=np.int32)
    def find(k):
        while parent[k]!=k:
            parent[k]=parent[parent[k]]
            k=parent[k]
        return k
    for row in np.flatnonzero(mask[:,0]&mask[:,-1]):
        a=find(int(lab[row,0]));b=find(int(lab[row,-1]))
        if a!=b:parent[max(a,b)]=min(a,b)
    roots=np.arange(count+1,dtype=np.int32)
    for k in range(1,count+1):roots[k]=find(k)
    return roots[lab],int(np.unique(roots[1:]).size)

def nearest_grid_indices(lat,lon):
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    if (lat.shape!=lon.shape or not np.isfinite(lat).all() or not np.isfinite(lon).all()
        or np.any(abs(lat)>90) or np.any(lon < -180) or np.any(lon>360)):
        raise ValueError("invalid supported point coordinates")
    row=np.clip(np.rint((lat+90-1/12)*6).astype(int),0,1079)
    col=np.rint(np.mod(lon,360)*6).astype(int)%2160
    return row,col

def summary(year,land,ice,lat):
    a=cell_area(lat);usable=icefree(land,ice)
    return {"age_ka":year,"land_cells":int(land.sum()),
        "land_area_km2":area(land,a),
        "land_ice_cells":int(np.count_nonzero(land&ice)),
        "land_ice_area_km2":area(land&ice,a),
        "floating_ice_cells":int(np.count_nonzero(~land&ice)),
        "ice_free_land_cells":int(usable.sum()),
        "ice_free_land_area_km2":area(usable,a)}

def transitions(land0,ice0,land21,ice21,lat):
    now=icefree(land0,ice0);old=icefree(land21,ice21);a=cell_area(lat)
    masks={
      "retained_ice_free_land":now&old,
      "ice_free_land_lost_0_to_21":now&~old,
      "ice_free_land_gained_0_to_21":~now&old,
      "new_emerged_land_from_0_ocean":~land0&land21,
      "new_ice_free_land_from_0_ocean":~land0&old,
      "land0_to_ice21":now&land21&ice21,
    }
    return {name:{"cells":int(v.sum()),"km2":area(v,a)} for name,v in masks.items()}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--plant-score-npz",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    q=p.parse_args();q.out.mkdir(parents=True,exist_ok=True)
    lat,lon,l0,i0=get_source(q.source_dir/"I6_C.VM5a_10min.0.nc",0)
    other_lat,other_lon,l21,i21=get_source(q.source_dir/"I6_C.VM5a_10min.21.nc",21)
    if not np.array_equal(lat,other_lat) or not np.array_equal(lon,other_lon):
        raise ValueError("published 0/21 ka geographical grids not aligned")
    free0,free21=icefree(l0,i0),icefree(l21,i21)
    comp0,n0=component_ids(free0);comp21,n21=component_ids(free21)
    with np.load(q.plant_score_npz,allow_pickle=False) as src:
        if src["scores_present"].shape!=src["scores_21ka"].shape!=(518,6337):
            raise ValueError("original 518 plant model 6337 site array changed")
        site_ids=src["global_support_indices"]
        r,c=nearest_grid_indices(src["support_latitude"],src["support_longitude"])
        x=comp0[r,c];y=comp21[r,c]
    land_summary={
      "schema":"chocho_ICE6GC0_21ka_source_verified_physical_mask_v01",
      "status":"PHYSICAL_LAND_ICE_TOPOLOGY_NOT_21KA_PLANT_HABITAT_OR_GENETICS",
      "source_md5":{"0ka":KNOWN[0][0],"21ka":KNOWN[21][0]},
      "source_sha256":{"0ka":KNOWN[0][1],"21ka":KNOWN[21][1]},
      "grid_cells":int(free0.size),
      "current":summary(0,l0,i0,lat),"LGM":summary(21,l21,i21,lat),
      "transitions":transitions(l0,i0,l21,i21,lat),
      "physical_land_component_count_4_neighbor_dateline_connected":{"0ka":n0,"21ka":n21},
      "native_WCVP_global_current_support_sites":int(len(r)),
      "support_ice_free_land_0ka":int(np.count_nonzero(x)),
      "support_ice_free_land_21ka":int(np.count_nonzero(y)),
      "support_retained_ice_free_land":int(np.count_nonzero((x>0)&(y>0))),
      "support_newly_blocked_0_to_21":int(np.count_nonzero((x>0)&(y==0))),
      "support_newly_ice_free_0_to_21":int(np.count_nonzero((x==0)&(y>0))),
      "limits":["PMIP4 source is reconstructed point-cell land/ice and uses nearest 10-minute cells, not proven exact coastlines",
        "The original 6337 WCVP support sites omit many newly exposed LGM continental-shelf locations",
        "A connected ice-free landmass does not establish a continuous climatically suitable host pathway",
        "No ancient plant occupancy, plant dispersal or insect genetics entered the physical source audit"]
    }
    np.savez_compressed(q.out/"physical_landice_by_6337_support_sites.npz",
       global_support_indices=site_ids,
       source_age0_ice_free=(x>0).astype("u1"),source_age21_ice_free=(y>0).astype("u1"),
       source_age0_land_component=x.astype("i4"),source_age21_land_component=y.astype("i4"))
    (q.out/"summary.json").write_text(json.dumps(land_summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(land_summary,sort_keys=True))
if __name__=="__main__":main()

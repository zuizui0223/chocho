"""Pure potential-climate habitat component metrics with a periodic longitude seam.

Connectivity is topology of discretized *climate-suitable ice-free land*, NOT
occupied plant habitat or dispersal history. Four neighbors are primary;
eight neighbors are deliberately permissive sensitivity.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import label

EARTH_R_KM=6371.0088

def regular_cell_areas(latitude_centers_deg, resolution_deg=0.5):
    lat=np.deg2rad(np.asarray(latitude_centers_deg,float))
    half=np.deg2rad(resolution_deg/2)
    if (lat.ndim!=1 or not np.isfinite(lat).all() or np.any(abs(lat)>np.pi/2+1e-8)
        or resolution_deg<=0):
        raise ValueError("invalid spherical raster latitudes")
    return (EARTH_R_KM**2*np.deg2rad(resolution_deg)
            *(np.sin(np.minimum(lat+half,np.pi/2))-np.sin(np.maximum(lat-half,-np.pi/2))))

def components(mask, connectivity=4):
    z=np.asarray(mask,bool)
    if z.ndim!=2 or z.shape[1]<3 or z.shape[0]<2 or connectivity not in (4,8):
        raise ValueError("invalid 2D potential habitat source/mode")
    neighborhood=(np.array([[0,1,0],[1,1,1],[0,1,0]],int)
                  if connectivity==4 else np.ones((3,3),int))
    labels,n=label(z,neighborhood)
    if n==0:return np.zeros(z.shape,dtype=np.int32),0
    parents=np.arange(n+1,dtype=np.int32)
    def find(k):
        while parents[k]!=k:
            parents[k]=parents[parents[k]];k=parents[k]
        return k
    def union(a,b):
        a=find(int(a));b=find(int(b))
        if a!=b:parents[max(a,b)]=min(a,b)
    rows=np.flatnonzero(z[:,0])
    for row in rows:
        for dr in ((0,) if connectivity==4 else (-1,0,1)):
            other=row+dr
            if 0<=other<len(z) and z[other,-1]:
                union(labels[row,0],labels[other,-1])
    lookup=np.arange(n+1,dtype=np.int32)
    for k in range(1,n+1):lookup[k]=find(k)
    labels=lookup[labels]
    return labels,int(np.unique(lookup[1:]).size)

def graph_metrics(mask,latitude_centers_deg,connectivity=4):
    mask=np.asarray(mask,bool)
    lat=np.asarray(latitude_centers_deg,float)
    if mask.ndim!=2 or len(lat)!=mask.shape[0]:
        raise ValueError("mask and latitude row coordinates misaligned")
    area=regular_cell_areas(lat)
    labels,n=components(mask,connectivity=connectivity)
    cells=int(mask.sum())
    area_total=float(np.sum(mask.sum(axis=1)*area))
    if n==0:
        return {"components":0,"suitable_grid_cells":0,"suitable_km2":0.,
                "largest_component_cells":0,"largest_component_km2":0.,
                "fraction_cells_in_largest_component":None,
                "fraction_area_in_largest_component":None}
    positive=labels[mask]
    weights=np.broadcast_to(area[:,None],mask.shape)[mask]
    counts=np.bincount(positive)
    area_by_comp=np.bincount(positive,weights=weights)
    largest_label=int(np.argmax(area_by_comp[1:]))+1
    return {"components":n,"suitable_grid_cells":cells,
            "suitable_km2":area_total,
            "largest_component_cells":int(counts[largest_label]),
            "largest_component_km2":float(area_by_comp[largest_label]),
            "fraction_cells_in_largest_component":float(counts[largest_label]/cells),
            "fraction_area_in_largest_component":float(area_by_comp[largest_label]/area_total)}

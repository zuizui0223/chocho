from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,geopandas as gpd,shapely
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
lock=O/'cell_land_overlap_lock.json';assert lock.exists()
c=pd.read_parquet(W/'cell_geography.parquet');assert len(c)==5282 and c.cell.is_unique
source=next(W.rglob('GSHHS_h_L1.shp'))
land=gpd.read_file(source).to_crs(6933)
geoms=land.geometry.to_numpy();bad=~shapely.is_valid(geoms)
geoms[bad]=shapely.make_valid(geoms[bad]);assert shapely.is_valid(geoms).all()
tree=shapely.STRtree(geoms);ids=land.sibling_id.astype(str).to_numpy()
rows=[];parts=[]
for i,r in enumerate(c.itertuples(index=False)):
    x,y=map(int,r.cell.split(':'));box=shapely.box(x*10000,y*10000,(x+1)*10000,(y+1)*10000)
    ix=tree.query(box,predicate='intersects');by={}
    for j in ix:
        z=shapely.intersection(box,geoms[j])
        if z.area>0:by.setdefault(ids[j],[]).append(z)
    areas={k:shapely.union_all(v).area for k,v in by.items()}
    total=shapely.union_all([z for v in by.values() for z in v]).area if by else 0.
    assert -1e-4<=total<=1e8+1e-4
    assert abs(sum(areas.values())-total)<1.,'Overlapping different landmass polygons'
    dominant=max(areas,key=areas.get) if areas else None
    centre=str(r.landmass_id) if pd.notna(r.landmass_id) else None
    centre_fraction=areas.get(centre,0)/total if total and centre is not None else np.nan
    rows.append(dict(cell=r.cell,centre_landmass_id=centre,land_fraction=total/1e8,
        intersecting_landmasses=len(areas),dominant_landmass_id=dominant,
        dominant_share_of_land=areas[dominant]/total if total else np.nan,
        centre_share_of_land=centre_fraction,
        centre_is_dominant=centre==dominant if centre and dominant else None,
        area_supported=bool(total>=5e7 and np.isfinite(centre_fraction) and centre_fraction>=.9)))
    parts.extend(dict(cell=r.cell,landmass_id=k,intersection_km2=v/1e6) for k,v in areas.items())
    if (i+1)%1000==0:print('CELLS',i+1,flush=True)
d=pd.DataFrame(rows);d.to_csv(O/'cell_land_overlap.csv',index=False)
pd.DataFrame(parts).to_parquet(W/'cell_land_overlap_parts.parquet',index=False)
audit=dict(status='AREA_OVERLAP_AUDITED_NOT_OCCURRENCE_ASSIGNMENT',cells=len(d),
    centre_unknown=int(d.centre_landmass_id.isna().sum()),multiple_landmasses=int(d.intersecting_landmasses.gt(1).sum()),
    no_land=int(d.land_fraction.eq(0).sum()),area_supported=int(d.area_supported.sum()),
    known_centre_not_dominant=int(d.centre_is_dominant.eq(False).sum()),repaired_source_geometries=int(bad.sum()),
    lock_sha256=hashlib.sha256(lock.read_bytes()).hexdigest(),shoreline_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    limits=json.loads(lock.read_text())['limits'])
(O/'cell_land_overlap_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit))

"""Acquire DSM pixels for fixed route samples; retain all missing values."""
from pathlib import Path
import os
os.environ.setdefault('GDAL_HTTP_MAX_RETRY','2');os.environ.setdefault('GDAL_HTTP_TIMEOUT','45')
import json,hashlib,argparse,concurrent.futures
import pandas as pd,numpy as np,requests,rasterio,planetary_computer
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
P=W/'route_terrain_by_tile';P.mkdir(exist_ok=True)
D=W/'route_terrain_values';D.mkdir(exist_ok=True)
parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int,default=0);args=parser.parse_args()
design=json.loads((O/'route_terrain_design_audit.json').read_text())
assert design['status']=='ROUTE_SAMPLE_DESIGN_COMPLETE_ELEVATIONS_NOT_ACQUIRED'
lock_digest=hashlib.sha256((O/'route_terrain_design_lock.json').read_bytes()).hexdigest();assert design['lock_sha256']==lock_digest
manifest=W/'route_terrain_partition_audit.json'
if manifest.exists():
    cached=json.loads(manifest.read_text());assert cached['lock_sha256']==lock_digest and cached['source_mapping_sha256']==design['source_mapping_sha256']
else:
    total=0
    for part in sorted((W/'route_terrain_design_parts').glob('part_*.parquet')):
        points=pd.read_parquet(part);total+=len(points)
        for tile,q in points.groupby('tile'):
            dest=P/tile;dest.mkdir(exist_ok=True)
            q[['route_index','longitude','latitude']].to_parquet(dest/part.name,index=False)
    assert total==design['planned_land_samples']
    manifest.write_text(json.dumps(dict(lock_sha256=lock_digest,source_mapping_sha256=design['source_mapping_sha256'],points=total)))
tiles=pd.read_csv(O/'route_terrain_required_tiles.csv')
assert len(tiles)==design['required_tiles']
jobs=list(tiles.itertuples(index=False));jobs=jobs[:args.limit] if args.limit else jobs
def acquire(job):
    tile=job.tile;receipt_path=D/(tile+'.json');data_path=D/(tile+'.parquet')
    q=pd.read_parquet(P/tile)
    assert len(q)==job.route_land_sample_points
    point_hash=hashlib.sha256(pd.util.hash_pandas_object(q,index=False).to_numpy().tobytes()).hexdigest()
    if receipt_path.exists():
        old=json.loads(receipt_path.read_text())
        if old.get('point_hash')==point_hash and old['status'] in ['COMPLETE','NO_TILE'] and data_path.exists():return old
    url='https://planetarycomputer.microsoft.com/api/stac/v1/collections/cop-dem-glo-90/items/'+tile
    receipt=dict(tile=tile,item_url=url,points=len(q),point_hash=point_hash,lock_sha256=lock_digest)
    try:
        cached=W/'dem_samples'/(tile+'.json')
        old=json.loads(cached.read_text()) if cached.exists() else {}
        if old.get('status')=='NO_TILE':
            q['elevation_m']=np.nan;q.to_parquet(data_path,index=False)
            receipt.update(status='NO_TILE',valid_points=0,metadata_source='prior_cell_source_404')
        else:
            if old.get('status')=='COMPLETE' and old.get('asset_url'):
                asset=old['asset_url'];date=old.get('item_datetime');metadata_source='prior_cell_STAC_receipt'
            else:
                response=requests.get(url,timeout=45)
                if response.status_code==404:
                    q['elevation_m']=np.nan;q.to_parquet(data_path,index=False)
                    receipt.update(status='NO_TILE',valid_points=0,metadata_source='new_STAC_404')
                    receipt_path.write_text(json.dumps(receipt));return receipt
                response.raise_for_status();item=response.json();asset=item['assets']['data']['href'];date=item.get('properties',{}).get('datetime');metadata_source='new_STAC'
            with rasterio.open(planetary_computer.sign_url(asset)) as ds:
                assert ds.crs.to_epsg()==4326
                a=ds.read(1,masked=True).astype('float32')
                rr,cc=rasterio.transform.rowcol(ds.transform,q.longitude.to_numpy(),q.latitude.to_numpy());rr=np.asarray(rr);cc=np.asarray(cc)
                for k in range(min(5,len(q))):assert (rr[k],cc[k])==ds.index(float(q.longitude.iloc[k]),float(q.latitude.iloc[k]))
                ok=(rr>=0)&(cc>=0)&(rr<a.shape[0])&(cc<a.shape[1]);v=np.full(len(q),np.nan)
                v[ok]=a[rr[ok],cc[ok]].filled(np.nan)
                q['elevation_m']=v;q.to_parquet(data_path,index=False)
                receipt.update(status='COMPLETE',valid_points=int(np.isfinite(v).sum()),outside_raster=int((~ok).sum()),asset_url=asset,item_datetime=date,metadata_source=metadata_source,raster_shape=list(ds.shape))
        receipt['data_sha256']=hashlib.sha256(data_path.read_bytes()).hexdigest()
    except Exception as exc:receipt.update(status='ERROR',error_type=type(exc).__name__)
    receipt_path.write_text(json.dumps(receipt));return receipt
counts={};valid=0
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    for i,result in enumerate(pool.map(acquire,jobs),1):
        counts[result['status']]=counts.get(result['status'],0)+1;valid+=result.get('valid_points',0)
        if i%10==0:print('route terrain tiles',i,'/',len(jobs),counts,flush=True)
summary=dict(status='ACQUISITION_PASS_FINISHED_NOT_FEATURE_VALIDATION',requested_tiles=len(jobs),all_design_tiles=len(tiles),counts=counts,valid_points=valid,partial_limit=args.limit,limits='DSM route samples only, not terrain resistance or model results. NO_TILE and failed retrieval remain missing; cache contains unsiged asset URLs, no access tokens.')
(O/('route_terrain_pilot_acquisition.json' if args.limit else 'route_terrain_acquisition.json')).write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from sklearn.neighbors import BallTree
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'transport_forecast_extension_lock.json',O/'port_connectivity_exact_matches.csv',W/'unctad_port_connectivity_2006_2021.parquet',
       O/'port_buffer_landmass_audit.csv',W/'ports_2019.geojson',W/'cell_geography.parquet']
lock=json.loads(files[0].read_text());m=pd.read_csv(files[1],keep_default_na=False)
series=pd.read_parquet(files[2]);geo=pd.read_csv(files[3],dtype={'single_landmass_id':str})
geo=geo[geo.status.eq('SINGLE_LANDMASS')]
m=m.merge(geo[['port_id','single_landmass_id']],on='port_id',validate='one_to_one')
features=json.loads(files[4].read_text())['features'];coords={f['properties']['OBJECTID']:f['geometry']['coordinates'][:2] for f in features}
m['longitude']=m.port_id.map(lambda pid:coords[pid][0]);m['latitude']=m.port_id.map(lambda pid:coords[pid][1])
cells=pd.read_parquet(files[5]);assert not cells.cell.duplicated().any()
outputs=[];receipts=[];portrows=[]
for period,start,end in [(0,2011,2014),(1,2015,2018)]:
    values=series[series.year.between(start,end)]
    index=values.groupby('series_code').agg(mean_connectivity=('value','mean'),quarters=('value','count'),max_year=('year','max')).reset_index()
    index=index[index.quarters.ge(13)]
    ports=m.merge(index,on='series_code',validate='one_to_one');ports['period']=period
    assert ports.mean_connectivity.ge(0).all() and ports.max_year.le(end).all()
    portrows.append(ports)
    result=cells[['cell','landmass_id','landmass_area_km2']].copy();result['period']=period
    result['port_id']=np.nan;result['log_nearest_indexed_port_km']=0.;result['indexed_port_support_unknown']=1
    result['log_nearest_port_connectivity']=0.;result['index_window_start']=start;result['index_window_end']=end
    for lid,cc in cells.dropna(subset=['landmass_id']).groupby('landmass_id'):
        pp=ports[ports.single_landmass_id.eq(lid)].sort_values('port_id').reset_index(drop=True)
        if pp.empty:continue
        tree=BallTree(np.radians(pp[['latitude','longitude']].values),metric='haversine')
        distance,ix=tree.query(np.radians(cc[['latitude','longitude']].values),k=1)
        selected=pp.iloc[ix.ravel()]
        result.loc[cc.index,'port_id']=selected.port_id.values
        result.loc[cc.index,'log_nearest_indexed_port_km']=np.log1p(distance.ravel()*6371.0088)
        result.loc[cc.index,'log_nearest_port_connectivity']=np.log1p(selected.mean_connectivity.values)
        result.loc[cc.index,'indexed_port_support_unknown']=0
    assert result[lock['features']].notna().all().all()
    assert result.loc[result.indexed_port_support_unknown.eq(1),'port_id'].isna().all()
    outputs.append(result)
    receipts.append(dict(period=period,years=[start,end],eligible_ports=len(ports),cells=len(result),
                         supported_cells=int(result.indexed_port_support_unknown.eq(0).sum()),
                         supported_landmasses=result.loc[result.indexed_port_support_unknown.eq(0),'landmass_id'].nunique()))
out=pd.concat(outputs,ignore_index=True);assert not out.duplicated(['cell','period']).any()
out.to_parquet(W/'transport_forecast_features.parquet',index=False)
pd.concat(portrows).to_parquet(W/'transport_forecast_port_windows.parquet',index=False)
audit=dict(status='FEATURE_JOIN_COMPLETE_NOT_MODEL_RESULT',periods=receipts,feature_columns=lock['features'],
           limits=lock['limits'],input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
(O/'transport_forecast_feature_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}))

from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
features=pd.read_parquet(W/'transport_forecast_features.parquet')
ports=pd.read_parquet(W/'transport_forecast_port_windows.parquet')
source=pd.read_parquet(W/'unctad_port_connectivity_2006_2021.parquet')
cells=pd.read_parquet(W/'cell_geography.parquet').set_index('cell')
assert len(features)==2*len(cells) and not features.duplicated(['cell','period']).any()
checked=0
for period,start,end in [(0,2011,2014),(1,2015,2018)]:
    f=features[features.period.eq(period)];p=ports[ports.period.eq(period)]
    values=source[source.year.between(start,end)].groupby('series_code').value.agg(['mean','count'])
    np.testing.assert_allclose(p.mean_connectivity,values.loc[p.series_code,'mean'])
    np.testing.assert_array_equal(p.quarters,values.loc[p.series_code,'count'])
    assert p.quarters.ge(13).all()
    known=f[f.indexed_port_support_unknown.eq(0)]
    joined=known.merge(p[['port_id','single_landmass_id','mean_connectivity']],on='port_id',validate='many_to_one')
    assert joined.landmass_id.eq(joined.single_landmass_id).all()
    np.testing.assert_allclose(joined.log_nearest_port_connectivity,np.log1p(joined.mean_connectivity))
    for lid,g in known.groupby('landmass_id'):
        pp=p[p.single_landmass_id.eq(lid)]
        cc=cells.loc[g.cell]
        lat1=np.radians(cc.latitude.values)[:,None];lon1=np.radians(cc.longitude.values)[:,None]
        lat2=np.radians(pp.latitude.values)[None,:];lon2=np.radians(pp.longitude.values)[None,:]
        aa=np.sin((lat1-lat2)/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin((lon1-lon2)/2)**2
        km=6371.0088*2*np.arcsin(np.sqrt(np.clip(aa,0,1)))
        np.testing.assert_allclose(np.expm1(g.log_nearest_indexed_port_km),km.min(axis=1),rtol=1e-9,atol=1e-8)
    unknown=f[f.indexed_port_support_unknown.eq(1)]
    assert unknown.port_id.isna().all() and unknown.log_nearest_indexed_port_km.eq(0).all() and unknown.log_nearest_port_connectivity.eq(0).all()
    assert not set(unknown.landmass_id.dropna())&set(p.single_landmass_id)
    checked+=len(f)
v=dict(status='PASS',feature_rows_checked=checked,checks=['baseline-only quarterly means and support','same-landmass port assignment','nearest distance independently recomputed','missing support retained with explicit flag'],
       limits='Does not establish historical port coordinates, crosswalk identity or forecast improvement.')
(O/'transport_forecast_feature_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8');print(json.dumps(v))

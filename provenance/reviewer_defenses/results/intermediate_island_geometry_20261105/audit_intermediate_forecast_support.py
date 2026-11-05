from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';F=R/'work/temporal/forecast';O=R/'outputs'
source=W/'intermediate_island_features.parquet'
f=pd.read_parquet(source)
assert not f.duplicated(['speciesKey','cell','period']).any()
assert (f.intermediate_islands.isna()==f.intermediate_identity_unknown.eq(1)).all()
np.testing.assert_allclose(f.log_intermediate_islands,np.log1p(f.intermediate_islands.fillna(0)))
old=pd.read_parquet(W/'sea_route_features.parquet')
q=f.merge(old[['speciesKey','cell','period','route_endpoint_unknown','sampled_marine_segments']],on=['speciesKey','cell','period'],validate='one_to_one')
assert len(q)==len(old)==len(f)
assert q.loc[q.intermediate_islands.gt(0),'sampled_marine_segments'].ge(2).all()
assert q.loc[q.route_endpoint_unknown.eq(1),'intermediate_identity_unknown'].eq(1).all()
p=pd.read_parquet(F/'species_cell_period_panel.parquet',columns=['speciesKey','species','cell','period','recording_probability10'])
d=pd.read_csv(O/'forecast_all_225_final_coverage.csv')[['species','taxon_key','final_disposition']]
assert len(d)==225 and d.taxon_key.is_unique
summary=[]
for k,label in [(0,'train'),(1,'test')]:
    a=p[p.period.eq(k)&p.recording_probability10.notna()]
    b=p[p.period.eq(k+1)&p.recording_probability10.notna()][['speciesKey','cell']]
    a=a.merge(b,on=['speciesKey','cell'],validate='one_to_one')
    t=a.merge(f,on=['speciesKey','cell','period'],how='left',validate='one_to_one')
    assert len(a)==len(t) and t.intermediate_identity_unknown.notna().all()
    for suffix,mask in [('evaluated',np.ones(len(t),bool)),('intermediate_positive',t.intermediate_islands.gt(0)),('geometry_unknown',t.intermediate_identity_unknown.eq(1))]:
        z=t[mask];counts=z.groupby('speciesKey').size()
        d[label+'_'+suffix+'_rows']=d.taxon_key.map(counts).fillna(0).astype(int)
        assert d[label+'_'+suffix+'_rows'].sum()==len(z)
        summary.append(dict(period=label,category=suffix,rows=len(z),species=z.speciesKey.nunique(),cells=z.cell.nunique(),routes=z.route_key.nunique()))
    if k==1:
        ref=pd.read_parquet(F/'predictions_sea_temporal.parquet',columns=['speciesKey','cell'])
        np.testing.assert_array_equal(a[['speciesKey','cell']],ref[['speciesKey','cell']])
d.to_csv(O/'intermediate_island_all225_support.csv',index=False)
pd.DataFrame(summary).to_csv(O/'intermediate_island_forecast_support.csv',index=False)
audit=dict(status='PASS_FEATURE_SUPPORT_NOT_MODEL_RESULT',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),summary=summary,checks=['225-candidate ledger retained','unique feature keys','unknown preserved','intermediate-positive implies multiple sea runs','all training/test eligible rows have a feature match','test keys exactly match B5 temporal predictions'],limits='No response association tested; positive means sampled geometric intermediate island, not occupancy or relay.')
(O/'intermediate_island_support_audit.json').write_text(json.dumps(audit,indent=2))
print(pd.DataFrame(summary).to_string(index=False))

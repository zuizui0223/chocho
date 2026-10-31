from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];S=R/'work/temporal/forecast';O=R/'outputs'
cols=['speciesKey','cell','period','recording_probability10','spatial_block','climate_id','latitude','longitude']
p=pd.read_parquet(S/'species_cell_period_panel.parquet',columns=cols)
a=p[p.period.eq(0)&p.recording_probability10.notna()]
b=p[p.period.eq(1)&p.recording_probability10.notna()][['speciesKey','cell']]
train=a.merge(b,on=['speciesKey','cell'],validate='one_to_one')
train['fold']=train.spatial_block.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)
test=pd.read_parquet(S/'predictions_spatial.parquet',columns=['speciesKey','cell','spatial_block','climate_id','spatial_fold','distance_to_training_cell_km'])
coords=p[p.period.eq(1)][['cell','latitude','longitude']].drop_duplicates('cell')
test=test.merge(coords,on='cell',validate='many_to_one')
assert (test.spatial_fold==test.spatial_block.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)).all()
checks=[]
for fold in range(5):
    tr=train[train.fold!=fold];te=test[test.spatial_fold==fold]
    assert set(tr.cell).isdisjoint(te.cell)
    assert set(tr.climate_id).isdisjoint(te.climate_id)
    assert set(tr.spatial_block).isdisjoint(te.spatial_block)
    trcoords=tr[['cell','latitude','longitude']].drop_duplicates('cell')
    sample=te.drop_duplicates('cell').sample(min(20,te.cell.nunique()),random_state=20261029)
    lat1=np.radians(sample.latitude.to_numpy())[:,None];lon1=np.radians(sample.longitude.to_numpy())[:,None]
    lat2=np.radians(trcoords.latitude.to_numpy())[None,:];lon2=np.radians(trcoords.longitude.to_numpy())[None,:]
    # Direct all-pairs haversine minimum, independent of the model's BallTree.
    h=np.sin((lat2-lat1)/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
    dist=2*6371.0088*np.arcsin(np.sqrt(np.clip(h,0,1)))
    err=float(abs(dist.min(axis=1)-sample.distance_to_training_cell_km.to_numpy()).max())
    assert err<1e-7
    checks.append({'fold':fold,'train_pairs':len(tr),'test_pairs':len(te),'buffer100_pairs':int(te.distance_to_training_cell_km.ge(100).sum()),'direct_distance_checks':len(sample),'max_distance_error_km':err})
receipt={'status':'PASS','fold_checks':checks,'buffer100_pairs':int(test.distance_to_training_cell_km.ge(100).sum()),'buffer100_species':int(test.loc[test.distance_to_training_cell_km.ge(100),'speciesKey'].nunique()),'buffer100_cells':int(test.loc[test.distance_to_training_cell_km.ge(100),'cell'].nunique()),'interpretation':'Validates spatial partitions and distances only, not predictive performance or generality'}
(O/'forecast_spatial_split_validation.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))

"""Compare coastal port-ring samples to same-island shoreline controls.
Refuse ecological summaries until every expected control sample is resolved.
"""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
design=W/'coastal_control_validated_design.parquet';digest=hashlib.sha256(design.read_bytes()).hexdigest()
p=pd.read_parquet(design);samples=[];missing=[];errors=[]
for f in (W/'coastal_control_vegetation_tiles').glob('*.json'):
    a=json.loads(f.read_text())
    if a.get('design_sha256')!=digest:continue
    if a['status']=='COMPLETE':
        v=np.array([s['landcover'] for s in a['samples']],dtype=np.uint8)
        assert hashlib.sha256(v.tobytes()).hexdigest()==a['values_sha256']
        samples.extend(a['samples'])
    elif a['status']=='SOURCE_NO_PUBLISHED_TILE':missing.extend(a['missing_sample_ids'])
    else:errors.append(f.name)
q=pd.DataFrame(samples,columns=['sample_id','landcover'])
assert not q.sample_id.duplicated().any()
assert not set(q.sample_id).intersection(missing)
expected=set(p.sample_id);resolved=set(q.sample_id).union(missing)
assert resolved<=expected
receipt=dict(status='COMPLETE_WITH_SOURCE_GAPS' if missing else 'COMPLETE',
             design_samples=len(p),retrieved_samples=len(q),missing_samples=len(missing),
             pending_samples=len(expected-resolved),error_tiles=errors,design_sha256=digest)
if errors or expected!=resolved:
    receipt['status']='ACQUISITION_INCOMPLETE_NO_COMPARISON'
    (O/'coastal_control_vegetation_status.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt));raise SystemExit(2)
far=p.merge(q,on='sample_id',how='left',validate='one_to_one');far['terrestrial_valid']=far.landcover.notna()&~far.landcover.isin([0,80])
far.to_parquet(W/'coastal_control_vegetation_joined.parquet',index=False)
near=pd.read_parquet(W/'port_vegetation_samples_joined.parquet')
probe=pd.read_parquet(W/'port_island_coastal_probe.parquet')
near=near.merge(probe,on='sample_id',how='inner',validate='one_to_one')
near=near[near.terrestrial_valid & near.coastal_probe_1km & near.radius_km.isin([1,5])]
far.landmass_id=far.landmass_id.astype(str);near.landmass_id=near.landmass_id.astype(str)
far_groups={key:g[g.terrestrial_valid] for key,g in far.groupby(['port_id','landmass_id'])}
features={10:'tree',20:'shrub',30:'grass',40:'crop',50:'built'};rows=[]
for (port,lid,radius),g in near.groupby(['port_id','landmass_id','radius_km']):
    f=far_groups.get((port,lid))
    if f is None or len(f)<3 or len(g)<3:continue
    row=dict(port_id=port,landmass_id=lid,near_radius_km=radius,near_points=len(g),far_points=len(f))
    for code,name in features.items():
        row[name+'_near']=g.landcover.eq(code).mean();row[name+'_far']=f.landcover.eq(code).mean()
        row[name+'_difference']=row[name+'_near']-row[name+'_far']
    rows.append(row)
paired=pd.DataFrame(rows);paired.to_csv(O/'coastal_control_vegetation_pairs.csv',index=False)
summary=[]
if len(paired):
    for radius,g in paired.groupby('near_radius_km'):
        means=g.groupby('landmass_id')[[name+'_difference' for name in features.values()]].mean()
        rng=np.random.default_rng(20261030);ix=rng.integers(len(means),size=(2000,len(means)))
        for name in features.values():
            values=means[name+'_difference'].to_numpy()
            lo,hi=np.quantile(values[ix].mean(axis=1),[.025,.975])
            summary.append(dict(near_radius_km=radius,vegetation=name,mean_difference=values.mean(),
                ci_low=lo if len(means)>=20 else np.nan,ci_high=hi if len(means)>=20 else np.nan,
                islands=len(means),ports=g.port_id.nunique(),pairs=len(g)))
pd.DataFrame(summary).to_csv(O/'coastal_control_vegetation_summary.csv',index=False)
receipt.update(pairs=len(paired),limits=['coastal probe near ports and250m shore-offset controls are not exact distance matches',
    'point frequencies, not buffer-area cover','static2021 vegetation cannot date introductions',
    'only islands supporting same-island controls; no global representativeness claim',
    'island-weighted bootstrap does not remove environmental confounding'])
(O/'coastal_control_vegetation_status.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt))

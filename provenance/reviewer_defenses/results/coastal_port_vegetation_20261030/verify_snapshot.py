from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
for name,digest in json.loads((R/'snapshot_manifest.json').read_text()).items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
p=pd.read_csv(R/'coastal_control_vegetation_pairs.csv',dtype={'landmass_id':str})
s=pd.read_csv(R/'coastal_control_vegetation_summary.csv')
assert not p.duplicated(['port_id','landmass_id','near_radius_km']).any()
assert (p[['near_points','far_points']]>=3).all().all()
for radius,g in p.groupby('near_radius_km'):
    means=g.groupby('landmass_id')[[name+'_difference' for name in ['tree','shrub','grass','crop','built']]].mean()
    rng=np.random.default_rng(20261030);ix=rng.integers(len(means),size=(2000,len(means)))
    for row in s[s.near_radius_km.eq(radius)].itertuples():
        name=row.vegetation
        assert np.allclose(g[name+'_near']-g[name+'_far'],g[name+'_difference'])
        values=means[name+'_difference'].to_numpy();lo,hi=np.quantile(values[ix].mean(axis=1),[.025,.975])
        assert np.isclose(values.mean(),row.mean_difference)
        assert np.isclose(lo,row.ci_low) and np.isclose(hi,row.ci_high)
        assert (len(means),g.port_id.nunique(),len(g))==(row.islands,row.ports,row.pairs)
print('PASS: snapshot hashes, 731 pair differences, all 10 means/intervals/support counts; not independent raster QA')

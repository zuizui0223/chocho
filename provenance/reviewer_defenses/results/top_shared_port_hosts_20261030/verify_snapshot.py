from pathlib import Path
import json,hashlib,pandas as pd,numpy as np
r=Path(__file__).resolve().parent
for n,h in json.loads((r/'snapshot_manifest.json').read_text()).items():
 assert hashlib.sha256((r/n).read_bytes()).hexdigest()==h,n
c=pd.read_csv(r/'top_shared_host_region_correlations.csv');j=pd.read_csv(r/'top_shared_host_region_join.csv');s=pd.read_csv(r/'top_shared_host_region_summary.csv');lock=json.loads((r/'top_shared_host_region_analysis_lock.json').read_text())
assert set(c.host_name)==set(lock['hosts'])
assert not j.duplicated(['host_name','species','region']).any()
assert j.host_ports.le(j.complete_ports).all() and j.cells.ge(5).all()
assert np.allclose(j.port_fraction,j.host_ports/j.complete_ports)
for z in c.itertuples():
 q=j[j.host_name.eq(z.host_name)&j.species.eq(z.species)]
 assert len(q)==z.regions and len(q)>=5
 assert np.isclose(np.corrcoef(q.port_fraction.rank(),q[z.outcome].rank())[0,1],z.rho)
for z in s.itertuples():
 q=c[c.host_name.eq(z.host_name)&c.outcome.eq(z.outcome)]
 assert len(q)==z.species and np.isclose(q.rho.mean(),z.mean_rho)
 assert q.rho.gt(0).sum()==z.positive_species
print('PASS all72correlations, regional fractions, support, summaries and snapshot hashes')

p=pd.read_csv(r/'top_shared_host_link_provenance.csv');s=pd.read_csv(r/'top_shared_host_source_sensitivity.csv')
x=c.merge(p,on=['host_name','species'],validate='many_to_one');assert len(x)==len(c)
for z in s.itertuples():
 q=x[x.host_name.eq(z.host_name)&x.outcome.eq(z.outcome)]
 if z.source=='frozen_HOSTS_only':q=q[q.in_frozen_HOSTS]
 assert len(q)==z.species and np.isclose(q.rho.mean(),z.mean_rho)
 assert q.rho.gt(0).sum()==z.positive_species
print('PASS all36source sensitivity summaries')

from pathlib import Path
import pandas as pd,json,hashlib
r=Path(__file__).resolve().parent
for n,h in json.loads((r/'snapshot_manifest.json').read_text()).items():
 assert hashlib.sha256((r/n).read_bytes()).hexdigest()==h,n
p=pd.read_csv(r/'port_alien_host_broad_all_ports.csv',keep_default_na=False,na_values=['']);c=pd.read_csv(r/'port_alien_host_broad_country_coverage.csv',keep_default_na=False,na_values=['']);s=pd.read_csv(r/'port_alien_host_broad_species.csv');a=json.loads((r/'port_alien_host_broad_validation.json').read_text())
assert len(p)==a['registered_ports'] and p.port_id.is_unique
q=p[p.query_complete];col='contemporary_not_native_in_catalogue'
assert len(q)==a['complete_ports'] and q[col].gt(0).sum()==a['ports_with_alien_candidates']
assert p.loc[~p.query_complete,col].isna().all()
assert s.ports.sum()==q[col].sum()==a['port_alien_host_pairs']
assert len(s)==a['alien_candidate_species']
for z in c.itertuples():
 x=q[q.country.eq(z.country)]
 assert len(x)==z.complete_ports
 assert x[col].gt(0).sum()==z.ports_with_alien_candidates
 assert (x[col]+x.native_in_catalogue).gt(0).sum()==z.ports_with_classified_hosts
print('PASS snapshot hashes; complete and pending scope; all country and species aggregate counts')

from pathlib import Path
import json,hashlib,pandas as pd
r=Path(__file__).resolve().parent
for n,h in json.loads((r/'snapshot_manifest.json').read_text()).items():
 assert hashlib.sha256((r/n).read_bytes()).hexdigest()==h,n
x=pd.read_csv(r/'port_alien_host_butterfly_corecording.csv');s=pd.read_csv(r/'port_alien_host_butterfly_all225_summary.csv')
a=json.loads((r/'port_alien_host_butterfly_corecording_audit.json').read_text())
assert len(s)==225 and s.species.is_unique and not x.duplicated(['species','port_id']).any()
assert len(x)==225*a['complete_plant_ports']
assert (x.both_recorded==(x.recorded_alien_candidate_host_species.gt(0)&x.butterfly_records.gt(0))).all()
for z in s.itertuples():
 q=x[x.species.eq(z.species)]
 assert len(q)==a['complete_plant_ports']
 assert q.both_recorded.sum()==z.ports_with_both
 assert q.butterfly_records.gt(0).sum()==z.ports_with_butterfly_records
 assert q.recorded_alien_candidate_host_species.gt(0).sum()==z.ports_with_alien_candidate_hosts
assert x.both_recorded.sum()==a['species_port_pairs_with_both']
assert s.ports_with_both.gt(0).sum()==a['species_with_both_at_same_port']
print('PASS snapshot hashes, all225 coverage, every aggregate and co-recording flag')

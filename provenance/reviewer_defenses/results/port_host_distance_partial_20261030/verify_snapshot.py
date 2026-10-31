from pathlib import Path
import json,hashlib,csv,math
r=Path(__file__).resolve().parent
for name,h in json.loads((r/'snapshot_manifest.json').read_text()).items():
 assert hashlib.sha256((r/name).read_bytes()).hexdigest()==h,name
a=json.loads((r/'port_host_distance_join_audit.json').read_text())
s=list(csv.DictReader((r/'port_host_distance_bin_summary.csv').open()))
assert sum(int(x['plant_records']) for x in s)==a['unique_records']
assert sum(int(x['nonnative_candidate_host_records']) for x in s)==a['host_origin_counts']['NONNATIVE_CANDIDATE']
for x in s:
 assert math.isclose(float(x['nonnative_candidate_fraction']),int(x['nonnative_candidate_host_records'])/int(x['classified_host_records']))
v=json.loads((r/'port_host_distance_validation.json').read_text())
assert v['snapshot_unique_records']==a['unique_records']
print('PASS frozen hashes, bin fractions, counts, validation snapshot identity; source distance checks reported separately')

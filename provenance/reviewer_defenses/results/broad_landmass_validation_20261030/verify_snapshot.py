"""Verify public bytes and aggregate consistency, not a source-data rerun."""
from pathlib import Path
import csv,json,hashlib,math
r=Path(__file__).resolve().parent
m=json.loads((r/'snapshot_manifest.json').read_text())
for name,expected in m.items():
 assert hashlib.sha256((r/name).read_bytes()).hexdigest()==expected,name
scores=list(csv.DictReader((r/'barrier_recomputed_scores.csv').open()))
receipts={mode:json.loads((r/f'forecast_barrier_{mode}_results.json').read_text()) for mode in ['temporal','spatial']}
j=receipts['temporal']
assert j['status']=='COMPLETE' and j['test_species']==155 and j['test_rows']==208482
for row in scores:
 if row['group']=='all':
  j=receipts[row['mode']]
  assert abs(float(row['score'])-j['results'][row['model']]['species_equal_'+row['metric']])<1e-12
ci=list(csv.DictReader((r/'barrier_paired_uncertainty.csv').open()))
for row in ci:
 values=[float(row[k]) for k in ['difference','ci_low','ci_high']]
 assert all(math.isfinite(v) for v in values)
 assert values[1]<=values[2]
 if row['group']=='all':
  j=receipts[row['mode']]
  metric='species_equal_'+row['metric']
  expected=j['results'][row['added']][metric]-j['results'][row['reference']][metric]
  assert abs(values[0]-expected)<1e-12
assert json.loads((r/'barrier_history_feature_validation.json').read_text())['sample_rows']==40
print('PASS',len(m),'hashes;',len(ci),'paired summaries; temporal/spatial score consistency. Not a scientific rerun.')

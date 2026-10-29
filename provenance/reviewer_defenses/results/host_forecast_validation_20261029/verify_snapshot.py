"""Verify published snapshot consistency; not an ecological model rerun."""
from pathlib import Path
import csv,hashlib,json,math
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'snapshot_manifest.json').read_text())
for item in manifest['files']:
    assert hashlib.sha256((ROOT/item['file']).read_bytes()).hexdigest()==item['sha256'],item['file']
def read_csv(name):
    with (ROOT/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
metrics={(r['mode'],r['model']):r for r in read_csv('forecast_validation_all_metrics.csv')}
for mode in manifest['completed_modes']:
    receipt=json.loads((ROOT/f'forecast_{mode}_results.json').read_text())
    assert receipt['status']=='COMPLETE'
    for model,values in receipt['results'].items():
        for metric in ['mse','mae']:
            key='species_equal_'+metric
            assert math.isclose(float(metrics[mode,model][key]),values[key],rel_tol=1e-12,abs_tol=1e-14)
contrasts=read_csv('forecast_paired_comparisons.csv')
for row in contrasts:
    mode=row['mode'];key='species_equal_'+row['metric']
    expected=float(metrics[mode,row['added']][key])-float(metrics[mode,row['reference']][key])
    assert math.isclose(expected,float(row['difference']),rel_tol=1e-10,abs_tol=1e-12)
    assert float(row['ci_low'])<=float(row['ci_high'])
census=read_csv('forecast_all_225_final_coverage.csv')
assert len(census)==225 and len({r['species'] for r in census})==225
for mode in manifest['completed_modes']:
    if mode=='spatial_buffer100':continue
    assert sum(int(r[mode+'_test_pairs']) for r in census)==int(metrics[mode,'M5_host_composition']['rows'])
print(f"PASS: {len(manifest['files'])} hashes; {len(metrics)} model summaries; {len(contrasts)} contrasts; 225-species census")

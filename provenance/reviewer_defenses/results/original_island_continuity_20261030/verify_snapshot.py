from pathlib import Path
import csv,hashlib,json,math,statistics
R=Path(__file__).resolve().parent
manifest=json.loads((R/'snapshot_manifest.json').read_text())
for name,digest in manifest.items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
def read(name):
    with (R/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
raw=read('global_original_continuity_species_counts.csv')
summary=read('global_original_continuity_distance_summary.csv')
for row in summary:
    q=[r for r in raw if r['cohort']==row['cohort'] and r['distance_band']==row['distance_band']]
    count=sum(int(r['events']) for r in q);total=sum(int(r['total_events']) for r in q)
    assert count==int(row['events']) and total==int(row['total_events'])
    assert len(q)==int(row['species'])
    assert math.isclose(count/total,float(row['event_weighted_share']),abs_tol=1e-12)
    assert math.isclose(statistics.mean(int(r['events'])/int(r['total_events']) for r in q),float(row['species_equal_share']),abs_tol=1e-12)
for cohort in {r['cohort'] for r in summary}:
    q=[r for r in summary if r['cohort']==cohort]
    assert sum(int(r['events']) for r in q)==int(q[0]['total_events'])
print('PASS: hashes and all 12 aggregate rows recomputed from species counts; raw geography not independently rerun')

from pathlib import Path
import csv,json,hashlib,math,statistics
R=Path(__file__).resolve().parent
manifest=json.loads((R/'snapshot_manifest.json').read_text(encoding='utf-8'))
for name,digest in manifest.items():
 assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
read=lambda name:list(csv.DictReader((R/name).open(encoding='utf-8',newline='')))
raw=read('island_within_species_landmass_correlations.csv');summaries=read('island_correlation_summary.csv')
for summary in summaries:
 rows=[x for x in raw if x['feature']==summary['feature']];bysp={}
 for row in rows:bysp.setdefault(row['speciesKey'],[]).append(float(row['rho']))
 species_means=[statistics.mean(v) for v in bysp.values()]
 assert len(bysp)==int(summary['species'])
 assert len(rows)==int(summary['species_island_groups'])
 assert len({x['landmass_id'] for x in rows})==int(summary['landmasses'])
 assert math.isclose(statistics.mean(species_means),float(summary['mean_within_island_rho_species_equal']),abs_tol=1e-12)
 assert math.isclose(statistics.median(species_means),float(summary['median_species_rho']),abs_tol=1e-12)
coefficients=read('port_host_adjusted_associations.csv');sensitivity=read('port_host_adjusted_associations_interval_sensitivity.csv')
assert len(coefficients)==len(sensitivity)==51
for row,alternative in zip(coefficients,sensitivity):
 assert row['subset']==alternative['subset'] and row['feature']==alternative['feature']
 beta=float(row['coefficient']);se=float(row['cluster_se'])
 assert math.isclose(float(row['ci_low']),beta-1.96*se,abs_tol=1e-12)
 assert math.isclose(float(row['ci_high']),beta+1.96*se,abs_tol=1e-12)
 assert math.isclose(float(alternative['coefficient']),beta,abs_tol=1e-12)
 assert float(alternative['cluster_t_low'])<=float(row['ci_low']) and float(alternative['cluster_t_high'])>=float(row['ci_high'])
audit=json.loads((R/'port_host_association_audit.json').read_text(encoding='utf-8'))
for fit,expected in zip(audit['fits'],[(49509,93,22,975,18),(5903,93,20,116,15),(22082,82,10,436,9)]):
 assert tuple(fit[x] for x in ['rows','species','landmasses','cells','spatial_blocks'])==expected
 assert fit['rank']==17 and fit['dropped_invariant_features']==['introduced_centered']
 assert fit['coefficient_lstsq_max_difference']<1e-9
 assert fit['negative_covariance_diagonal']==0
print('PASS:',len(manifest),'hashes;',len(summaries),'correlation summaries;',len(coefficients),'coefficient/interval pairs; 3 support audits')

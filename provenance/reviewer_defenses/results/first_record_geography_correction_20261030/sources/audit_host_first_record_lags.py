from pathlib import Path
import sys,json,re
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs';S=R/'work/remote_sensing/host_source/extracted'
sys.path.insert(0,str(R/'work/audit/uk-snapshot/scripts'))
from analyze_bce_augmented_plant_network_climate_null_20261019 import source_augmented_model
model,regions,climate,provenance=source_augmented_model(S)
targets=set(pd.read_csv(O/'forecast_candidate_disposition.csv').species)
names={};conflicts=[]
raw=pd.read_csv(S/'inputs/original/insect_host_accepted.csv',dtype=str)
for h,n in raw[['accepted_plant_name_id','accepted_name']].drop_duplicates().itertuples(index=False,name=None):
 if h in names and names[h]!=n:conflicts.append(h)
 names[h]=n
bce=json.loads((S/'inputs/bce/accepted_host_edges.json').read_text())
for rec in bce['per_butterfly']:
 for item in rec['details']:
  if item.get('match_status')!='ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY':continue
  h=str(item['accepted_plant_name_id']);n=item['accepted_name']
  if h in names and names[h]!=n:conflicts.append(h)
  names[h]=n
for h in conflicts:names.pop(h,None)
edges=[];unknown=[]
for sp,hosts in zip(model['species'],model['species_hosts']):
 if sp not in targets:continue
 for host in hosts:
  if str(host) in names:edges.append({'species':sp,'host_name':names[str(host)]})
  else:unknown.append({'species':sp,'host_id':str(host)})
edges=pd.DataFrame(edges).drop_duplicates()
assert len(targets)==225
d=pd.read_excel(W/'GlobalAlienSpeciesFirstRecordDatabase_v3.1_freedata.xlsx',sheet_name='FirstRecords',dtype=str)
def interval(value):
 s=str(value).strip()
 if re.fullmatch(r'\d{4}',s):return int(s),int(s),'exact'
 m=re.fullmatch(r'(\d{4})\s*[-–]\s*(\d{4})',s)
 if m:
  lo,hi=map(int,m.groups())
  if lo<=hi:return lo,hi,'interval'
 if re.fullmatch(r'\d{3}0s',s):return int(s[:4]),int(s[:4])+9,'decade'
 return np.nan,np.nan,'unresolved'
ints=d.FirstRecord_orig.map(interval)
d[['lower','upper','date_type']]=pd.DataFrame(ints.tolist(),index=d.index)
but=d[d.TaxonName.isin(targets)].copy().rename(columns={'TaxonName':'species'})
plant=d[d.TaxonName.isin(set(edges.host_name))].copy().rename(columns={'TaxonName':'host_name'})
linked=but.merge(edges,on='species').merge(plant,on=['Region','host_name'],suffixes=('_butterfly','_plant'))
linked['lag_lower']=linked.lower_butterfly-linked.upper_plant
linked['lag_upper']=linked.upper_butterfly-linked.lower_plant
linked['ordering']=np.select([linked.lag_lower.isna()|linked.lag_upper.isna(),linked.lag_lower.gt(0),linked.lag_upper.lt(0)],['unresolved_date','plant_recorded_first','butterfly_recorded_first'],default='overlapping_intervals')
# Source-linked edge/date rows stay local: avoid republishing restricted host catalogue.
linked.to_parquet(W/'host_butterfly_first_record_pairs_source_unfiltered.parquet',index=False)
# The cited primary paper covers South African sub-Antarctic islands, not Canada.
# Keep source geography intact and exclude the conflicting join; do not guess a relabel.
conflict=linked.Region.eq('Canada') & (linked.Source_butterfly.eq('Greve et al. (2017)') | linked.Source_plant.eq('Greve et al. (2017)'))
quarantined_pairs=int(conflict.sum())
linked.loc[conflict].to_parquet(W/'host_butterfly_first_record_pairs_geography_quarantine.parquet',index=False)
linked=linked.loc[~conflict].copy()
linked.to_parquet(W/'host_butterfly_first_record_pairs.parquet',index=False)
by=linked.groupby(['species','Region']).agg(host_pairs=('host_name','size'),dated_pairs=('lag_lower','count'),earliest_possible_lag=('lag_lower','min'),latest_possible_lag=('lag_upper','max')).reset_index()
by.to_csv(O/'first_record_lag_butterfly_region_summary.csv',index=False)
coverage=[]
for sp in sorted(targets):
 b=but[but.species==sp];q=linked[linked.species==sp]
 coverage.append({'species':sp,'butterfly_first_record_rows':len(b),'host_names_in_catalogue':edges.loc[edges.species==sp,'host_name'].nunique(),'same_region_host_pairs':len(q),'dated_host_pairs':int(q.lag_lower.notna().sum())})
pd.DataFrame(coverage).to_csv(O/'first_record_lag_all225_disposition.csv',index=False)
receipt={'status':'DESCRIPTIVE_SOURCE_PAIRING_COMPLETE','candidate_species':225,'catalogue_named_edges':len(edges),'unmapped_host_ids':len(unknown),'name_conflicts':len(set(conflicts)),'butterfly_source_rows':len(but),'plant_source_rows_matching_host_names':len(plant),'paired_rows':len(linked),'paired_butterfly_species':linked.species.nunique(),'paired_regions':linked.Region.nunique(),'butterfly_region_groups':len(by),'date_order_counts':linked.ordering.value_counts().to_dict(),'limits':'First-record dates only, not independent introduction/arrival dates; unresolved dates excluded from ordering; exact source Region/name matching, incomplete catalogue and recording coverage. No causal or co-transport inference.'}
receipt.update(status='KNOWN_GEOGRAPHY_CONFLICT_SCREENED_OTHER_SOURCES_NOT_FULLY_VERIFIED',quarantined_region_source_pairs=quarantined_pairs,geography_evidence='https://journals.abcjournal.aosis.co.za/index.php/abc/article/view/2143/2182')
(O/'first_record_lag_audit.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2),flush=True)

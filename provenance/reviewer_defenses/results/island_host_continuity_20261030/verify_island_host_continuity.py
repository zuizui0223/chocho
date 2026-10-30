from pathlib import Path
import json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[2]; O=R/'outputs'; W=R/'work/barriers'
d=pd.read_csv(O/'island_host_continuity_events.csv',dtype={'landmass_id':str})
y=pd.read_parquet(W/'global_original_species_island_year.parquet');y.landmass_id=y.landmass_id.astype(str)
raw=pd.read_csv(O/'global_original_sequence_observation_opportunity.csv',dtype={'landmass_id':str})
raw=raw[raw.first_record_year.between(2011,2017)]
keys=['speciesKey','landmass_id']
assert set(map(tuple,d[keys].values))==set(map(tuple,raw[keys].values))
z=d[keys+['first_record_year']].merge(y[keys+['year']],on=keys,validate='one_to_many')
z['lag']=z.year-z.first_record_year;z=z[z.lag.between(1,5)]
counts=z.groupby(keys).year.nunique();year5=z[z.lag.eq(5)].groupby(keys).size()
ix=pd.MultiIndex.from_frame(d[keys])
np.testing.assert_array_equal(d.followup_record_years,counts.reindex(ix,fill_value=0))
np.testing.assert_array_equal(d.record_in_year5,year5.reindex(ix,fill_value=0).gt(0))
np.testing.assert_array_equal(d.any_followup_record,d.followup_record_years.gt(0))
p=pd.read_csv(O/'island_host_continuity_within_species.csv')
metrics=['any_followup_record','followup_record_years','record_in_year5']
for m in metrics:
    pivot=d.groupby(['species','host_group'])[m].mean().unstack()
    v=(pivot.ALIEN_CANDIDATE_RECORDED-pivot.NATIVE_ONLY_CLASSIFIED).dropna()
    assert set(v.index)==set(p.species)
    np.testing.assert_allclose(p[m+'_difference'],v.reindex(p.species))
c=pd.read_csv(O/'island_host_continuity_all225.csv');assert len(c)==225 and c.species.nunique()==225
assert c.index_events.sum()==len(d)
v={'status':'PASS','events_checked':len(d),'paired_species_checked':len(p),
   'minimum_events_per_group_in_paired_comparison':int(p[['alien_events','native_events']].min().min()),
   'paired_species_at_least3_events_each':int((p.alien_events.ge(3)&p.native_events.ge(3)).sum()),
   'checks':['exact event cohort equality','all followup counts recomputed via year join','all year5 flags','all paired differences','225 candidate accounting'],
   'not_verified':'Underlying occurrence truth, causal effect, transport pathway, establishment; host classification relies on source inventory audit.'}
(O/'island_host_continuity_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8');print(json.dumps(v))

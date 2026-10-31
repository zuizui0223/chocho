from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
prefixes={'Earlier, all':'island_host_family_null_earlier','New, all':'island_host_family_null',
          'Earlier, common59':'island_host_family_null_earlier_common','New, common59':'island_host_family_null_common'}
rows=[];audits={};draws={}
for label,prefix in prefixes.items():
    a=json.loads((O/(prefix+'_audit.json')).read_text());d=pd.read_csv(O/(prefix+'_draws.csv'))
    np.testing.assert_allclose(d.species_equal_fraction.mean(),a['null_mean_species_equal'])
    assert len(d)==199 and a['degree_checks_all_draws'] and a['family_degree_checks_all_draws']
    rows.append(dict(cohort=label,species=a.get('both_island_species',a.get('both_new_species')),
                     species_pairs=a.get('both_island_species_pairs',a.get('both_new_species_pairs')),
                     observed_species_equal=a['observed_species_equal_fraction'],null_species_equal=a['null_mean_species_equal'],
                     excess_species_equal=a['observed_species_equal_fraction']-a['null_mean_species_equal']))
    audits[label]=a;draws[label]=d
s=pd.DataFrame(rows);s.to_csv(O/'island_host_early_new_summary.csv',index=False)
old=pd.read_csv(O/'island_host_family_null_earlier_common_observed_pairs.csv',dtype={'island_a':str,'island_b':str})
new=pd.read_csv(O/'island_host_family_null_common_observed_pairs.csv',dtype={'island_a':str,'island_b':str})
assert set(old.species)==set(new.species) and old.species.nunique()==59
history=pd.read_csv(O/'global_original_species_island_history.csv',dtype={'landmass_id':str})
candidate=pd.read_csv(O/'forecast_candidate_disposition.csv')[['species','taxon_key']]
assert len(candidate)==225
history=history.merge(candidate,left_on='speciesKey',right_on='taxon_key',validate='many_to_one')
for df,mask in [(old,history.first_record_year.le(2010)),(new,history.first_record_year.between(2011,2022))]:
    valid=set(map(tuple,history.loc[mask,['species','landmass_id']].values))
    assert all((r.species,r.island_a) in valid and (r.species,r.island_b) in valid for r in df.itertuples())
candidate['common_period_comparison']=candidate.species.isin(old.species)
candidate.to_csv(O/'island_host_early_new_all225.csv',index=False)
a=draws['Earlier, common59'];b=draws['New, common59']
np.testing.assert_array_equal(a.draw,b.draw);np.testing.assert_array_equal(a.original_edge_retention,b.original_edge_retention)
delta=b.species_equal_fraction-a.species_equal_fraction
observed=audits['New, common59']['observed_species_equal_fraction']-audits['Earlier, common59']['observed_species_equal_fraction']
pd.DataFrame({'draw':a.draw,'new_minus_earlier':delta}).to_csv(O/'island_host_early_new_null_differences.csv',index=False)
v=dict(status='PASS_DESCRIPTIVE_COMPARISON_NOT_FORECAST',common_species=59,
       observed_new_minus_earlier=observed,null_difference_mean=float(delta.mean()),
       null_difference_quantiles=delta.quantile([.025,.5,.975]).to_dict(),
       enrichment_difference=observed-float(delta.mean()),draws_at_least_observed_difference=int(delta.ge(observed).sum()),
       checks=['both periods have identical59 species','all included species-island records satisfy period rules','225 candidate accounting',
               '199 draw summary recomputed','identical draw indices and original-edge-retention trace'],
       limits='Same deterministic seed/input edge order/algorithm used for paired draws. Time coverage, island pairs and observation effort differ. Plant data pooled. No claim of temporal increase or causal mechanism.')
(O/'island_host_early_new_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8')
fig,axs=plt.subplots(1,2,figsize=(10,4.8))
q=s[s.cohort.str.contains('common59')]
xx=np.arange(2)
axs[0].bar(xx-.17,q.observed_species_equal*100,width=.34,label='Observed',color='#177e89')
axs[0].bar(xx+.17,q.null_species_equal*100,width=.34,label='Null mean',color='#b8cece')
axs[0].set_xticks(xx,['Earlier records','New records']);axs[0].set_ylabel('Matching shared host (%)');axs[0].set_title('Same 59 species, equal species weights');axs[0].legend()
axs[1].hist(delta*100,bins=15,color='#b8cece',edgecolor='white');axs[1].axvline(observed*100,color='#b84b36',lw=2,label=f'Observed {observed*100:.1f} pp')
axs[1].set_xlabel('New minus earlier matching (percentage points)');axs[1].set_ylabel('Randomized networks');axs[1].legend();axs[1].set_title('Paired family-preserving null comparison')
fig.text(.02,.02,'Retrospective comparison; earlier 2001–2010, new first records 2011–2022; plants pooled 2001–2024.',fontsize=8)
fig.tight_layout(rect=[0,.06,1,1]);fig.savefig(O/'island_host_early_new.png',dpi=180);plt.close(fig)
print(s.to_string(index=False));print(json.dumps(v))

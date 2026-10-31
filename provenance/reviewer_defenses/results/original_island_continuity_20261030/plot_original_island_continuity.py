from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
source=O/'global_original_sequence_observation_opportunity.csv'
p=pd.read_csv(source,dtype={'landmass_id':str,'prior_landmass_id':str})
assert len(p)==12028 and not p.duplicated(['speciesKey','landmass_id']).any()
assert p.prior_record_year.lt(p.first_record_year).all()
bins=[0,10,100,500,np.inf];labels=['10 km未満','10–100 km','100–500 km','500 km以上']
p['distance_band']=pd.cut(p.nearest_prior_distance_km,bins,right=False,labels=labels)
assert p.distance_band.notna().all()
cohorts=[('all','全記録',p),('prior_opportunity','過去の観測被覆あり',p[p.prior_opportunity]),
         ('both_period_opportunity','前後の観測被覆あり',p[p.both_period_opportunity])]
rows=[];origins=[];species_rows=[]
for name,label,q in cohorts:
    for band in labels:
        flag=q.distance_band.eq(band)
        per_species=flag.groupby(q.speciesKey).mean()
        for key,part in q.groupby('speciesKey'):
            species_rows.append(dict(cohort=name,speciesKey=int(key),distance_band=band,
                events=int(part.distance_band.eq(band).sum()),total_events=len(part)))
        rows.append(dict(cohort=name,distance_band=band,events=int(flag.sum()),total_events=len(q),
            event_weighted_share=flag.mean(),species_equal_share=per_species.mean(),
            species=q.speciesKey.nunique(),islands=q.landmass_id.nunique()))
    for origin,r in q.groupby('prior_landmass_is_island'):
        origins.append(dict(cohort=name,prior_record_on_island=bool(origin),events=len(r),
            near100_events=int(r.nearest_prior_distance_km.lt(100).sum()),
            far500_events=int(r.nearest_prior_distance_km.ge(500).sum()),
            rerecorded_events=int(r.re_recorded_within2years.sum())))
s=pd.DataFrame(rows);s.to_csv(O/'global_original_continuity_distance_summary.csv',index=False)
pd.DataFrame(species_rows).to_csv(O/'global_original_continuity_species_counts.csv',index=False)
pd.DataFrame(origins).to_csv(O/'global_original_continuity_origin_summary.csv',index=False)
plt.rcParams['font.family']='Meiryo'
fig,axes=plt.subplots(1,2,figsize=(12,5.8),sharey=True)
colors=['#196b76','#71b7b2','#e5ba72','#b85b46']
for ax,weight,title in zip(axes,['event_weighted_share','species_equal_share'],['記録イベントを等重み','チョウ種を等重み']):
    bottom=np.zeros(3)
    for band,color in zip(labels,colors):
        values=np.array([s[(s.cohort==name)&(s.distance_band==band)][weight].iloc[0]*100 for name,_,_ in cohorts])
        ax.bar(range(3),values,bottom=bottom,color=color,label=band,width=.65)
        for i,value in enumerate(values):
            if value>=7:ax.text(i,bottom[i]+value/2,f'{value:.1f}%',ha='center',va='center',fontsize=9,color='white' if color in [colors[0],colors[3]] else '#222222')
        bottom+=values
    ax.set_xticks(range(3),[f'{label}\n{len(q):,}件 / {q.speciesKey.nunique()}種' for _,label,q in cohorts],fontsize=9)
    ax.set_title(title);ax.set_ylim(0,100);ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('以前の最寄り記録からの距離別割合（%）')
fig.suptitle('島の新規記録は、以前の分布記録からどれだけ離れていたか',fontsize=15)
handles,legendlabels=axes[0].get_legend_handles_labels();fig.legend(handles,legendlabels,loc='lower center',bbox_to_anchor=(.5,.15),ncol=4,frameon=False)
fig.text(.035,.035,'2011–2022年の島ごとの初記録。2001年以降・前年までの元座標記録を参照。\n観測被覆あり：過去5年に他の候補チョウ5種以上・2年以上。前後：さらに後2年とも記録あり。\n記録上の近接性であり、移動経路・自然分散・到来年の証明ではない。',fontsize=9)
fig.tight_layout(rect=[0,.23,1,.93]);fig.savefig(O/'global_original_island_continuity.png',dpi=170)
assert np.allclose(s.groupby('cohort').event_weighted_share.sum(),1)
assert np.allclose(s.groupby('cohort').species_equal_share.sum(),1)
receipt=dict(status='SUMMARY_RECONCILED',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
             events=len(p),species=int(p.speciesKey.nunique()),islands=int(p.landmass_id.nunique()),
             checks=['unique species-island events','strict earlier-year reference','all distances assigned',
                     'both weighting schemes sum to one'],
             limits=['no inference of observed movement','first record not true arrival','no absence from missing records'])
(O/'global_original_continuity_figure_audit.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(s.to_dict('records'),ensure_ascii=True))

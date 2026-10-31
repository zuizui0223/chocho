from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
p=pd.read_parquet(R/'work/temporal/forecast/species_cell_period_panel.parquet')
a=p[p.period.eq(1)].copy();b=p[p.period.eq(2)][['speciesKey','cell','recording_probability10','days']].rename(columns={'recording_probability10':'later_probability','days':'later_days'})
a=a.merge(b,on=['speciesKey','cell'],validate='one_to_one').dropna(subset=['recording_probability10','later_probability'])
a['change']=a.later_probability-a.recording_probability10
g=pd.read_parquet(W/'cell_geography.parquet')
a=a.merge(g[['cell','landmass_id','landmass_area_km2','distance_port_2019_km']],on='cell',validate='many_to_one')
q=pd.read_parquet(W/'recent_island_features.parquet');a=a.merge(q[['speciesKey','cell','period','log_recent_other_island_distance','no_recent_other_island']],on=['speciesKey','cell','period'],validate='one_to_one')
q=pd.read_parquet(W/'sea_route_features.parquet');a=a.merge(q[['speciesKey','cell','period','sampled_max_marine_gap_km','route_endpoint_unknown']],on=['speciesKey','cell','period'],validate='one_to_one')
a['recent_island_distance_km']=np.expm1(a.log_recent_other_island_distance)
a.loc[a.no_recent_other_island.eq(1),'recent_island_distance_km']=np.nan
a.loc[a.route_endpoint_unknown.eq(1),'sampled_max_marine_gap_km']=np.nan
a=a[a.landmass_area_km2.lt(7e6)].copy()
assert not a.duplicated(['speciesKey','cell']).any()
features={'distance_port_2019_km':'港までの距離（2019年版）','recent_island_distance_km':'直近に記録された別島までの距離','sampled_max_marine_gap_km':'最寄り既往記録からの最大連続海上区間'}
bins=[0,10,25,50,100,250,500,np.inf];labels=['0–10','10–25','25–50','50–100','100–250','250–500','500以上']
curves=[];corr=[]
for feature in features:
 v=a.dropna(subset=[feature]).copy();v['distance_bin']=pd.cut(v[feature],bins,labels=labels,right=False,include_lowest=True)
 for group,t in [('all_islands',v),('excluding_GB_Ireland',v[~v.landmass_id.isin(['15','26'])])]:
  for band,h in t.groupby('distance_bin',observed=True):
   bysp=h.groupby('speciesKey').change.mean()
   curves.append(dict(feature=feature,group=group,distance_bin=str(band),mean_change=bysp.mean(),species=len(bysp),cells=h.cell.nunique(),rows=len(h),landmasses=h.landmass_id.nunique()))
 for (sp,island),h in v.groupby(['speciesKey','landmass_id']):
  if len(h)<10 or h[feature].nunique()<3 or h.change.nunique()<3:continue
  rho=float(spearmanr(h[feature],h.change).statistic)
  corr.append(dict(feature=feature,speciesKey=sp,landmass_id=island,rows=len(h),rho=rho))
c=pd.DataFrame(curves);c.to_csv(O/'island_correlation_distance_bins.csv',index=False)
d=pd.DataFrame(corr);d.to_csv(O/'island_within_species_landmass_correlations.csv',index=False)
s=[]
for f,h in d.groupby('feature'):
 z=h.groupby('speciesKey').rho.mean()
 s.append(dict(feature=f,species=len(z),species_island_groups=len(h),landmasses=h.landmass_id.nunique(),mean_within_island_rho_species_equal=z.mean(),median_species_rho=z.median(),positive_species=int(z.gt(0).sum()),negative_species=int(z.lt(0).sum())))
pd.DataFrame(s).to_csv(O/'island_correlation_summary.csv',index=False)
plt.rcParams['font.family']='Meiryo';fig,axes=plt.subplots(1,3,figsize=(16,6),sharey=True)
for ax,(f,label) in zip(axes,features.items()):
 for group,color,text in [('all_islands','#347b9b','島側すべて'),('excluding_GB_Ireland','#bd783b','GB島・アイルランド島を除く')]:
  h=c[(c.feature==f)&(c.group==group)].set_index('distance_bin').reindex(labels)
  ax.plot(range(len(labels)),h.mean_change,marker='o',color=color,label=text)
 ax.axhline(0,color='gray',ls='--',lw=1);ax.set_title(label,fontsize=11);ax.set_xticks(range(len(labels)),labels,rotation=45);ax.set_xlabel('距離区分（km）');ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('記録確率の変化（後期 − 前期）');axes[0].legend(fontsize=9)
fig.suptitle('島での記録変化と距離の関連：まず未調整のパターンを確認',fontsize=15)
fig.text(.04,.02,'2015–2018年 → 2019–2022年。各区分は種ごとの平均を等重みで集計。区分ごとに種・島の構成が異なる。\n観測努力・気候・土地利用は未調整。港は2019年版の位置。変化は到来・定着を意味しない。',fontsize=10)
fig.tight_layout(rect=[0,.14,1,.94]);fig.savefig(O/'island_correlation_distance_patterns.png',dpi=160)
(O/'island_correlation_audit.json').write_text(json.dumps(dict(status='DESCRIPTIVE_ONLY',rows=len(a),species=int(a.speciesKey.nunique()),cells=int(a.cell.nunique()),landmasses=int(a.landmass_id.nunique()),missing_features={f:int(a[f].isna().sum()) for f in features},correlation_min_rows=10,correlation_min_unique_values=3,adjusted=False,limits=['distance bins differ in composition','within-species-landmass correlations not adjusted for environment or observation effort','regional host covariate is constant within species-region; cannot test local host aggregation','2019 port positions contemporaneous, not historical forecast inputs']),indent=2),encoding='utf-8')
print(pd.DataFrame(s).to_string(index=False))

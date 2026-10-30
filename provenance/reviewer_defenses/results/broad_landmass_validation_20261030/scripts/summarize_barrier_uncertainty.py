from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast'
geo=pd.read_parquet(R/'work/barriers/cell_geography.parquet');rows=[]
for mode in ['temporal','spatial']:
 receipt=O/f'forecast_barrier_{mode}_results.json'
 if not receipt.exists() or json.loads(receipt.read_text())['status']!='COMPLETE':continue
 p=pd.read_parquet(F/f'predictions_barrier_{mode}.parquet').merge(geo,on='cell',validate='many_to_one')
 groups={'all':p,'continental':p[p.landmass_area_km2.ge(7e6)],'noncontinental':p[p.landmass_area_km2.lt(7e6)],'islands_excluding_two':p[p.landmass_area_km2.lt(7e6)&~p.landmass_id.isin(['15','26'])]}
 for group,q in groups.items():
  species=sorted(q.species.unique());blocks=sorted(q.spatial_block.unique());rng=np.random.default_rng(20261030)
  count=pd.crosstab(q.species,q.spatial_block).reindex(index=species,columns=blocks,fill_value=0).to_numpy()
  sw=rng.multinomial(len(species),np.full(len(species),1/len(species)),size=2000)
  bw=rng.multinomial(len(blocks),np.full(len(blocks),1/len(blocks)),size=2000)
  denom=bw@count.T;weights=sw*(denom>0)
  for reference in ['M5_host_composition','B1_history']:
   for metric in ['mse','mae']:
    ea=q.B2_landmass-q.target;er=q[reference]-q.target
    delta=ea**2-er**2 if metric=='mse' else abs(ea)-abs(er)
    sums=q[['species','spatial_block']].assign(delta=delta).groupby(['species','spatial_block']).delta.sum().unstack(fill_value=0).reindex(index=species,columns=blocks,fill_value=0).to_numpy()
    means=np.divide(bw@sums.T,denom,out=np.zeros(denom.shape),where=denom>0)
    samples=(means*weights).sum(axis=1)/weights.sum(axis=1)
    low,high=np.quantile(samples,[.025,.975]);sp=delta.groupby(q.species).mean()
    rows.append(dict(mode=mode,group=group,reference=reference,added='B2_landmass',metric=metric,difference=sp.mean(),ci_low=low,ci_high=high,species=len(species),blocks=len(blocks),cells=q.cell.nunique(),rows=len(q),species_improved=int((sp<0).sum())))
result=pd.DataFrame(rows);result.to_csv(O/'barrier_paired_uncertainty.csv',index=False)
plt.rcParams['font.family']='Meiryo'
for mode in result['mode'].unique():
 fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
 labels=['全体 (3,495セル)','大陸側 (2,423セル)','島側 (993セル)','島側：GB島・アイルランド島除外* (134セル)'];groups=['all','continental','noncontinental','islands_excluding_two']
 for ax,metric in zip(axes,['mse','mae']):
  for ref,offset,color,marker,label in [('M5_host_composition',-.12,'#416a8c','o','B2 − 既存M5'),('B1_history',.12,'#b67632','s','B2 − 記録距離・密度B1')]:
   for i,g in enumerate(groups):
    a=result[(result['mode']==mode)&(result.group==g)&(result.reference==ref)&(result.metric==metric)].iloc[0]
    ax.plot([a.ci_low,a.ci_high],[i+offset]*2,color=color,lw=2)
    ax.scatter(a.difference,i+offset,color=color,marker=marker,s=28,label=label if i==0 else None)
  ax.axvline(0,color='#777',ls='--',lw=1);ax.set_yticks(range(4),labels);ax.set_xlabel(metric.upper()+' の差（負 = 誤差が小さい）');ax.spines[['top','right']].set_visible(False);ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2))
 axes[0].invert_yaxis();axes[1].legend(loc='best',fontsize=9)
 fig.suptitle(('時間外' if mode=='temporal' else '空間外')+'検証：陸塊情報を含むモデルの予測誤差',fontsize=14)
 fig.text(.04,.035,'95%区間：種×空間ブロックの2,000回再標本化。固定予測の不確実性であり、再学習・観測誤差は含まない。\n*事後的な部分集合。目的変数は標準化した記録確率で、到来・定着の確率ではない。',fontsize=9)
 fig.tight_layout(rect=[0,.13,1,.93]);fig.savefig(O/f'barrier_{mode}_uncertainty.png',dpi=160);plt.close(fig)
print(result[(result['mode']=='temporal')&(result.metric=='mse')].to_string(index=False))

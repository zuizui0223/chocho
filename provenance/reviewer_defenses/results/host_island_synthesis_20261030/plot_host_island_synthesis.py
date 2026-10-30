"""Synthesize verified positive results; preserve distinct estimands and cohorts."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
sources=['forecast_temporal_results.json','forecast_spatial_results.json','island_host_network_audit.json','island_host_family_null_audit.json','island_host_family_null_draws.csv']
v=[json.loads((O/f).read_text()) for f in sources[:4]]
improvements=[]
for a in v[:2]:
    assert a['status']=='COMPLETE'
    m=a['results'];improvements.append(100*(1-m['M5_host_composition']['species_equal_mse']/m['M3_both']['species_equal_mse']))
n=v[2];null=v[3];draws=pd.read_csv(O/sources[4])
np.testing.assert_allclose(draws.fraction_matched.mean(),null['null_mean_fraction'])
assert len(draws)==199 and null['candidate_species']==225
rows=[('A','temporal',improvements[0],'percent MSE reduction'),('A','spatial',improvements[1],'percent MSE reduction'),
      ('B','raw',n['spearman_host_new'],'rank correlation'),('B','distance_and_earlier_adjusted',n['partial_rank_host_new'],'partial rank correlation'),
      ('C','catalogue',null['observed_fraction']*100,'percent matched'),('C','family_null_mean',null['null_mean_fraction']*100,'percent matched')]
pd.DataFrame(rows,columns=['panel','comparison','value','unit']).to_csv(O/'host_island_synthesis_values.csv',index=False)
plt.rcParams.update({'font.family':'Yu Gothic','font.size':11,'axes.unicode_minus':False})
fig,axes=plt.subplots(1,3,figsize=(13,5.4));fig.subplots_adjust(left=.07,right=.97,bottom=.29,top=.73,wspace=.42)
labels=[['時間外','空間外'],['単純相関','距離・以前の\nチョウ相を調整'],['実際の\n食草対応','科構成を保つ\n置換の平均']]
values=[improvements,[n['spearman_host_new'],n['partial_rank_host_new']],[null['observed_fraction']*100,null['null_mean_fraction']*100]]
titles=['A　食草情報の予測価値\n155種','B　食草と新規チョウ相の関連\n51島・1,275島ペア','C　具体的な食草関係との対応\n82種・1,163種×島ペア']
units=['MSE改善率（%）','順位相関','共通する対応食草がある割合（%）']
for ax,lab,vals,title,unit,limit in zip(axes,labels,values,titles,units,[5,.4,22]):
    bars=ax.bar([0,1],vals,width=.55,color=['#285d78','#a3bbc8'],edgecolor='#285d78')
    for b,x in zip(bars,vals):ax.text(b.get_x()+b.get_width()/2,x+limit*.03,f'{x:.3f}' if unit=='順位相関' else f'{x:.2f}%',ha='center')
    ax.set(xticks=[0,1],xticklabels=lab,ylim=(0,limit),ylabel=unit,title=title)
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('食草情報の予測価値と、島々の食草・チョウ分布の対応',fontsize=16,y=.96)
fig.text(.04,.08,'候補225種を維持し、解析ごとに利用可能な対象を使用。3パネルは異なる対象・指標の結果。\nB・Cの植物記録は2001–2024年集約、新規チョウは2011–2022年の初記録。到来年ではない。\nCは199回の構造的置換の平均。輸送経路、先行する食草の効果、島間中継を実証した図ではない。',fontsize=10)
fig.savefig(O/'host_island_synthesis.png',dpi=170);plt.close(fig)
(O/'host_island_synthesis_audit.json').write_text(json.dumps(dict(status='SOURCE_VALUES_RECONCILED',
    sources={s:hashlib.sha256((O/s).read_bytes()).hexdigest() for s in sources},
    limits='Distinct cohorts and estimands; not independent replications of a causal pathway. Positive synthesis retains qualifications in accompanying report.'),indent=2))
print(rows)

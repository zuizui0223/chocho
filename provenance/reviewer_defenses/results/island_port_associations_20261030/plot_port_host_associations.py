from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
d=pd.read_csv(O/'port_host_adjusted_associations_interval_sensitivity.csv')
plt.rcParams['font.family']='Meiryo';fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
labels=['島内比較：93種・22島','GB島・アイルランド島除外：93種・20島','提供者優占度差≤0.1：82種・10島']
for ax,f,title in zip(axes,['log_port_distance','port_x_introduced'],['港までの距離との関連\n（導入食草割合50%での係数）','導入食草割合による関連の違い\n（港距離 × 導入食草割合）']):
 for i,sub in enumerate(['all_islands','excluding_GB_Ireland','stable_provider']):
  r=d[(d.subset==sub)&(d.feature==f)].iloc[0];ax.plot([r.cluster_t_low,r.cluster_t_high],[i,i],color='#427a93',lw=2);ax.scatter(r.coefficient,i,color='#427a93')
 ax.axvline(0,color='gray',ls='--');ax.set_yticks(range(3),labels);ax.set_title(title,fontsize=11);ax.set_xlabel('後期記録確率との調整後の関連係数');ax.spines[['top','right']].set_visible(False)
axes[0].invert_yaxis();fig.suptitle('港・食草との関連：同じ種・島の中で比較',fontsize=14)
fig.text(.03,.02,'距離は ln(1+km)。種×島固定効果、初期値・観測努力・気候・土地利用を調整。\n区間：二方向クラスタ分散＋t近似。因果効果ではない。優占度が近くても提供者自体は入れ替わり得る。',fontsize=9)
fig.tight_layout(rect=[0,.15,1,.91]);fig.savefig(O/'port_host_adjusted_associations.png',dpi=160)

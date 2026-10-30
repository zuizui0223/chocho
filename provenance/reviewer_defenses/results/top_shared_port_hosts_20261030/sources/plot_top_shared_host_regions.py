from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
c=pd.read_csv(O/'top_shared_host_region_correlations.csv');c=c[c.outcome.eq('change')]
order=['Medicago sativa','Robinia pseudoacacia','Plantago lanceolata','Capsella bursa-pastoris','Hesperis matronalis','Malus domestica']
plt.rcParams.update({'font.family':'Yu Gothic','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(10,6))
for i,host in enumerate(order):
    q=c[c.host_name.eq(host)].sort_values('species');offset=np.linspace(-.13,.13,len(q)) if len(q)>1 else np.zeros(len(q))
    ax.scatter(q.rho,i+offset,s=45,color='#397a85',alpha=.8,zorder=3)
    ax.scatter(q.rho.mean(),i,marker='D',s=60,color='#c16a2b',zorder=4)
    ax.text(1.05,i,f'{len(q)}種',va='center')
ax.axvline(0,color='#888888',lw=1,linestyle='--')
ax.set_yticks(range(len(order)),order);ax.invert_yaxis();ax.set_xlim(-1,1.23)
ax.set_xticks([-1,-.5,0,.5,1]);ax.set_xlabel('港での食草記録割合とチョウの記録確率変化：種別の順位相関 ρ')
ax.set_title('上位の共通食草：チョウ種によって関連の方向が異なる',loc='left',pad=18)
fig.text(.03,.075,'青丸：チョウ1種の地域間相関　橙の菱形：食草ごとの相関平均。各相関は5–12地域。',fontsize=10)
fig.text(.03,.035,'植物は2001–2024年の記録を集約。未補正の探索的照合で、先行する食草導入や拡大の因果効果ではない。',fontsize=10)
fig.subplots_adjust(left=.28,right=.95,bottom=.20,top=.87)
fig.savefig(O/'top_shared_host_region_correlations.png',dpi=180)
print('Saved24species-host correlations; means are descriptive, no independence-based intervals')

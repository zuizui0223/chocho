from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
plt.rcParams.update({'font.family':'Yu Gothic','font.size':11,'axes.unicode_minus':False})
fig,axes=plt.subplots(1,2,figsize=(12,5),sharex=True,sharey=True)
fig.subplots_adjust(left=.22,right=.96,bottom=.28,top=.75,wspace=.17)
refs=['M5_host_composition','B5_sea_terrain','T1_port_connectivity'];rows=[]
for ax,mode,title in zip(axes,['temporal','spatial'],['時間外','空間外＋時間外']):
    d=pd.read_csv(O/f'joint_{mode}_paired_uncertainty.csv')
    d=d[d.cohort.eq('strict_islands')].set_index('reference').loc[refs]
    assert (d.species==95).all() and (d.rows==41928).all()
    for i,(ref,r) in enumerate(d.iterrows()):
        ax.plot([r.interval_low,r.interval_high],[i,i],color='#285d78',lw=2)
        ax.plot(r.mse_difference,i,'o',color='#285d78',ms=7)
        ax.text(r.mse_difference,i-.18,f'{r.mse_difference:+.6f}',ha='center',fontsize=10)
        rows.append(dict(mode=mode,reference=ref,difference=r.mse_difference,low=r.interval_low,high=r.interval_high))
    ax.axvline(0,color='#444444',lw=1)
    ax.set(yticks=range(3),yticklabels=['既存の食草モデル M5','自然地理を含む B5','港指標を含む T1'],ylim=(2.5,-.6),xlim=(-.0075,.002),
           title=title,xlabel='統合J1 − 比較モデルのMSE\n負：J1の誤差が小さい')
    ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('島での統合モデル：改善方向は共通するが、差の区間はゼロをまたぐ',fontsize=15,y=.95)
fig.text(.04,.08,'原座標帰属が一致する25島・95種・41,928行。同じ評価行で比較。\n線：種×空間ブロック2,000回再抽出の中央95%区間。予測は固定し、モデル再学習や資料誤差は含まない。\n港指数の支持は5島に偏る。自然地理の代理指標と記録確率の比較であり、輸送・中継の因果効果ではない。',fontsize=10)
fig.savefig(O/'joint_island_comparison.png',dpi=170);plt.close(fig)
pd.DataFrame(rows).to_csv(O/'joint_island_figure_values.csv',index=False)

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs'
plt.rcParams.update({'font.family':'Yu Gothic','font.size':11,'axes.unicode_minus':False})
t=pd.read_csv(O/'transport_temporal_comparison.csv');s=pd.read_csv(O/'transport_spatial_comparison.csv')
rows=[]
for mode,data,cohort,label in [('temporal',t,'all155','時間外・155種'),('spatial',s,'all155','空間外・155種'),
    ('spatial',s,'at_least_100km_from_training','空間外・学習地点から100 km以上'),('spatial',s,'strict_islands','空間外・25島')]:
    q=data[data.cohort.eq(cohort)&data.model.eq('T1_port_connectivity')]
    assert len(q)==1
    rows.append(dict(mode=mode,cohort=cohort,label=label,improvement_pct=q.improvement_vs_M5_pct.iloc[0],rows=q.rows.iloc[0]))
data=pd.DataFrame(rows);data.to_csv(O/'transport_external_figure_values.csv',index=False)
fig,(ax,bx)=plt.subplots(1,2,figsize=(13,5),gridspec_kw={'width_ratios':[1.4,1]})
fig.subplots_adjust(left=.25,right=.96,bottom=.25,top=.79,wspace=.72)
for i,r in data.iterrows():
    ax.plot(r.improvement_pct,i,'o',color='#285d78',markersize=8)
    ax.text(r.improvement_pct+(.025 if r.improvement_pct>=0 else -.025),i,f'{r.improvement_pct:+.3f}%',
            va='center',ha='left' if r.improvement_pct>=0 else 'right')
ax.set(yticks=range(4),yticklabels=data.label,xlim=(-.85,.5),ylim=(3.5,-.5),xlabel='M5に対するMSE改善率（%）\n負：悪化　／　正：改善',title='港の接続性を加えたT1')
ax.axvline(0,color='#444444',lw=1)
sup=pd.read_csv(O/'transport_evaluation_support.csv').set_index('cohort')
for i,(key,label) in enumerate([('strict_islands','島全体'),('strict_islands_without_GB_Ireland','英国・アイルランド以外')]):
    r=sup.loc[key];pct=100*r.supported_landmasses/r.landmasses
    bx.barh(i,100,color='#e5e5e5',height=.5)
    bx.barh(i,pct,color='#285d78',height=.5)
    bx.text(0,i+.37,f'{int(r.supported_landmasses)}/{int(r.landmasses)}島・{int(r.supported_rows):,}/{int(r.rows):,}行',fontsize=10,va='top')
bx.set(yticks=[0,1],yticklabels=['島全体','英国・アイルランド\n以外'],ylim=(1.8,-.6),xlim=(0,100),xlabel='指数を接続できた島の割合（%）',title='港指数の被覆')
for a in [ax,bx]:
    a.spines[['top','right']].set_visible(False)
fig.suptitle('港の接続性の追加改善は、空間外検証で再現しなかった',fontsize=16,x=.5,y=.96)
fig.text(.03,.08,'候補225種、既存の評価対象155種を維持。M5＝気候・土地利用・食草等。点は種等重みMSEの比較。\n記録確率の遡及的予測であり、輸送の因果効果ではない。被覆の少ない島への一般化はできない。',fontsize=10)
fig.savefig(O/'transport_external_comparison.png',dpi=170)
plt.close(fig)

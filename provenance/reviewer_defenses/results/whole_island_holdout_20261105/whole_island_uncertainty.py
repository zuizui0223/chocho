from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast'
source=F/'predictions_joint_island.parquet'
receipt=json.loads((O/'whole_island_holdout_validation.json').read_text())
assert receipt['status']=='PASS_SAVED_PREDICTION_REPLAY'
assert hashlib.sha256(source.read_bytes()).hexdigest()==receipt['prediction_sha256']
design=dict(status='POST_POINT_ESTIMATE_UNCERTAINTY_DESIGN',seed=20261105,draws=2000,
    method='Independent multinomial resampling of species and island labels; fixed paired predictions; percentile intervals. Species-equal and island-equal estimands retained separately.',
    limitations='Not a model-refit bootstrap or proof of independent islands/species. Geographic and phylogenetic dependencies may remain. This specification follows inspection of point estimates, not pre-result preregistration.')
lock=O/'whole_island_uncertainty_design.json'
if lock.exists():assert json.loads(lock.read_text())==design
else:lock.write_text(json.dumps(design,indent=2))
p=pd.read_parquet(source)
comparisons=[('B5_sea_terrain','M5_host_composition','Natural geography'),('J1_natural_transport','M5_host_composition','Geography + ports'),('J1_natural_transport','B5_sea_terrain','Ports added to geography')]
rows=[];influence=[]
for cohort,q in [('all_islands',p),('without_GB_Ireland',p[~p.landmass_id.isin(['15','26'])])]:
    sp=sorted(q.species.unique());isl=sorted(q.landmass_id.unique());rng=np.random.default_rng(design['seed'])
    def matrix(frame,col):
        return frame.groupby(['species','landmass_id'])[col].sum().unstack(fill_value=0).reindex(index=sp,columns=isl,fill_value=0).to_numpy()
    counts=matrix(q.assign(n=1),'n')
    sw=rng.multinomial(len(sp),np.ones(len(sp))/len(sp),size=design['draws'])
    iw=rng.multinomial(len(isl),np.ones(len(isl))/len(isl),size=design['draws'])
    for model,ref,label in comparisons:
        d=(q[model]-q.target)**2-(q[ref]-q.target)**2
        sums=matrix(q.assign(delta=d),'delta')
        for weighting,group in [('species_equal','species'),('island_equal','landmass_id')]:
            point=float(d.groupby(q[group]).mean().mean())
            if weighting=='species_equal':
                denom=iw@counts.T;num=iw@sums.T;weights=sw*(denom>0)
                direct=(sums.sum(axis=1)/counts.sum(axis=1)).mean()
            else:
                denom=sw@counts;num=sw@sums;weights=iw*(denom>0)
                direct=(sums.sum(axis=0)/counts.sum(axis=0)).mean()
            np.testing.assert_allclose(point,direct,atol=1e-15)
            means=np.divide(num,denom,out=np.zeros_like(num),where=denom>0)
            assert (weights.sum(axis=1)>0).all()
            samples=(means*weights).sum(axis=1)/weights.sum(axis=1)
            lo,hi=np.quantile(samples,[.025,.975])
            rows.append(dict(cohort=cohort,weighting=weighting,model=model,reference=ref,label=label,mse_difference=point,interval_low=lo,interval_high=hi,rows=len(q),species=len(sp),islands=len(isl),draws=design['draws']))
            for omitted in isl:
                keep=q.landmass_id.ne(omitted)
                influence.append(dict(cohort=cohort,weighting=weighting,model=model,reference=ref,omitted_island=omitted,mse_difference=float(d[keep].groupby(q.loc[keep,group]).mean().mean())))
out=pd.DataFrame(rows);out.to_csv(O/'whole_island_paired_uncertainty.csv',index=False)
pd.DataFrame(influence).to_csv(O/'whole_island_leave_one_out.csv',index=False)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True)
for i,cohort in enumerate(['all_islands','without_GB_Ireland']):
    for j,weighting in enumerate(['species_equal','island_equal']):
        ax=axes[i,j];s=out[(out.cohort==cohort)&(out.weighting==weighting)]
        for y,(_,row) in enumerate(s.iterrows()):
            color='#245980' if y<2 else '#a06320'
            ax.plot([row.interval_low,row.interval_high],[y,y],color=color,lw=2)
            ax.plot(row.mse_difference,y,'o' if y<2 else 's',color=color,ms=7)
        ax.axvline(0,color='#555555',ls='--',lw=1)
        ax.set_yticks(range(3),s.label);ax.invert_yaxis();ax.set_ylim(2.6,-.6)
        ax.set_title(('25 islands' if i==0 else '23 islands (excluding GB / Ireland)')+'\n'+('Species-equal' if j==0 else 'Island-equal'),fontsize=11)
        ax.grid(axis='x',alpha=.2)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        if i==1:ax.set_xlabel('MSE difference (added model - reference)\nNegative = lower prediction error')
fig.suptitle('Whole-island holdout: geography and incremental port information',fontsize=15,y=.98)
fig.subplots_adjust(left=.23,right=.98,bottom=.19,top=.86,wspace=.72,hspace=.62)
fig.text(.03,.035,'95% percentile ranges: 2,000 crossed species x island resamples of fixed predictions.\n95 species; 41,812 rows overall / 3,248 excluding GB and Ireland. Historical predictors retained; not record-free island forecasts.',fontsize=9)
fig.savefig(O/'whole_island_uncertainty.png',dpi=180);plt.close(fig)
lines=['# 島外予測の不確実性と港の追加価値','', '既存225候補・155種の主解析を変更せず、固定済みの25島・95種の追加評価を再標本化した。点推定を確認した後の追加診断であり、事前登録解析とは区別する。種と島を独立に復元抽出する交差再標本化2,000回。予測値は固定し、モデルの再学習はしていない。','', '|対象|重み|比較|MSE差|95%再標本化区間|','|---|---|---|---:|---|']
for x in rows:lines.append(f"|{x['cohort']}|{x['weighting']}|{x['label']}|{x['mse_difference']:+.6f}|[{x['interval_low']:+.6f}, {x['interval_high']:+.6f}]|")
lines+=['','MSE差は追加モデル−比較モデル。負なら改善。種等重みは各種内の全評価行を平均、島等重みは各島内の全評価行を平均する。再標本化では反対側のクラスタ重みも適用し、記録が残らない種／島はその反復の平均から除く。','', '区間は固定された予測誤差の標本構成への感度を表す。島間の地理依存、種間の系統依存、学習データやモデル推定の不確実性を全て解消する区間ではない。過去の記録は説明変数に保持される。港の輸送指標の被覆不足や後年の港位置情報も残る。自然地理群の改善を特定の海峡・山地の因果効果と同一視しない。','', '図：whole_island_uncertainty.png。島を一つずつ除外した全結果：whole_island_leave_one_out.csv。']
if ((out.interval_low<0)&(out.interval_high>0)).all():
    lines+=['','今回の全12比較で95%再標本化区間はゼロを含んだ。自然地理の点推定での改善は保持するが、島や種を越えた確実な改善とまでは判断しない。港追加の点推定は種等重みと島等重みで符号が異なり、港の追加価値も未確定。既存の植物・港・チョウの組成相関を否定する検定ではない。']
(O/'whole_island_uncertainty_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(out[['cohort','weighting','label','mse_difference','interval_low','interval_high']].to_string(index=False))

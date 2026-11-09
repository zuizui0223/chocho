from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
O=Path(__file__).resolve().parents[2]/'outputs'
frames=[];states={}
for mode in ['temporal','spatial']:
    receipt=O/f'forecast_intermediate_{mode}_results.json'
    states[mode]=json.loads(receipt.read_text())['status'] if receipt.exists() else 'NOT_COMPLETE'
    if states[mode]=='COMPLETE':
        assert json.loads((O/f'intermediate_{mode}_validation.json').read_text())['status']=='PASS_PAIRED_PREDICTION_REPLAY'
        frames.append(pd.read_csv(O/f'intermediate_{mode}_comparison.csv'))
d=pd.concat(frames,ignore_index=True)
lines=['# 中間島情報の追加予測価値','', '原225候補・適格155種を保持し、B5（食草・気候・土地利用・過去の記録接続性・海・セル内地形）に、幾何的な中間島数と判定不明フラグだけを追加したB6を比較した。RF設定と評価行を固定し、比較B5/M5のキー・応答・予測値の一致、特徴量ハッシュ、MSE/MAEの再計算を検証。','', '|評価|対象|比較対象|MSE改善率|MSE差の95%再標本化区間|','|---|---|---|---:|---|']
for row in d.itertuples():lines.append(f'|{row.mode}|{row.cohort}|{row.reference}|{row.improvement_pct:+.3f}%|[{row.interval_low:+.6f}, {row.interval_high:+.6f}]|')
lines+=['','MSE差はB6−比較対象で、負が改善。改善率は正が改善。区間は固定予測誤差の種×空間ブロック再標本化2,000回で、モデル再学習の不確実性や独立した島の反復ではない。','', '時間外評価の島部分集計は、島全体を学習から除いた評価ではない。中間島は直線経路の幾何情報で、中間島の占有・食草利用・実際の中継移動を示さない。海上距離などを既に含むB5に対する追加価値と、食草モデルM5との差を分ける。','', '評価状況：'+json.dumps(states)+'. 未完了の評価は結果として扱わない。']
if ((d.interval_low<0)&(d.interval_high>0)).all():lines+=['','現時点の完成済み比較では全区間がゼロを含み、中間島情報の追加効果は未確定。小さな改善の点推定を保持するが、島による中継機構の確認には使わない。']
(O/'intermediate_island_forecast_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
s=d[d.reference.eq('B5_sea_terrain')]
fig,ax=plt.subplots(figsize=(8,4.8))
labels={'all155':'All 155 species','area_supported_islands':'Area-supported islands','islands_without_GB_Ireland':'Islands excluding GB / Ireland','at_least100km':'>=100km from training cells'}
for y,r in enumerate(s.itertuples()):
    color='#245980' if r.mode=='temporal' else '#a06320'
    ax.plot([r.interval_low,r.interval_high],[y,y],color=color,lw=2);ax.plot(r.mse_difference,y,'o',color=color)
ax.set_yticks(range(len(s)),[f'{r.mode}: {labels[r.cohort]}' for r in s.itertuples()]);ax.invert_yaxis();ax.axvline(0,color='#555555',ls='--')
ax.set_title('Intermediate-island information added to natural geography')
ax.set_xlabel('MSE difference: B6 - B5 (negative = lower error)')
ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.spines[['top','right']].set_visible(False)
fig.subplots_adjust(left=.40,right=.96,bottom=.27,top=.86)
fig.text(.03,.06,'95% percentile ranges: 2,000 species x spatial-block resamples; fixed predictions.\nGeographic potential only, not observed stepping-stone use.',fontsize=9)
fig.savefig(O/'intermediate_island_forecast.png',dpi=180)
print(json.dumps(states))

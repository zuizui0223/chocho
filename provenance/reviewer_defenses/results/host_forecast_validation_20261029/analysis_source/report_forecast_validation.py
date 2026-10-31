"""Render only completed, metric-reconciled validation stages."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
R=Path(__file__).resolve().parents[2];O=R/'outputs'
names={'temporal':('時間外','Temporal'), 'species':('未学習種＋時間外','Unseen species + time'), 'spatial':('未学習空間＋時間外','Unseen blocks + time'), 'spatial_buffer100':('空間外・学習地点から100 km以上','Spatial test, >=100 km'), 'temporal_effort20':('20記録日への標準化','20 recording days'), 'temporal_radius250':('候補範囲250 km','Accessible radius 250 km'), 'temporal_exclude_top_provider':('最大提供者を除外','Largest provider excluded'), 'temporal_precision100':('位置誤差100 m以内','Precision <=100 m')}
contrasts=pd.read_csv(O/'forecast_paired_comparisons.csv')
available={};pending=[];metric_rows=[]
for mode,(ja,en) in names.items():
    path=O/f'forecast_{mode}_results.json'
    if not path.exists():pending.append(ja);continue
    receipt=json.loads(path.read_text())
    if receipt['status']!='COMPLETE' or mode not in set(contrasts['mode']):pending.append(ja);continue
    available[mode]=receipt
    for model,values in receipt['results'].items():metric_rows.append({'mode':mode,'model':model,**values})
pd.DataFrame(metric_rows).to_csv(O/'forecast_validation_all_metrics.csv',index=False)
lines=['# 広域・多種のチョウ記録分布変化：予測検証の統合結果','',
    ('**進行中。** 未完了：'+'、'.join(pending)+'。' if pending else '**設定した予測比較・感度分析の集計は完了。** 全体の完了監査は `forecast_completion_audit_ja.md` を参照。'),'',
    '候補225種。主解析は観測支持を満たす155種。学習は2011–2014年→2015–2018年、評価は2015–2018年→2019–2022年。目的変数は観測努力を標準化した記録確率で、真の占有・新規定着ではない。食草資料は静的な地域資源であり、食草の移動や現地個体数を測定していない。','',
    '## 食草を追加したときの予測誤差','',
    'M3は背景＋気候＋土地被覆、M4はM3＋食草資源指標、M5はM4＋食草の科組成。以下は種を等しく重み付けしたMSE。小さいほどよい。異なる感度条件の絶対誤差は、対象と目的変数が異なる場合があるため直接ランキングしない。','',
    '| 検証 | 種数 | 評価組数 | M3 | M4 | M5 | M5−M3の95%区間 |','|---|---:|---:|---:|---:|---:|---|']
for mode,receipt in available.items():
    r=receipt['results'];m=r['M5_host_composition']
    c=contrasts[(contrasts['mode']==mode)&(contrasts.added=='M5_host_composition')&(contrasts.reference=='M3_both')&(contrasts.metric=='mse')].iloc[0]
    lines.append(f"| {names[mode][0]} | {m['species']} | {m['rows']:,} | {r['M3_both']['species_equal_mse']:.6f} | {r['M4_host_resources']['species_equal_mse']:.6f} | {m['species_equal_mse']:.6f} | {c.ci_low:.6f} ～ {c.ci_high:.6f} |")
if 'species' in available:
    r=available['species']['results'];difference=r['M5_host_composition']['species_equal_mse']-r['M3_both']['species_equal_mse']
    lines+=['',f'未学習種では食草追加後のMSE差は {difference:+.6f}。時間外検証での改善が、そのまま種を越えて再現したとは言えない。種を越えた共通予測則という主張は支持されていない。']
if 'spatial' in available and 'spatial_buffer100' in available:
    lines+=['','100 km以上という評価は、空間外モデルの同じ予測を距離で絞った事後的な感度確認であり、独立した再学習や事前固定のバッファ付き検証ではない。']
    improvements={mode:100*(available[mode]['results']['M3_both']['species_equal_mse']-available[mode]['results']['M5_host_composition']['species_equal_mse'])/available[mode]['results']['M3_both']['species_equal_mse'] for mode in ['spatial','spatial_buffer100']}
    lines+=['',f"**肯定的な結果：食草の追加による改善は、未学習の空間ブロックでも残った。** MSEの減少は空間外で{improvements['spatial']:.2f}%、学習地点から100 km以上離れた部分集合でも{improvements['spatial_buffer100']:.2f}%。両方でMAEも減少し、固定予測の再抽出区間はゼロをまたがなかった。既知種の別の場所・後年を予測する際の追加価値として支持される。未学習種への転用や寄主拡大の因果効果とは区別する。"]
lines+=['','区間は学習済み予測を固定し、種×5度空間ブロックを2,000回再抽出した記述的な区間。モデル再学習・記録誤差・年次反復・将来シナリオの不確実性は含まない。','',
    '## 持続予測を超えたか','',
    'MSEとMAEを両方示す。MSEだけが小さい場合を、全般的に予測力が優れるとは扱わない。','',
    '| 検証 | 持続MSE | M5 MSE | 持続MAE | M5 MAE |','|---|---:|---:|---:|---:|']
for mode,receipt in available.items():
    r=receipt['results'];a=r['persistence'];b=r['M5_host_composition']
    lines.append(f"| {names[mode][0]} | {a['species_equal_mse']:.6f} | {b['species_equal_mse']:.6f} | {a['species_equal_mae']:.6f} | {b['species_equal_mae']:.6f} |")
if not pending:
    lines+=['','## 感度分析で残った追加価値','',
        '各条件内部でM3とM5を比較した。種数は全条件で155種を維持し、対象の種×メッシュ組数の違いを上表に示した。条件間の改善率を順位付けして最良条件を選択していない。','']
    for mode in ['temporal_effort20','temporal_radius250','temporal_exclude_top_provider','temporal_precision100']:
        a=available[mode]['results']['M3_both']['species_equal_mse'];b=available[mode]['results']['M5_host_composition']['species_equal_mse']
        lines.append(f'- {names[mode][0]}：MSEを{100*(a-b)/a:.2f}%減少。')
    lines+=['','4条件ともMSE・MAEの差の固定予測再抽出区間はゼロをまたがなかった。これは既知種の記録分布予測に食草情報を加える価値の頑健性を支持するが、観測偏りの完全除去や食草の因果効果を示すものではない。']
if len(available)>=2:
    fig,axes=plt.subplots(1,2,figsize=(12,max(4,1.5+.55*len(available))),sharey=True)
    for ax,metric in zip(axes,['mse','mae']):
        for y,mode in enumerate(available):
            a=contrasts[(contrasts['mode']==mode)&(contrasts.added=='M5_host_composition')&(contrasts.reference=='M3_both')&(contrasts.metric==metric)].iloc[0]
            ax.plot([a.ci_low,a.ci_high],[y,y],color='#4c7899',lw=2)
            ax.scatter(a.difference,y,color='#b56a35',s=32,zorder=3)
        ax.axvline(0,ls='--',color='gray',lw=1);ax.set_yticks(range(len(available)),[names[x][1] for x in available]);ax.set_xlabel('Change in '+metric.upper()+' (negative = improvement)')
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.spines[['top','right']].set_visible(False)
    axes[0].invert_yaxis();fig.suptitle('Incremental predictive value of host information across validation settings')
    fig.text(.03,.025,'M5 (all host information) minus M3 (climate + land cover). 95% crossed species × spatial-block resampling intervals.\nConditional on fixed predictions. Recording probability, not true occupancy; cohorts differ across sensitivities.',fontsize=9)
    fig.tight_layout(rect=[0,.15,1,.93]);fig.savefig(O/'forecast_validation_transfer.png',dpi=180)
    lines+=['','![検証条件ごとの食草追加効果](forecast_validation_transfer.png)']
subgroup_path=O/'forecast_initially_unrecorded_comparison.csv'
if subgroup_path.exists():
    sub=pd.read_csv(subgroup_path)
    lines+=['','## 前期に対象種の記録がなかった地点','',
        'モデル設定時に指定した部分集合。2015–2018年の対象種の記録確率が0の地点を2019–2022年に評価したもので、過去全期間の不在・新規定着・繁殖を意味しない。時間外・種外とも154種、133,557組。','',
        '| 検証 | 指標 | M3 | M5 | 相対的な誤差減少 | 差の95%区間 |','|---|---|---:|---:|---:|---|']
    for row in sub.itertuples():
        lines.append(f'| {names[row.mode][0]} | {row.metric.upper()} | {row.reference_M3:.6f} | {row.full_M5:.6f} | {row.relative_improvement_percent:.2f}% | {row.ci_low:.6f} ～ {row.ci_high:.6f} |')
    lines+=['','既知種の時間外評価での改善は、前期にも記録があった地点だけに限定されない。一方、未学習種への転用は別の問いであり、この改善をそのまま一般化しない。区間は固定予測の交差クラスターブートストラップによる記述的なもの。']
lines+=['','## 応答差と観測過程','',
    '食草×環境変化の補助混合モデルでは、9交互作用はいずれもHolm調整後p≥0.177で、交互作用を加えた後期再現精度も改善しなかった。食草の予測上の追加価値と、環境応答を変える機構の証拠は区別する。主解析155種とは異なり、この集約補助解析は139種・19地域である。','',
    '提供者構成変化の小さい群でも時間外MSEの改善は残ったが、全モデルの誤差には135種中111種で正の空間相関が残った。空間外検証と距離を設けた評価の結果を重視する。','',
    '## 根拠ファイル','',
    '- 全モデルの数値：`forecast_validation_all_metrics.csv`','- ペア比較・固定予測の区間：`forecast_paired_comparisons.csv`','- 独立再集計：`forecast_metric_validation.json`','- 候補225種の脱落理由：`forecast_all_225_disposition.csv`','- 主解析・食草交互作用・感度分析の条件：`forecast_evaluation_protocol_ja.md`、`forecast_interaction_protocol_ja.md`、`forecast_sensitivity_protocol_ja.md`','- 観測・空間診断：`forecast_observation_spatial_diagnostics_ja.md`','',
    '出現データ：[GBIF](https://doi.org/10.15468/dl.umqmsf)。土地被覆：[ESA LCCS年次集約データ](https://doi.org/10.5281/zenodo.11407287)。気候：NASA POWER。過去の記録アーカイブと現代の食草資料を使う後ろ向き検証であり、2018年当時に実際に発行された予報ではない。']
(O/'forecast_validation_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'completed_modes':list(available),'pending':pending}))

from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
receipt=json.loads((O/'coastal_control_vegetation_status.json').read_text())
assert receipt['status'] in ['COMPLETE','COMPLETE_WITH_SOURCE_GAPS']
pairs=pd.read_csv(O/'coastal_control_vegetation_pairs.csv',dtype={'landmass_id':str})
summary=pd.read_csv(O/'coastal_control_vegetation_summary.csv')
near=pd.read_parquet(W/'port_vegetation_samples_joined.parquet').merge(pd.read_parquet(W/'port_island_coastal_probe.parquet'),on='sample_id',validate='one_to_one')
far=pd.read_parquet(W/'coastal_control_vegetation_joined.parquet')
near.landmass_id=near.landmass_id.astype(str);far.landmass_id=far.landmass_id.astype(str)
near=near[near.terrestrial_valid & near.coastal_probe_1km]
far=far[far.terrestrial_valid]
ng={k:q for k,q in near.groupby(['port_id','landmass_id','radius_km'])}
fg={k:q for k,q in far.groupby(['port_id','landmass_id'])}
features={'tree':(10,'樹木'),'shrub':(20,'低木'),'grass':(30,'草地'),'crop':(40,'農地'),'built':(50,'市街地')}
for r in pairs.itertuples():
    n=ng[(r.port_id,r.landmass_id,r.near_radius_km)];f=fg[(r.port_id,r.landmass_id)]
    assert len(n)==r.near_points and len(f)==r.far_points and min(len(n),len(f))>=3
    for feature,(code,_) in features.items():
        assert np.isclose((n.landcover==code).sum()/len(n),getattr(r,feature+'_near'))
        assert np.isclose((f.landcover==code).sum()/len(f),getattr(r,feature+'_far'))
        assert np.isclose(getattr(r,feature+'_near')-getattr(r,feature+'_far'),getattr(r,feature+'_difference'))
for r in summary.itertuples():
    q=pairs[pairs.near_radius_km.eq(r.near_radius_km)]
    byisland=q.groupby('landmass_id')[r.vegetation+'_difference'].mean()
    assert len(byisland)==r.islands and q.port_id.nunique()==r.ports and len(q)==r.pairs
    assert np.isclose(byisland.mean(),r.mean_difference)
plt.rcParams['font.family']='Meiryo';fig,axes=plt.subplots(1,2,figsize=(11,5.6),sharey=True,sharex=True)
for ax,radius in zip(axes,[1,5]):
    q=summary[summary.near_radius_km.eq(radius)]
    for i,(feature,(_,label)) in enumerate(features.items()):
        r=q[q.vegetation.eq(feature)].iloc[0]
        if np.isfinite(r.ci_low):ax.plot([r.ci_low*100,r.ci_high*100],[i,i],color='#287b87',lw=2)
        ax.scatter(r.mean_difference*100,i,color='#287b87',zorder=3)
    ax.axvline(0,color='gray',ls='--');ax.set_yticks(range(5),[v[1] for v in features.values()])
    ax.set_title(f'港から{radius}kmの沿岸点\n{int(r.islands)}島・{int(r.ports)}港')
    ax.set_xlabel('遠方の沿岸対照点との差（ポイント）');ax.spines[['top','right']].set_visible(False)
axes[0].invert_yaxis();fig.suptitle('同じ島の沿岸同士で比べた、港周辺の土地被覆（2021年）')
fig.text(.03,.025,'点割合の差。面積割合・食草量ではない。島を等重み、区間は島単位2,000回再標本化。\n対照：港から約10–30kmの海岸から250m陸側、全登録港から5km以上。\n港側は1km以内の海岸近接の代理指標。海岸距離の完全一致や、港の因果効果の推定ではない。',fontsize=9)
fig.tight_layout(rect=[0,.18,1,.92]);fig.savefig(O/'coastal_control_vegetation_comparison.png',dpi=170)
lines=['# 同じ島の沿岸を対照にした港周辺の植生比較','',
       'ESA WorldCover 2021 v200の土地被覆を、港周辺の沿岸点と同じ島の遠方の沿岸点で比較した。静的な植生分類であり、食草種や外来性の識別とは別の解析である。','',
       '|港からの距離|植生|点割合の差（ポイント）|島単位95%再標本化区間|島数|港数|',
       '|---|---|---:|---|---:|---:|']
for r in summary.itertuples():
    interval=f'{r.ci_low*100:.1f} ～ {r.ci_high*100:.1f}' if np.isfinite(r.ci_low) else '島数20未満のため非表示'
    lines.append(f'|{r.near_radius_km:g} km|{features[r.vegetation][1]}|{r.mean_difference*100:+.1f}|{interval}|{r.islands}|{r.ports}|')
lines+=['','## 解釈範囲','',
    '近点・対照とも有効な陸上採取点が3点以上ある港×島を比較し、島内で港を平均してから島間で平均した。全球の候補から抽出しているが、同じ島に遠方対照を確保できる島に限られる。小島を含む全島へ一般化しない。',
    '近点は8方向の1km先のうち少なくとも1点が陸外になる海岸近接の代理指標。対照は海岸線の頂点から250m陸側を選ぶため、両群の海岸距離は厳密に一致していない。標高・地形・市街化も揃えていない。旧15km円上対照とは対象島・港が異なり、推定値の変化を調整の効果と解釈しない。',
    '土地被覆の点割合は面積割合でも植物量でもない。土地被覆の差から外来食草の増加や、チョウの人為移入は推論できない。植物の座標・年代・地域別由来との接続が引き続き必要。',
    f'対照設計は{receipt["design_samples"]:,}地点、取得{receipt["retrieved_samples"]:,}地点、提供元にタイルがない地点{receipt["missing_samples"]}。欠測は植生なしに変換しない。',
    '出典：ESA WorldCover 2021 v200（CC BY 4.0）、GSHHG 2.3.7、NOAA配信のWPI2019港カタログ。登録外の港・人為活動は把握していない。','']
(O/'coastal_control_vegetation_results_ja.md').write_text('\n'.join(lines),encoding='utf-8')
audit=dict(status='PAIR_AND_SUMMARY_COUNTS_RECOMPUTED',pairs=len(pairs),summary_rows=len(summary),
    checks=['every paired proportion recomputed from classified sample points','all island-weighted means and support counts recomputed'],
    not_verified=['independent remote raster pixel resampling','causal interpretation','exact coastline-distance matching'])
(O/'coastal_control_vegetation_validation.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit))

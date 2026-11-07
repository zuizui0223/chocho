from pathlib import Path
import hashlib,json,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast'
models=['M5_host_composition','B5_sea_terrain','T1_port_connectivity','J1_natural_transport']
comparisons=[(models[1],models[0]),(models[3],models[0]),(models[3],models[1]),(models[3],models[2])]
rows=[];sources=[];uncertainty=[]
for variant in ['temporal','temporal_effort20','temporal_exclude_top_provider']:
    receipt_path=O/f'forecast_joint_{variant}_results.json'
    receipt=json.loads(receipt_path.read_text());assert receipt['status']=='COMPLETE'
    path=F/f'predictions_joint_{variant}.parquet';p=pd.read_parquet(path)
    assert p.species.nunique()==155 and not p.duplicated(['speciesKey','cell']).any()
    assert np.isfinite(p[models+['target']].to_numpy()).all()
    for prefix,col in [('',models[0]),('transport_',models[2])]:
        ref=pd.read_parquet(F/f'predictions_{prefix}{variant}.parquet')
        np.testing.assert_array_equal(p[['speciesKey','cell','target',col]],ref[['speciesKey','cell','target',col]])
    if receipt['sensitivity']:
        lock=O/f"joint_{receipt['sensitivity']}_lock.json"
        assert hashlib.sha256(lock.read_bytes()).hexdigest()==receipt['joint_sensitivity_lock_sha256']
        panel=F/f"panel_{receipt['sensitivity']}.parquet"
        assert hashlib.sha256(panel.read_bytes()).hexdigest()==receipt['input_panel_sha256']
    for m in models:
        p[m+'_mse']=(p[m]-p.target)**2
        p[m+'_mae']=abs(p[m]-p.target)
        for metric in ['mse','mae']:
            np.testing.assert_allclose(p.groupby('species')[m+'_'+metric].mean().mean(),receipt['results'][m]['species_equal_'+metric],rtol=1e-12)
    scores=p.groupby('species')[[m+'_'+metric for m in models for metric in ['mse','mae']]].mean().mean()
    species=sorted(p.species.unique());blocks=sorted(p.spatial_block.unique())
    counts=pd.crosstab(p.species,p.spatial_block).reindex(index=species,columns=blocks,fill_value=0).to_numpy()
    rng=np.random.default_rng(20261105)
    sw=rng.multinomial(len(species),np.full(len(species),1/len(species)),size=2000)
    bw=rng.multinomial(len(blocks),np.full(len(blocks),1/len(blocks)),size=2000)
    denom=bw@counts.T;weights=sw*(denom>0)
    assert (weights.sum(axis=1)>0).all()
    for model,ref in comparisons:
        rows.append(dict(variant=variant,model=model,reference=ref,rows=len(p),species=p.species.nunique(),model_mse=scores[model+'_mse'],reference_mse=scores[ref+'_mse'],improvement_pct=100*(1-scores[model+'_mse']/scores[ref+'_mse']),model_mae=scores[model+'_mae'],reference_mae=scores[ref+'_mae']))
        delta=p[model+'_mse']-p[ref+'_mse']
        sums=p[['species','spatial_block']].assign(delta=delta).groupby(['species','spatial_block']).delta.sum().unstack(fill_value=0).reindex(index=species,columns=blocks,fill_value=0).to_numpy()
        point=delta.groupby(p.species).mean().mean()
        np.testing.assert_allclose(point,(sums.sum(axis=1)/counts.sum(axis=1)).mean(),atol=1e-15)
        means=np.divide(bw@sums.T,denom,out=np.zeros(denom.shape),where=denom>0)
        draws=(means*weights).sum(axis=1)/weights.sum(axis=1)
        low,high=np.quantile(draws,[.025,.975])
        uncertainty.append(dict(variant=variant,model=model,reference=ref,mse_difference=point,interval_low=low,interval_high=high,draws=2000,seed=20261105,spatial_blocks=len(blocks)))
    sources.append(dict(variant=variant,prediction_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),receipt_sha256=hashlib.sha256(receipt_path.read_bytes()).hexdigest()))
out=pd.DataFrame(rows);out.to_csv(O/'joint_observation_sensitivity_comparison.csv',index=False)
pd.DataFrame(uncertainty).to_csv(O/'joint_observation_sensitivity_uncertainty.csv',index=False)
(O/'joint_observation_sensitivity_validation.json').write_text(json.dumps(dict(status='PASS_THREE_VARIANT_PAIRED_REPLAY',sources=sources,checks=['155 species each variant','M5/T1 exact matching keys targets and predictions','sensitivity lock and source panel hashes','independent MSE/MAE replay'],limits='Comparison within each variant only; recording probability targets differ across variants. No model-refit uncertainty or complete provider-removal geography.'),indent=2))
lines=['# 自然地理・輸送モデルの観測努力と提供者感度','', '225候補を維持し、各条件で適格な155種を全て評価した。同じ条件のM5とT1の予測値・評価キー・応答の一致を確認し、B5/J1を固定RF設定で再学習した。値は種等重みMSEの改善率。正なら比較モデルより誤差が小さい。','', '|条件|追加モデル|比較対象|評価行|改善率|','|---|---|---|---:|---:|']
for x in rows:lines.append(f"|{x['variant']}|{x['model']}|{x['reference']}|{x['rows']}|{x['improvement_pct']:+.3f}%|")
lines+=['','条件ごとに応答となる標準化記録確率が異なるので、条件間の絶対MSE差を生物学的効果と解釈しない。最大提供者除外は観測努力・応答の変更であり、過去距離・候補領域・食草カタログ・地理の全入力を提供者なしで再構築した検証ではない。','', '本表は時間外の全155種の比較。島保持・空間外・種外の感度を追加実施したものではない。港距離と指数は実際の船舶運搬や食草利用を識別しない。','', '## 固定予測誤差の不確実性','', '種×空間ブロックを独立に再標本化した2,000回の95%百分位区間。通常解析の点推定を見た後の診断で、モデル再学習・系統や地理依存を完全に扱う区間ではない。負のMSE差が改善。','', '|条件|追加モデル|比較対象|MSE差|95%区間|','|---|---|---|---:|---|']
for x in uncertainty:lines.append(f"|{x['variant']}|{x['model']}|{x['reference']}|{x['mse_difference']:+.6f}|[{x['interval_low']:+.6f}, {x['interval_high']:+.6f}]|")
if all(x['interval_low']<0<x['interval_high'] for x in uncertainty):
    lines+=['','全12比較で区間はゼロを含んだ。全155種で一貫した追加改善は確認できず、点推定の符号も観測条件によって変わる。島部分集団の改善を広域全体に一般化しない。この結果は、独立した食草情報の予測改善や島間の組成相関を否定する検定ではない。']
(O/'joint_observation_sensitivity_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(out[['variant','model','reference','rows','improvement_pct']].to_string(index=False))

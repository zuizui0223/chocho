from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[2]; O=R/'outputs'; W=R/'work/barriers'; F=R/'work/temporal/forecast'
r=json.loads((O/'forecast_joint_island_results.json').read_text())
assert r['status']=='COMPLETE'
assert hashlib.sha256((O/'whole_island_holdout_lock.json').read_bytes()).hexdigest()==r['island_lock_sha256']
p=pd.read_parquet(F/'predictions_joint_island.parquet')
k=pd.read_parquet(W/'whole_island_holdout_test_keys.parquet')
keys=['speciesKey','cell']
assert not p.duplicated(keys).any()
pd.testing.assert_frame_equal(p[keys+['landmass_id','island_fold']].sort_values(keys).reset_index(drop=True),k.sort_values(keys).reset_index(drop=True))
assert (p.island_fold==p.landmass_id.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)).all()
panel=pd.read_parquet(F/'species_cell_period_panel.parquet',columns=keys+['period','recording_probability10'])
target=panel[panel.period.eq(2)][keys+['recording_probability10']]
q=p.merge(target,on=keys,validate='one_to_one')
np.testing.assert_array_equal(q.target,q.recording_probability10)
train=panel[panel.period.eq(0)&panel.recording_probability10.notna()].merge(panel[panel.period.eq(1)&panel.recording_probability10.notna()][keys],on=keys,validate='one_to_one')
parts=pd.read_parquet(W/'cell_land_overlap_parts.parquet')
for f in r['fold_support']:
    held=set(p.loc[p.island_fold.eq(f['fold']),'landmass_id'])
    touched=set(parts.loc[parts.landmass_id.isin(held),'cell'])
    assert int(train.cell.isin(touched).sum())==f['excluded_training_rows']
    assert int(p.island_fold.eq(f['fold']).sum())==f['test_rows']
models=['M5_host_composition','B5_sea_terrain','J1_natural_transport']
for m in models:
    assert np.isfinite(p[m]).all() and p[m].between(0,1).all()
    p[m+'_mse']=(p[m]-p.target)**2
    p[m+'_mae']=abs(p[m]-p.target)
    for metric in ['mse','mae']:
        np.testing.assert_allclose(p.groupby('species')[m+'_'+metric].mean().mean(),r['results'][m]['species_equal_'+metric],rtol=1e-12)
rows=[]
for name,sub in [('all_islands',p),('without_GB_Ireland',p[~p.landmass_id.isin(['15','26'])])]:
    for weighting,group in [('species_equal','species'),('island_equal','landmass_id')]:
        scores=sub.groupby(group)[[m+'_'+v for m in models for v in ['mse','mae']]].mean().mean()
        for m in models:
            rows.append(dict(cohort=name,weighting=weighting,model=m,rows=len(sub),species=sub.species.nunique(),islands=sub.landmass_id.nunique(),mse=scores[m+'_mse'],mae=scores[m+'_mae'],improvement_vs_M5_pct=100*(1-scores[m+'_mse']/scores[models[0]+'_mse'])))
out=pd.DataFrame(rows);out.to_csv(O/'whole_island_holdout_comparison.csv',index=False)
audit=dict(status='PASS_SAVED_PREDICTION_REPLAY',rows=len(p),species=p.species.nunique(),islands=p.landmass_id.nunique(),checks=['lock hash','unique test keys and locked folds','targets equal original period2','polygon intersection exclusion counts reconstructed','MSE and MAE agree with receipt'],prediction_sha256=hashlib.sha256((F/'predictions_joint_island.parquet').read_bytes()).hexdigest(),limits='Checks saved predictions and fit-code exclusion construction, not a second model refit. Point estimates only; no uncertainty or independent island sampling claim.')
(O/'whole_island_holdout_validation.json').write_text(json.dumps(audit,indent=2))
lines=['# 島全体を学習応答から除外した予測検証','', '225候補・155種の主解析を保持した追加評価。対象は25島・95種・41,812行。2019–2022年の標準化された記録確率を予測する。各島の面積に交差する全セルを学習から除き、M5（食草・気候・土地利用など）、B5（M5＋過去の記録接続性・海・地形）、J1（B5＋港・輸送指標）を各foldで再学習した。','', '|対象|重み|モデル|MSE|M5比改善率|','|---|---|---|---:|---:|']
for x in rows:lines.append(f"|{x['cohort']}|{x['weighting']}|{x['model']}|{x['mse']:.6f}|{x['improvement_vs_M5_pct']:+.2f}%|")
lines+=['','島等重みは各島内の行平均誤差を平均する。種等重みとは異なる評価量であり、相互に代用しない。英国・アイルランドを除く結果も全対象を保持して集計した。','', 'これは過去に未観測の島への予測ではない。予測時点までの島内記録は説明変数に利用する。港指標の被覆は英国・アイルランドに偏り、港位置・DEMの年代整合にも制約がある。自然地理群は記録接続性・地形・海を合わせたもので、個別の障壁効果や船による持ち込みを同定しない。改善率は点推定で、不確実性の検証は別途必要。','', '検証：固定テストキー、元パネルの応答、fold別の面積交差による学習除外行数、予測値からの種等重みMSE/MAEを独立に再計算した。']
(O/'whole_island_holdout_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(out.to_string(index=False))

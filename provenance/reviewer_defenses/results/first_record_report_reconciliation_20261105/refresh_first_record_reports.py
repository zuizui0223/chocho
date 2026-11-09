from pathlib import Path
import json, hashlib
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2]; O=R/'outputs'; W=R/'work/barriers'
p=pd.read_parquet(W/'host_butterfly_first_record_pairs.parquet')
e=pd.read_csv(O/'first_record_lag_event_summary.csv')
assert not e.duplicated(['species','region']).any()
assert len(e)==len(p.groupby(['species','Region']))==17
for row in e.itertuples():
    q=p[p.species.eq(row.species)&p.Region.eq(row.region)]
    assert len(q)==row.host_pairs
    dated=q[q.lower_plant.notna()&q.upper_plant.notna()&q.lower_butterfly.notna()&q.upper_butterfly.notna()]
    assert (q.lower_plant.notna()&q.upper_plant.notna()).sum()==row.dated_host_pairs
    assert row.any_host_recorded_before==q.ordering.eq('plant_recorded_first').any()
    assert row.any_host_recorded_after==q.ordering.eq('butterfly_recorded_first').any()
    if len(dated):
        np.testing.assert_allclose(row.earliest_host_lag_lower,dated.lower_butterfly.iloc[0]-dated.upper_plant.min())
        np.testing.assert_allclose(row.earliest_host_lag_upper,dated.upper_butterfly.iloc[0]-dated.lower_plant.min())
    else: assert pd.isna(row.earliest_host_lag_lower) and pd.isna(row.earliest_host_lag_upper)
t=pd.read_csv(O/'top_shared_host_first_record_geography_screened.csv')
for row in t.itertuples():
    q=p[p.host_name.eq(row.host_name)]
    assert len(q)==row.pairs
    for col,label in [('plant_recorded_first','plant_recorded_first'),('butterfly_recorded_first','butterfly_recorded_first'),('overlapping_intervals','overlapping_intervals'),('unresolved','unresolved_date')]:
        assert getattr(row,col)==q.ordering.eq(label).sum()
lines=['# 共通上位食草の初記録順序（地域不整合除外後）','','上位6食草に対応する現行19組のうち16組で植物が先に記録された。チョウ先行1組、区間重複1組、年代未確定1組。Canada/Greveの2組を隔離済み。残存する原典の全件確認や導入年代の確定ではない。','','|食草|組数|植物先行|チョウ先行|区間重複|年代未確定|','|---|---:|---:|---:|---:|---:|']
for x in t.itertuples():lines.append(f'|{x.host_name}|{x.pairs}|{x.plant_recorded_first}|{x.butterfly_recorded_first}|{x.overlapping_intervals}|{x.unresolved}|')
lines+=['','初記録の順序であり、船舶による持ち込み・同時輸送・食草への追随を識別しない。対応なしは年代資料の被覆不足で、食草利用の否定ではない。現在の港周辺記録とは時空間単位が異なる。','出典：Alien Species First Records v3.1、DOI 10.5281/zenodo.10039630。乱数を含むFirstRecordではなくFirstRecord_origの年代区間を使用。']
(O/'top_shared_host_first_record_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
f=O/'terrain_and_record_lag_results_ja.md';s=f.read_text(encoding='utf-8')
start=s.index('# 地形と初記録時間差の現状'); s=s[start:]
a=s.index('## 時間差：'); b=s.index('DEM出典：',a)
new='## 時間差：地域不整合を除外した現行集計\n\n2026-11-05再照合。独立初記録資料の265組、チョウ9種・9地域・17種×地域が対応した。植物先行178組、チョウ先行68組、年代未確定16組、区間重複3組。Canada/Greveの地域不整合35組を隔離しており、旧300組を現行値として使わない。\n\n種×地域単位では13/17群に植物先行記録、8/17群にチョウ先行記録があり、7群は両方を含む。食草が複数あるため、265組を独立した到来イベントとみなさない。最も早い対応食草との時間差図は14群が年代比較可能、3群が年代不足。線は原資料の年代幅であり、統計的信頼区間ではない。\n\n時間差は初記録の差で、実際の導入・到来・定着の待ち時間ではない。未収録の在来食草や栽培植物、記録開始の差が残る。原典全行の確認は未完了。この9種を広域225候補／155解析種の代表値として扱わず、予測セルへの時間差説明変数としても使用していない。\n\n出典：Alien Species First Records v3.1、DOI 10.5281/zenodo.10039630。FirstRecord_origの年代区間を使用。元資料の地域名は予測セル・評価島と同一視しない。メタデータとワークブックの利用条件の差は出典監査に保持する。\n\n'
f.write_text(s[:a]+new+s[b:],encoding='utf-8')
audit=dict(status='PASS_CURRENT_DATE_REPORT_REPLAY',pairs=len(p),groups=len(e),dated_groups=int(e.earliest_host_lag_lower.notna().sum()),top_host_pairs=int(t.pairs.sum()),pair_sha256=hashlib.sha256((W/'host_butterfly_first_record_pairs.parquet').read_bytes()).hexdigest(),limits='Report/interval reconciliation only; no new independent introduction-date verification.')
(O/'first_record_report_replay.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit))

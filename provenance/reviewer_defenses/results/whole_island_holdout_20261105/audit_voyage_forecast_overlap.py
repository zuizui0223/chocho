from pathlib import Path
import subprocess,hashlib,json,io
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';repo=R/'work/remote-sensing-runner'
sha='bb9f7d08e09416635725359eb4c7c3c48b866a33'
prefix='provenance/reviewer_defenses/results/island_observational_robustness_20261031/'
sources={}
def read(name):
    raw=subprocess.check_output(['git','show',f'{sha}:{prefix}{name}'],cwd=repo)
    sources[name]=dict(commit=sha,path=prefix+name,sha256=hashlib.sha256(raw).hexdigest())
    return json.loads(raw) if name.endswith('.json') else pd.read_csv(io.BytesIO(raw),dtype={'island_a':str,'island_b':str})
ais=read('ais_2020_wpi51_island_5km.audit.json')
voc=read('voc_57669_17c_interisland_routes_v1.audit.json')
a=read('ais_2020_wpi51_island_5km.island_pairs.csv')
v=read('voc_57669_17c_interisland_routes_v1.positive_historic_links.csv')
local_hash=hashlib.sha256((O/'island_host_network_pairs.csv').read_bytes()).hexdigest()
assert ais['island_pair_source_sha256']==local_hash==voc['original_source']['original_host_butterfly_51_island_pair_sha256']
a=a[a.has_positive_2020_AIS_observed_from_to_pair.eq(1)].copy()
assert len(a)==ais['source_undirected_island_pairs_observed']==12
assert len(v)==voc['historical_17c_positive_source_island_pairs']==9
p=pd.read_parquet(R/'work/temporal/forecast/predictions_joint_island.parquet')
assert p.landmass_id.nunique()==25
held=set(p.landmass_id)
rows=[];summary=[]
for name,edges,period,admissible,reason in [('AIS2020',a,'2020',False,'Postdates2011–2014/2015–2018 predictor windows; overlaps response2019–2022'),('VOC17c',v,'1595–1698',False,'Predates predictors but is sparse archival history, not measured2011/2015 transport intensity')]:
    nodes=set(edges.island_a)|set(edges.island_b)
    overlap=nodes&held
    both=edges.island_a.isin(held)&edges.island_b.isin(held)
    q=p[p.landmass_id.isin(overlap)]
    summary.append(dict(source=name,period=period,positive_links=len(edges),positive_endpoint_islands=len(nodes),holdout_islands_with_positive_endpoint=len(overlap),links_both_endpoints_in_holdout=int(both.sum()),holdout_rows=len(q),holdout_species=q.species.nunique(),holdout_rows_excluding_GB_Ireland=int((~q.landmass_id.isin(['15','26'])).sum()),admissible_as_period_specific_shipping_predictor=admissible,reason=reason))
    for lid in sorted(held):
        sub=p[p.landmass_id.eq(lid)]
        rows.append(dict(source=name,landmass_id=lid,rows=len(sub),species=sub.species.nunique(),positive_route_endpoint=lid in nodes,interpretation='POSITIVE_ROUTE_DOCUMENTED' if lid in nodes else 'NO_POSITIVE_MATCH_UNKNOWN_NOT_NO_SHIPPING'))
    edges.assign(both_endpoints_in_holdout=both).to_csv(O/f'{name.lower()}_forecast_positive_route_overlap.csv',index=False)
pd.DataFrame(rows).to_csv(O/'voyage_forecast_island_coverage.csv',index=False)
pd.DataFrame(summary).to_csv(O/'voyage_forecast_overlap_summary.csv',index=False)
audit=dict(status='SOURCE_PINNED_SPATIAL_TEMPORAL_OVERLAP_AUDIT_NOT_NEW_MODEL',source_commit=sha,sources=sources,shared_original_pair_sha256=local_hash,summary=summary,limits=['Same original island-pair source hash verifies shared landmass ID lineage; not a new coastline reconstruction','Unmatched island/edge is unknown, never zero shipping','Historical voyage rows are not independent vessels or introductions','Neither source identifies ship-borne host establishment; no outcome-selected model fitted'])
(O/'voyage_forecast_overlap_audit.json').write_text(json.dumps(audit,indent=2))
lines=['# 実航路と島外予測対象の照合','', 'PR #117の固定コミット '+sha+' に保存された実航路の集計を、25島・95種の予測評価に結合した。元の51島ペア表のSHA256がローカルと両航路受領記録で完全一致することを確認して、陸塊IDの系統を照合した。','', '|資料|実航路の島数|25評価島との共通島|両端が評価島の航路|該当予測行|該当種|','|---|---:|---:|---:|---:|---:|']
for x in summary:lines.append(f"|{x['source']}|{x['positive_endpoint_islands']}|{x['holdout_islands_with_positive_endpoint']}|{x['links_both_endpoints_in_holdout']}|{x['holdout_rows']}|{x['holdout_species']}|")
lines+=['','2020年AISは現代の実航路の存在を示すが、2015–2018年説明変数より後で2019–2022年応答と重なる。したがって固定済みの予測モデルへ事前輸送量として追加しない。1595–1698年VOC史料は歴史的接続の記述に使えるが、現代の輸送強度の代用ではない。ここでの非採用は資料を否定するものではなく、固定予測期間・対象に対する判断である。','', '未照合の島・航路を無輸送とした対照群は作らない。表は実航路の端点との一致であって、食草の船舶持ち込み・チョウ到来の経路、直前の生息地を同定していない。現在の指数モデルと実航路モデルを同じ検証として数えない。','', '全225候補と155種の主解析は保持し、航路と一致する部分だけで主解析を置換しない。対応する植物種と年代の確認、期間整合した航路情報は別の未解決事項として残る。']
for x in summary:lines.append(f"\n{x['source']}端点と一致する評価行は英国・アイルランド除外後{x['holdout_rows_excluding_GB_Ireland']}行。多数の評価行を独立した島の数と解釈しない。")
(O/'voyage_forecast_overlap_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=True))

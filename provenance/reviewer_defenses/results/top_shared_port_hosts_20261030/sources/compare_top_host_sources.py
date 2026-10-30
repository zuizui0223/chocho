from pathlib import Path
import json
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs'
c=pd.read_csv(O/'top_shared_host_region_correlations.csv')
p=pd.read_csv(O/'top_shared_host_link_provenance.csv')
x=c.merge(p,on=['host_name','species'],validate='many_to_one')
assert len(x)==len(c) and x.in_frozen_HOSTS.notna().all()
rows=[]
for source,q in [('all_catalogue_links',x),('frozen_HOSTS_only',x[x.in_frozen_HOSTS])]:
    for (host,outcome),g in q.groupby(['host_name','outcome']):
        rows.append(dict(source=source,host_name=host,outcome=outcome,species=len(g),mean_rho=g.rho.mean(),positive_species=int(g.rho.gt(0).sum())))
s=pd.DataFrame(rows);s.to_csv(O/'top_shared_host_source_sensitivity.csv',index=False)
# Recompute retained species sets separately to check every summary.
allowed=set(map(tuple,p.loc[p.in_frozen_HOSTS,['host_name','species']].to_numpy()))
for z in s.itertuples():
    q=c[c.host_name.eq(z.host_name)&c.outcome.eq(z.outcome)]
    if z.source=='frozen_HOSTS_only':q=q[[tuple(v) in allowed for v in q[['host_name','species']].to_numpy()]]
    assert len(q)==z.species and np.isclose(q.rho.mean(),z.mean_rho)
    assert q.rho.gt(0).sum()==z.positive_species
receipt=dict(status='PASS',source_summary_rows=len(s),all_links=len(p),original_HOSTS_links=int(p.in_frozen_HOSTS.sum()),added_candidate_links=int(p.BCE_added_candidate.sum()),limits=['posthoc source sensitivity','does not independently validate feeding','same frozen outcomes and regional support; only link inclusion changes'])
(O/'top_shared_host_source_sensitivity_validation.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
table=s[s.outcome.eq('change')].pivot(index='host_name',columns='source',values=['species','mean_rho'])
lines=['# 共通上位食草：食草資料による感度比較','', '上位6種の49利用関係は42組が既存HOSTS表、7組がBCE追加候補だった。地域数等の条件を満たして相関を計算できた24組のうち、追加候補は4組。これらを外して相関の種間平均を再集計した。新しいモデルを学習した結果ではない。','', '|食草|全資料：種数／平均ρ|既存HOSTSのみ：種数／平均ρ|','|---|---:|---:|']
for name,z in table.iterrows():
    lines.append(f'|{name}|{int(z["species","all_catalogue_links"])} / {z["mean_rho","all_catalogue_links"]:+.3f}|{int(z["species","frozen_HOSTS_only"])} / {z["mean_rho","frozen_HOSTS_only"]:+.3f}|')
lines+=['','アルファルファは既存HOSTSのみでも平均がゼロ付近で、強い共通の正方向は現れなかった。Hesperis matronalisも同様。ナズナ・ヘラオオバコの正の相関と、ハリエンジュの負の相関は、今回のBCE追加候補の有無では変わらない。ただしこれらは各2種しか解析できていない。','', '既存HOSTS由来であることは、その地域での幼虫の利用や発育成功を確認したことを意味しない。植物は2001–2024年の集約記録であり、2015–2018年から2019–2022年のチョウ記録変化に先行した資源とは限らない。結果は記録された地域内の未補正の探索的関連である。']
(O/'top_shared_host_source_sensitivity_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(receipt))

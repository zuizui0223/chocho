from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
d=pd.read_csv(O/'island_host_network_pairs.csv',dtype={'island_a':str,'island_b':str})
i=pd.read_csv(O/'island_host_network_islands.csv',dtype={'landmass_id':str})
ids=sorted(i.loc[i.eligible,'landmass_id']);pos={lid:j for j,lid in enumerate(ids)}
assert len(d)==len(ids)*(len(ids)-1)//2 and not d.duplicated(['island_a','island_b']).any()
b=pd.read_csv(O/'global_original_species_island_history.csv',dtype={'landmass_id':str})
for metric,mask in [('new_butterfly_jaccard',b.first_record_year.between(2011,2022)),('earlier_butterfly_jaccard',b.first_record_year.le(2010))]:
    m=pd.crosstab(b.loc[mask,'landmass_id'],b.loc[mask,'speciesKey']).reindex(ids,fill_value=0).gt(0).astype(int).values
    shared=m@m.T; union=m.sum(axis=1)[:,None]+m.sum(axis=1)[None,:]-shared
    expected=np.array([shared[pos[a],pos[z]]/union[pos[a],pos[z]] if union[pos[a],pos[z]] else np.nan for a,z in zip(d.island_a,d.island_b)])
    np.testing.assert_allclose(d[metric],expected,equal_nan=True)
def partial(q):
    x=q.host_jaccard.rank().values;y=q.new_butterfly_jaccard.rank().values
    design=np.column_stack([np.ones(len(q)),q.port_distance_km.rank(),q.earlier_butterfly_jaccard.rank()])
    proj=np.linalg.pinv(design)
    return float(np.corrcoef(x-design@(proj@x),y-design@(proj@y))[0,1])
audit=json.loads((O/'island_host_network_audit.json').read_text())
q=d.dropna(subset=['new_butterfly_jaccard'])
np.testing.assert_allclose(partial(q),audit['partial_rank_host_new'])
leave=[]
for lid in ids:
    sub=q[(q.island_a!=lid)&(q.island_b!=lid)]
    leave.append(dict(omitted_island=lid,pairs=len(sub),rho=spearmanr(sub.host_jaccard,sub.new_butterfly_jaccard).statistic,partial_rank=partial(sub)))
leave=pd.DataFrame(leave);leave.to_csv(O/'island_host_network_leave_one_island.csv',index=False)
c=pd.read_csv(O/'island_host_network_all225.csv');assert len(c)==225 and c.species.nunique()==225
assert c.new_records_on_eligible_islands.sum()==i.loc[i.eligible,'new_butterfly_species'].sum()
v=dict(status='PASS',pairs_verified=len(d),butterfly_similarity_recomputed='binary matrix cross-product',
       partial_rank_recomputed='pandas ranks and pseudoinverse',
       leave_one_island_rho_range=[float(leave.rho.min()),float(leave.rho.max())],
       leave_one_island_partial_range=[float(leave.partial_rank.min()),float(leave.partial_rank.max())],
       caveat='Leave-one-island is post-result influence analysis, not a confidence interval, external validation, or causal adjustment. Host states rely on inventory classification audit.')
(O/'island_host_network_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8');print(json.dumps(v))

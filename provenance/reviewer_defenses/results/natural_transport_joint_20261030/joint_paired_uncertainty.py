from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,sys
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast';W=R/'work/barriers'
mode=sys.argv[1]
assert json.loads((O/f'forecast_joint_{mode}_results.json').read_text())['status']=='COMPLETE'
p=pd.read_parquet(F/f'predictions_joint_{mode}.parquet')
p=p.merge(pd.read_parquet(W/'cell_geography.parquet')[['cell','landmass_id','landmass_area_km2']],on='cell',validate='many_to_one')
p=p.merge(pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']],on='cell',how='left',validate='many_to_one')
island=p.strict_single.fillna(False)&p.landmass_area_km2.lt(7e6)
groups={'all155':p,'strict_islands':p[island],'strict_islands_without_GB_Ireland':p[island&~p.landmass_id.isin(['15','26'])]}
if mode=='spatial':groups['at_least_100km']=p[p.distance_to_training_cell_km.ge(100)]
rows=[]
for group,q in groups.items():
    species=sorted(q.species.unique());blocks=sorted(q.spatial_block.unique());rng=np.random.default_rng(20261030)
    count=pd.crosstab(q.species,q.spatial_block).reindex(index=species,columns=blocks,fill_value=0).to_numpy()
    sw=rng.multinomial(len(species),np.full(len(species),1/len(species)),size=2000)
    bw=rng.multinomial(len(blocks),np.full(len(blocks),1/len(blocks)),size=2000)
    denom=bw@count.T;weights=sw*(denom>0)
    assert (weights.sum(axis=1)>0).all()
    for ref in ['M5_host_composition','B5_sea_terrain','T1_port_connectivity']:
        delta=(q.J1_natural_transport-q.target)**2-(q[ref]-q.target)**2
        sums=q[['species','spatial_block']].assign(delta=delta).groupby(['species','spatial_block']).delta.sum().unstack(fill_value=0).reindex(index=species,columns=blocks,fill_value=0).to_numpy()
        np.testing.assert_allclose((sums.sum(axis=1)/count.sum(axis=1)).mean(),delta.groupby(q.species).mean().mean())
        means=np.divide(bw@sums.T,denom,out=np.zeros(denom.shape),where=denom>0)
        samples=(means*weights).sum(axis=1)/weights.sum(axis=1)
        lo,hi=np.quantile(samples,[.025,.975])
        rows.append(dict(mode=mode,cohort=group,reference=ref,mse_difference=delta.groupby(q.species).mean().mean(),interval_low=lo,interval_high=hi,species=len(species),blocks=len(blocks),rows=len(q),draws=2000))
out=pd.DataFrame(rows);out.to_csv(O/f'joint_{mode}_paired_uncertainty.csv',index=False)
print(out.to_string(index=False))

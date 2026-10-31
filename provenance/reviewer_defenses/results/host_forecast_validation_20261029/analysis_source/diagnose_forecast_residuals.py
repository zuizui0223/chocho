"""Post-hoc provider strata and spatial residual diagnostics, no model tuning."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy import sparse
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[2];S=R/'work/temporal/forecast';O=R/'outputs'
p=pd.read_parquet(S/'predictions_temporal.parquet')
t=pd.read_csv(O/'forecast_test_provider_turnover.csv')
p=p.merge(t,on='cell',validate='many_to_one')
assert p.provider_turnover.notna().all()
p['provider_stratum']=pd.cut(p.provider_turnover,[-1e-12,.25,.5,1.000001],labels=['low_0_025','middle_025_05','high_05_1'])
rows=[]
models=['persistence','M0_background','M3_both','M4_host_resources','M5_host_composition']
for stratum,q in p.groupby('provider_stratum',observed=True):
    for name in models:
        e=q[name]-q.target
        rows.append(dict(stratum=stratum,model=name,rows=len(q),species=q.species.nunique(),cells=q.cell.nunique(),mse=e.pow(2).groupby(q.species).mean().mean(),mae=e.abs().groupby(q.species).mean().mean()))
pd.DataFrame(rows).to_csv(O/'forecast_provider_stratified_metrics.csv',index=False)
coords=pd.read_parquet(S/'candidate_cells.parquet',columns=['cell','latitude','longitude'])
p=p.merge(coords,on='cell',validate='many_to_one')
rng=np.random.default_rng(20261029);diagnostics=[]
for species,q in p.groupby('species'):
    lat=np.radians(q.latitude.to_numpy());lon=np.radians(q.longitude.to_numpy())
    xyz=np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    dist,idx=cKDTree(xyz).query(xyz,k=min(9,len(q)),distance_upper_bound=2*np.sin(100/(2*6371.0088)))
    if dist.ndim==1:continue
    ii=np.repeat(np.arange(len(q)),idx.shape[1]);jj=idx.ravel();valid=(ii!=jj)&(jj<len(q))&np.isfinite(dist.ravel())
    matrix=sparse.csr_matrix((np.ones(valid.sum()),(ii[valid],jj[valid])),shape=(len(q),len(q)))
    matrix=matrix.maximum(matrix.T)
    active=np.asarray(matrix.sum(axis=1)).ravel()>0
    matrix=matrix[active][:,active]
    if active.sum()<30:continue
    matrix=sparse.diags(1/np.asarray(matrix.sum(axis=1)).ravel())@matrix
    for name in ['M0_background','M5_host_composition']:
        residual=(q.target-q[name]).to_numpy()[active];z=residual-residual.mean();den=z@z
        if den==0:continue
        value=float(z@(matrix@z)/den)
        null=np.array([perm@(matrix@perm)/den for perm in (rng.permutation(z) for _ in range(499))])
        # One-sided positive autocorrelation test, declared diagnostic.
        pv=(1+np.count_nonzero(null>=value))/500
        diagnostics.append(dict(species=species,model=name,cells=int(active.sum()),excluded_isolated_cells=int((~active).sum()),moran_i=value,p_positive=pv))
d=pd.DataFrame(diagnostics)
for name,inds in d.groupby('model').groups.items():d.loc[inds,'p_fdr']=multipletests(d.loc[inds,'p_positive'],method='fdr_bh')[1]
d.to_csv(O/'forecast_spatial_residual_diagnostics.csv',index=False)
summary={'status':'COMPLETE','provider_strata':'Post-hoc total-variation cutoffs <=0.25, <=0.5, >0.5; future provider composition used for descriptive scoring only',
    'spatial_graph':'Per species: symmetrized 8-nearest-neighbour graph, great-circle distance <=100km; drop isolated cells; row standardize; at least30 connected cells',
    'permutations':499,'test':'One-sided positive Moran autocorrelation; BH across tested species separately per model; residual-permutation diagnostic, not proof of absence of source bias',
    'models':{name:{'tested_species':len(q),'median_moran_i':float(q.moran_i.median()),'positive_fdr05_species':int(((q.moran_i>0)&(q.p_fdr<.05)).sum())} for name,q in d.groupby('model')}}
(O/'forecast_residual_diagnostics.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

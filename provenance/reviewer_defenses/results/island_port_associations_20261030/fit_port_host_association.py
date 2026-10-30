from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups
R=Path(__file__).resolve().parents[2];O=R/'outputs'
p=pd.read_parquet(R/'work/temporal/forecast/species_cell_period_panel.parquet')
a=p[p.period.eq(1)].copy();b=p[p.period.eq(2)][['speciesKey','cell','recording_probability10','days','provider_dominance']].rename(columns={'recording_probability10':'target','days':'future_days','provider_dominance':'future_provider_dominance'})
a=a.merge(b,on=['speciesKey','cell'],validate='one_to_one')
g=pd.read_parquet(R/'work/barriers/cell_geography.parquet');a=a.merge(g[['cell','landmass_id','landmass_area_km2','distance_port_2019_km']],on='cell',validate='many_to_one')
a=a[a.landmass_area_km2.lt(7e6)].copy()
a['log_port_distance']=np.log1p(a.distance_port_2019_km);a['introduced_centered']=a.introduced_share-.5
a['port_x_introduced']=a.log_port_distance*a.introduced_centered
a['log_days']=np.log1p(a.days);a['log_future_days']=np.log1p(a.future_days)
features=['log_port_distance','introduced_centered','port_x_introduced','recording_probability10','log_days','log_future_days','provider_dominance','future_provider_dominance','temperature','precipitation','temperature_trend_per_year','precipitation_trend_per_year','cropland','built','tree','cropland_trend_per_year','built_trend_per_year','tree_trend_per_year']
rows=[];audits=[]
for name,mask in [('all_islands',np.ones(len(a),bool)),('excluding_GB_Ireland',~a.landmass_id.isin(['15','26'])),('stable_provider',abs(a.provider_dominance-a.future_provider_dominance)<=.1)]:
 q=a.loc[mask].dropna(subset=features+['target']).copy();before=len(q)
 q['stratum']=q.speciesKey.astype(str)+'_'+q.landmass_id.astype(str)
 q=q[q.groupby('stratum').cell.transform('size')>=2].copy()
 weights=q.speciesKey.map(1/q.speciesKey.value_counts()).to_numpy();weights*=len(weights)/weights.sum()
 v=q[features+['target']].astype(float)
 # Constant weights within species-landmass: ordinary stratum centering equals weighted centering.
 z=v-v.groupby(q.stratum).transform('mean')
 kept=[f for f in features if z[f].std()>1e-10];scales=z[kept].std();X=(z[kept]/scales).to_numpy();y=z.target.to_numpy()
 rank=np.linalg.matrix_rank(X);assert rank==len(kept)
 fit=sm.WLS(y,X,weights=weights).fit()
 sp=pd.factorize(q.speciesKey)[0];block=pd.factorize(q.spatial_block)[0]
 cov,_,_=cov_cluster_2groups(fit,sp,block,use_correction=True)
 # Account for absorbed stratum intercepts in finite-sample residual degrees of freedom.
 n=len(q);k=len(kept);groups=q.stratum.nunique();cov*= (n-k)/(n-k-groups)
 diagonal=np.diag(cov);se=np.sqrt(np.maximum(diagonal,0))/scales.to_numpy();beta=fit.params/scales.to_numpy()
 for i,f in enumerate(kept):
  rows.append(dict(subset=name,feature=f,coefficient=beta[i],cluster_se=se[i] if diagonal[i]>=0 else np.nan,ci_low=beta[i]-1.96*se[i] if diagonal[i]>=0 else np.nan,ci_high=beta[i]+1.96*se[i] if diagonal[i]>=0 else np.nan))
 # Independent coefficient solution and within-stratum centering checks.
 direct=np.linalg.lstsq(X*np.sqrt(weights[:,None]),y*np.sqrt(weights),rcond=None)[0]
 assert np.max(abs(direct-fit.params))<1e-9
 assert z.groupby(q.stratum).mean().abs().max().max()<1e-8
 audits.append(dict(subset=name,complete_rows_before_singleton_removal=before,rows=len(q),species=int(q.speciesKey.nunique()),landmasses=int(q.landmass_id.nunique()),cells=int(q.cell.nunique()),spatial_blocks=int(q.spatial_block.nunique()),strata=int(groups),rank=int(rank),condition=float(np.linalg.cond(X)),dropped_invariant_features=[f for f in features if f not in kept],negative_covariance_diagonal=int((diagonal<0).sum()),covariance_min_eigenvalue=float(np.linalg.eigvalsh(cov).min()),coefficient_lstsq_max_difference=float(np.max(abs(direct-fit.params)))))
pd.DataFrame(rows).to_csv(O/'port_host_adjusted_associations.csv',index=False)
(O/'port_host_association_audit.json').write_text(json.dumps({'status':'EXPLORATORY_FIT','fits':audits,'limits':['two-way clustering does not guarantee removal of all spatial dependence','normal approximation with limited geographic clusters','static regional hosts do not measure local resource supply','2019 ports contemporaneous; no historical forecast or causal attribution']},indent=2),encoding='utf-8')
print(pd.DataFrame(rows).query("feature in ['log_port_distance','port_x_introduced']").to_string(index=False))
print(json.dumps(audits))

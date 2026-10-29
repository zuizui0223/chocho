from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,warnings
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[2];S=R/'work/temporal/forecast';O=R/'outputs'
env=['temperature','cropland','built','precipitation','tree']
hosts=['global_hosts','available_share','host_family_Poaceae']
cols=['species','region','cell','period','recording_probability10','days','provider_dominance']+env+hosts
p=pd.read_parquet(S/'species_cell_period_panel.parquet',columns=cols)
tables=[]
for period in [0,1]:
    a=p[p.period.eq(period)].merge(p[p.period.eq(period+1)],on=['species','region','cell'],suffixes=('_before','_after'),validate='one_to_one')
    a=a[a.recording_probability10_before.notna()&a.recording_probability10_after.notna()].copy()
    a['change']=a.recording_probability10_after-a.recording_probability10_before
    a['baseline']=a.recording_probability10_before
    for c in env:a['delta_'+c]=a[c+'_after']-a[c+'_before']
    a['effort_change']=np.log1p(a.days_after)-np.log1p(a.days_before)
    a['dominance_change']=a.provider_dominance_after-a.provider_dominance_before
    a['breadth']=np.log1p(a.global_hosts_before)
    a['availability']=a.available_share_before;a['grass_host_share']=a.host_family_Poaceae_before
    variables=['change','baseline','effort_change','dominance_change','breadth','availability','grass_host_share']+['delta_'+c for c in env]
    g=a.groupby(['species','region'])[variables].mean().join(a.groupby(['species','region']).size().rename('cells')).reset_index()
    g['transition']=period;tables.append(g[g.cells>=10].copy())
train,test=tables
predictors=[x for x in variables if x!='change']
scale={}
for c in predictors:
    mean=float(train[c].mean());sd=float(train[c].std());assert sd>0
    scale[c]={'mean':mean,'sd':sd}
    for t in tables:t['z_'+c]=(t[c]-mean)/sd
base='change ~ '+' + '.join('z_'+c for c in predictors)+' + C(region)'
interactions=[f'z_delta_{e}:z_{h}' for e in ['temperature','cropland','built'] for h in ['breadth','availability','grass_host_share']]
results={};allcoef=[]
for name,formula in [('additive',base),('interactions',base+' + '+' + '.join(interactions))]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        model=smf.mixedlm(formula,train,groups=train.species).fit(reml=False,method='powell',maxiter=200,disp=False)
    with warnings.catch_warnings(record=True) as check_warnings:
        warnings.simplefilter('always')
        alternate=smf.mixedlm(formula,train,groups=train.species).fit(reml=False,method='lbfgs',start_params=model.params,maxiter=200,disp=False)
    coef=pd.DataFrame({'term':model.params.index,'estimate':model.params.values,'se':model.bse.values,'p':model.pvalues.values,'ci_low':model.conf_int()[0].values,'ci_high':model.conf_int()[1].values})
    coef['model']=name;allcoef.append(coef)
    eligible=test[test.region.isin(train.region)&test.species.isin(train.species)].copy()
    prediction=model.predict(eligible)
    # Include learned species random intercepts only for species observed in training.
    try:prediction+=eligible.species.map({k:float(v.iloc[0]) for k,v in model.random_effects.items()}).to_numpy()
    except ValueError:pass
    eligible['prediction']=prediction;eligible.to_csv(O/f'host_interaction_{name}_hindcast.csv',index=False)
    results[name]={'converged':bool(model.converged),'loglikelihood':float(model.llf),'aic':float(model.aic),'random_intercept_variance':float(model.cov_re.iloc[0,0]),'warnings':[str(x.message) for x in w], 'test_rows':len(eligible),'mse':float(np.mean((prediction-eligible.change)**2)),'mae':float(np.mean(abs(prediction-eligible.change)))}
    results[name].update(residual_variance=float(model.scale),alternate_optimizer_converged=bool(alternate.converged),alternate_loglikelihood=float(alternate.llf),optimizer_max_fixed_effect_difference=float(abs(model.fe_params-alternate.fe_params).max()),alternate_warnings=[str(x.message) for x in check_warnings])
    print(name,results[name],flush=True)
    if name=='interactions':
        selected=coef[coef.term.isin(interactions)].copy()
        selected['p_holm']=multipletests(selected.p,method='holm')[1]
        selected.to_csv(O/'host_environment_interactions.csv',index=False)
pd.concat(allcoef).to_csv(O/'host_interaction_all_coefficients.csv',index=False)
pd.concat(tables).to_parquet(S/'host_interaction_aggregates.parquet',index=False)
receipt={'status':'COMPLETE' if all(x['converged'] for x in results.values()) else 'CONVERGENCE_WARNING','train_rows':len(train),'train_species':int(train.species.nunique()),'train_regions':int(train.region.nunique()),'test_rows_before_category_filter':len(test),'results':results,'scaling':scale,'interpretation':'Exploratory aggregate mixed-model associations using realised subsequent environment; not strict forecast or causal effects'}
(O/'host_interaction_results.json').write_text(json.dumps(receipt,indent=2))

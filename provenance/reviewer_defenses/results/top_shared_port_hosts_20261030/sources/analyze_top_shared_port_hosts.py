from pathlib import Path
import json
import pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
lock=json.loads((O/'top_shared_host_region_analysis_lock.json').read_text())
ports=pd.read_csv(O/'port_flora_regional_support.csv');ports=ports[ports.region_count.eq(1)].rename(columns={'candidate_regions':'region'})
denom=ports.groupby('region').port_id.nunique()
hosts=pd.read_parquet(W/'port_named_host_regional_status.parquet')
hosts=hosts[hosts.regional_catalogue_status.eq('contemporary_not_native_in_catalogue')].merge(ports[['port_id','region']],on='port_id',validate='many_to_one')
edges=pd.read_parquet(W/'all225_named_host_edges_local.parquet')
p=pd.read_parquet(R/'work/temporal/forecast/predictions_temporal.parquet')
assert not p.duplicated(['species','cell']).any()
p['change']=p.target-p.baseline_probability
p['residual_M3']=p.target-p.M3_both;p['residual_M5']=p.target-p.M5_host_composition
regional=p.groupby(['species','region']).agg(cells=('cell','size'),change=('change','mean'),residual_M3=('residual_M3','mean'),residual_M5=('residual_M5','mean')).reset_index()
regional=regional[regional.cells.ge(lock['minimum_cells_per_species_region'])]
rows=[];checks=[];linked=[]
for host in lock['hosts']:
    count=hosts[hosts.host_name.eq(host)].groupby('region').port_id.nunique()
    exposure=pd.DataFrame({'host_ports':count,'complete_ports':denom.reindex(count.index)})
    exposure['port_fraction']=exposure.host_ports/exposure.complete_ports
    for sp in sorted(set(edges.loc[edges.host_name.eq(host),'species'])):
        q=regional[regional.species.eq(sp)].merge(exposure.reset_index(),on='region',validate='one_to_one')
        q['host_name']=host;linked.append(q)
        for outcome in ['change','residual_M3','residual_M5']:
            eligible=len(q)>=lock['minimum_regions_per_species_correlation'] and q.port_fraction.nunique()>1 and q[outcome].nunique()>1
            checks.append(dict(host_name=host,species=sp,outcome=outcome,regions=len(q),status='ELIGIBLE' if eligible else 'INSUFFICIENT_REGIONS_OR_VARIATION'))
            if eligible:rows.append(dict(host_name=host,species=sp,outcome=outcome,regions=len(q),rho=float(spearmanr(q.port_fraction,q[outcome]).statistic)))
c=pd.DataFrame(rows);pd.DataFrame(checks).to_csv(O/'top_shared_host_region_support.csv',index=False)
pd.concat(linked,ignore_index=True).to_csv(O/'top_shared_host_region_join.csv',index=False)
c.to_csv(O/'top_shared_host_region_correlations.csv',index=False)
s=c.groupby(['host_name','outcome']).agg(species=('species','size'),mean_rho=('rho','mean'),median_rho=('rho','median'),positive_species=('rho',lambda x:int(x.gt(0).sum()))).reset_index()
s.to_csv(O/'top_shared_host_region_summary.csv',index=False)
print(s.to_string(index=False))

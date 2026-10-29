from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
from scipy.stats import hypergeom
from sklearn.neighbors import BallTree
from pyproj import Transformer
R=Path(__file__).resolve().parents[2];W=R/'work/temporal';S=W/'forecast';O=R/'outputs'
lock=json.loads((O/'forecast_model_lock.json').read_text());assert lock['status']=='LOCKED_BEFORE_MODEL_FIT'
cells=pd.read_parquet(S/'candidate_cells.parquet');env=pd.read_parquet(S/'cell_environment_period.parquet')
base=env[env.period==0].set_index('cell')
cells=cells[(cells.center_region_matches==1)&~cells.region_boundary_flag&cells.cell.map(base.unclassified).le(.05)].copy()
grid=pd.read_parquet(S/'cell_climate_grid.parquet');cells=cells.merge(grid[['cell','climate_id','climate_lat','climate_lon']],on='cell',validate='one_to_one')
cells['spatial_block']=np.floor(cells.climate_lat/5).astype(int).astype(str)+'_'+np.floor(cells.climate_lon/5).astype(int).astype(str)
disp=pd.read_csv(O/'forecast_candidate_disposition.csv');eligible=disp[disp.stage_disposition.str.startswith('BASELINE_ELIGIBLE')]
days=pd.read_parquet(S/'species_recording_days.parquet');days['period']=((days.year-2011)//4).astype(int)
eff=pd.read_parquet(S/'target_group_recording_days.parquet');eff['period']=((eff.year-2011)//4).astype(int)
denom=eff.groupby(['cell','period']).size().rename('days').reset_index()
numer=days.groupby(['speciesKey','cell','period']).size().rename('success').reset_index()
quarters=eff.assign(quarter=(eff.month.astype(int)-1)//3).groupby(['cell','period','quarter']).size().unstack(fill_value=0).reindex(columns=range(4),fill_value=0)
quarters=quarters.div(quarters.sum(axis=1),axis=0);quarters.columns=['season_q'+str(i) for i in range(4)];quarters=quarters.reset_index()
provider=eff.groupby(['cell','period','datasetKey']).size().rename('n').reset_index();dominance=(provider.groupby(['cell','period']).n.max()/provider.groupby(['cell','period']).n.sum()).rename('provider_dominance').reset_index()
env=env.merge(denom,on=['cell','period'],how='left',validate='one_to_one').merge(quarters,on=['cell','period'],how='left',validate='one_to_one').merge(dominance,on=['cell','period'],how='left',validate='one_to_one')
past=pd.read_parquet(W/'species_year_cells_10km.parquet',columns=['speciesKey','year','cell_x_10km','cell_y_10km']);past=past[past.year<=2014].drop_duplicates(['speciesKey','cell_x_10km','cell_y_10km'])
inv=Transformer.from_crs(6933,4326,always_xy=True)
target_radians=np.radians(cells[['latitude','longitude']].values)
rows=[]
for sp in eligible.itertuples():
 g=past[past.speciesKey==sp.taxon_key];lon,lat=inv.transform((g.cell_x_10km.values+.5)*10000,(g.cell_y_10km.values+.5)*10000)
 tree=BallTree(np.radians(np.column_stack([lat,lon])),metric='haversine');distance=tree.query(target_radians,k=1)[0].ravel()*6371.0088
 chosen=cells[distance<=500].copy();chosen['distance_km']=distance[distance<=500];chosen['speciesKey']=sp.taxon_key;chosen['species']=sp.species;rows.append(chosen)
domain=pd.concat(rows,ignore_index=True)
metrics=pd.read_parquet(W/'host_region_composition_metrics.parquet').drop(columns=['composition_hash'])
counts=pd.read_csv(R/'work/remote_sensing/host_region_counts.csv')[['species','region','native_host_count','introduced_host_count','introduced_share']]
metrics=metrics.merge(counts,on=['species','region'],validate='one_to_one');metrics.introduced_share=metrics.introduced_share.fillna(0)
family=pd.read_parquet(W/'host_region_family_composition.parquet');fs=pd.read_csv(O/'temporal_host_family_support.csv');chosen_families=fs[(fs.butterfly_species>=10)&(fs.regions>=10)].family.tolist()
fp=family[family.family.isin(chosen_families)].pivot(index=['species','region'],columns='family',values='share').fillna(0);fp.columns=['host_family_'+x for x in fp.columns]
metrics=metrics.merge(fp.reset_index(),on=['species','region'],how='left',validate='one_to_one').fillna({c:0 for c in fp.columns})
domain=domain.merge(metrics,on=['species','region'],how='left',validate='many_to_one');assert domain.global_hosts.notna().all()
panels=[]
for period in [0,1,2]:
 p=domain.copy();p['period']=period;p=p.merge(env,on=['cell','period'],how='left',validate='many_to_one').merge(numer,on=['speciesKey','cell','period'],how='left',validate='one_to_one')
 p.success=p.success.fillna(0);p['days']=p.days.fillna(0)
 assert (p.success<=p.days).all()
 p['recording_probability10']=np.nan;ok=p.days>=10
 p.loc[ok,'recording_probability10']=hypergeom.sf(0,p.loc[ok,'days'].astype(int),p.loc[ok,'success'].astype(int),10)
 assert p.loc[ok,'recording_probability10'].between(0,1).all()
 panels.append(p)
panel=pd.concat(panels,ignore_index=True);panel.to_parquet(S/'species_cell_period_panel.parquet',index=False,compression='zstd')
manifest=dict(status='PASS',candidate_species=155,primary_species=int(domain.species.nunique()),primary_cells=int(domain.cell.nunique()),species_cell_pairs=len(domain),period_rows=len(panel),host_family_features=chosen_families,model_lock_sha256=hashlib.sha256((O/'forecast_model_lock.json').read_bytes()).hexdigest(),response='Rarefied recording probability in10 provider-days, not true occupancy',holdout_labels='Built mechanically but no model fit/score/tuning yet')
(O/'forecast_panel_audit.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2),flush=True)

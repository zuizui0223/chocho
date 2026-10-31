"""Rebuild responses and effort covariates, holding the primary domain fixed.
Sensitivity definitions are post-primary-analysis robustness checks, not tuning.
"""
from pathlib import Path
import json,sqlite3
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
R=Path(__file__).resolve().parents[2]; W=R/'work/temporal'; S=W/'forecast'; O=R/'outputs'
primary=S/'species_cell_period_panel.parquet'
top=json.loads((O/'forecast_recording_day_audit.json').read_text())['top_provider']
audit={'status':'IN_PROGRESS','definitions':'Post-primary robustness; identical fixed candidate species/domain; no holdout tuning', 'variants':{}}
def save_panel(label,denom=None,numer=None,quarter=None,dominance=None,n=10,radius=None):
    p=pd.read_parquet(primary)
    if radius is not None:p=p[p.distance_km<=radius].copy()
    if denom is not None:
        p=p.drop(columns=['days','success','provider_dominance']+['season_q'+str(i) for i in range(4)])
        for table,keys in [(denom,['cell','period']),(numer,['speciesKey','cell','period']),(quarter,['cell','period']),(dominance,['cell','period'])]:
            p=p.merge(table,on=keys,how='left',validate='many_to_one')
        p[['days','success']]=p[['days','success']].fillna(0)
    assert (p.success<=p.days).all()
    ok=p.days>=n;p['recording_probability10']=np.nan
    p.loc[ok,'recording_probability10']=hypergeom.sf(0,p.loc[ok,'days'].astype(int),p.loc[ok,'success'].astype(int),n)
    assert p.loc[ok,'recording_probability10'].between(0,1).all()
    p.to_parquet(S/f'panel_{label}.parquet',index=False,compression='zstd')
    audit['variants'][label]={'rows':len(p),'rarefaction_n':n,'supported_period_rows':int(ok.sum()),'species':int(p.species.nunique()),'fixed_domain':True}
    (O/'forecast_sensitivity_preparation.json').write_text(json.dumps(audit,indent=2))
    print('PREPARED',label,audit['variants'][label],flush=True)

save_panel('effort20',n=20)
save_panel('radius250',radius=250)
d=pd.read_parquet(S/'species_recording_days.parquet')
d=d[d.datasetKey!=top].copy();d['period']=((d.year-2011)//4).astype(int)
e=d[['cell','period','datasetKey','date','month']].drop_duplicates(['cell','datasetKey','date'])
den=e.groupby(['cell','period']).size().rename('days').reset_index()
num=d.groupby(['speciesKey','cell','period']).size().rename('success').reset_index()
q=e.assign(quarter=(e.month.astype(int)-1)//3).groupby(['cell','period','quarter']).size().unstack(fill_value=0).reindex(columns=range(4),fill_value=0)
q=q.div(q.sum(axis=1),axis=0);q.columns=['season_q'+str(i) for i in range(4)];q=q.reset_index()
g=e.groupby(['cell','period','datasetKey']).size().rename('n').reset_index()
dom=(g.groupby(['cell','period']).n.max()/g.groupby(['cell','period']).n.sum()).rename('provider_dominance').reset_index()
del d,e,g
save_panel('exclude_top_provider',den,num,q,dom)
audit['excluded_provider']=top

# SQLite bounds memory and ensures duplicates spanning source partitions collapse.
db=S/'precision100_recording_days.sqlite'
con=sqlite3.connect(db)
con.execute('CREATE TABLE IF NOT EXISTS units (speciesKey INTEGER,cell TEXT,provider TEXT,date TEXT,period INTEGER,quarter INTEGER,PRIMARY KEY(speciesKey,cell,provider,date)) WITHOUT ROWID')
con.execute('DELETE FROM units');con.commit()
cells=pd.read_parquet(S/'candidate_cells.parquet',columns=['cell_x_10km','cell_y_10km'])
allowed=set(zip(cells.cell_x_10km,cells.cell_y_10km))
source=json.loads((O/'temporal_source_audit.json').read_text())
for i,item in enumerate(source['parts']):
    d=pd.read_parquet(W/'occurrence_parts'/item['file'],columns=['speciesKey','year','month','day','datasetKey','cell_x_10km','cell_y_10km','coordinateUncertaintyInMeters'])
    d=d[d.year.between(2011,2022)&d.coordinateUncertaintyInMeters.le(100)].copy()
    d=d.loc[np.array([xy in allowed for xy in zip(d.cell_x_10km,d.cell_y_10km)],dtype=bool)].copy()
    if d.empty:continue
    dates=pd.to_datetime(d[['year','month','day']].apply(pd.to_numeric,errors='coerce').astype('float64'),errors='coerce')
    d=d[dates.notna()].copy();d['date']=dates[dates.notna()].dt.strftime('%Y-%m-%d')
    d['cell']=d.cell_x_10km.astype(str)+':'+d.cell_y_10km.astype(str)
    d['period']=((d.year-2011)//4).astype(int);d['quarter']=((d.month-1)//3).astype(int)
    values=((int(a),str(b),str(c),str(dt),int(e),int(f)) for a,b,c,dt,e,f in d[['speciesKey','cell','datasetKey','date','period','quarter']].drop_duplicates().itertuples(index=False,name=None))
    con.executemany('INSERT OR IGNORE INTO units VALUES (?,?,?,?,?,?)',values);con.commit()
    if i%10==0:print('PRECISION_SOURCE_PART',i+1,flush=True)
con.execute('DROP TABLE IF EXISTS effort')
con.execute('CREATE TABLE effort AS SELECT DISTINCT cell,provider,date,period,quarter FROM units');con.commit()
den=pd.read_sql_query('SELECT cell,period,count(*) AS days FROM effort GROUP BY cell,period',con)
num=pd.read_sql_query('SELECT speciesKey,cell,period,count(*) AS success FROM units GROUP BY speciesKey,cell,period',con)
q=pd.read_sql_query('SELECT cell,period,quarter,count(*) AS n FROM effort GROUP BY cell,period,quarter',con).pivot(index=['cell','period'],columns='quarter',values='n').fillna(0).reindex(columns=range(4),fill_value=0)
q=q.div(q.sum(axis=1),axis=0);q.columns=['season_q'+str(i) for i in range(4)];q=q.reset_index()
g=pd.read_sql_query('SELECT cell,period,provider,count(*) AS n FROM effort GROUP BY cell,period,provider',con)
dom=(g.groupby(['cell','period']).n.max()/g.groupby(['cell','period']).n.sum()).rename('provider_dominance').reset_index()
audit['precision100_unique_species_days']=int(con.execute('SELECT count(*) FROM units').fetchone()[0]);con.close()
save_panel('precision100',den,num,q,dom)
audit['status']='COMPLETE';(O/'forecast_sensitivity_preparation.json').write_text(json.dumps(audit,indent=2))
print('COMPLETE',flush=True)

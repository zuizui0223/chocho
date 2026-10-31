from pathlib import Path
import json,math,sqlite3
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[2];S=R/'work/temporal/forecast';O=R/'outputs'
cols=['speciesKey','species','cell','period','days','success','recording_probability10','distance_km']
base=pd.read_parquet(S/'species_cell_period_panel.parquet',columns=cols).set_index(['speciesKey','cell','period'])
disp=pd.read_csv(O/'forecast_candidate_disposition.csv');audit={}
for label in ['effort20','radius250','exclude_top_provider','precision100']:
    d=pd.read_parquet(S/f'panel_{label}.parquet',columns=cols).set_index(['speciesKey','cell','period'])
    assert d.index.is_unique and d.index.isin(base.index).all()
    b=base.loc[d.index]
    assert (d.days<=b.days).all() and (d.success<=b.success).all()
    n=20 if label=='effort20' else 10
    assert np.array_equal(d.recording_probability10.notna(),d.days>=n)
    if label=='radius250':
        pd.testing.assert_frame_equal(d.sort_index(),base[base.distance_km<=250].sort_index())
    if label=='effort20':
        assert (d.days==b.days).all() and (d.success==b.success).all()
        ok=d.recording_probability10.notna()
        assert (d.loc[ok,'recording_probability10']>=b.loc[ok,'recording_probability10']-1e-12).all()
    samples=d[d.days>=n].sample(100,random_state=20261029)
    for row in samples.itertuples():
        N,K=int(row.days),int(row.success)
        expected=1-(math.comb(N-K,n)/math.comb(N,n) if N-K>=n else 0)
        assert abs(expected-row.recording_probability10)<1e-12
    if label=='precision100':
        con=sqlite3.connect(S/'precision100_recording_days.sqlite')
        for (sp,cell,period),row in samples.iterrows():
            success=con.execute('SELECT count(*) FROM units WHERE speciesKey=? AND cell=? AND period=?',(int(sp),cell,int(period))).fetchone()[0]
            days=con.execute('SELECT count(*) FROM effort WHERE cell=? AND period=?',(cell,int(period))).fetchone()[0]
            assert success==row.success and days==row.days
        con.close()
    d=d.reset_index();counts={}
    for before,prefix in [(0,'train'),(1,'test')]:
        a=d[d.period.eq(before)&d.recording_probability10.notna()]
        b=d[d.period.eq(before+1)&d.recording_probability10.notna()]
        q=a.merge(b,on=['speciesKey','species','cell'],validate='one_to_one')
        group=q.groupby('species').size()
        disp[label+'_'+prefix+'_pairs']=disp.species.map(group).fillna(0).astype(int)
        counts[prefix+'_rows']=len(q);counts[prefix+'_species']=q.species.nunique()
    audit[label]={'status':'PASS','independent_combination_formula_checks':100,'source_sql_count_checks':100 if label=='precision100' else 0,**counts}
disp.to_csv(O/'forecast_sensitivity_all_225_support.csv',index=False)
(O/'forecast_sensitivity_validation.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))

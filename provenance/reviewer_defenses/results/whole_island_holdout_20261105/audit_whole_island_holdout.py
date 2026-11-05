from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers';F=R/'work/temporal/forecast'
p=pd.read_parquet(F/'species_cell_period_panel.parquet',columns=['speciesKey','species','cell','period','recording_probability10','spatial_block'])
def transition(k):
    a=p[p.period.eq(k)&p.recording_probability10.notna()]
    b=p[p.period.eq(k+1)&p.recording_probability10.notna()][['speciesKey','cell']]
    return a.merge(b,on=['speciesKey','cell'],validate='one_to_one')
tr=transition(0);te=transition(1)
geo=pd.read_parquet(W/'cell_geography.parquet')[['cell','landmass_id','landmass_area_km2']]
strict=pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']]
area=pd.read_csv(O/'cell_land_overlap.csv')[['cell','area_supported']]
parts=pd.read_parquet(W/'cell_land_overlap_parts.parquet')
te=te.merge(geo,on='cell',validate='many_to_one').merge(strict,on='cell',how='left',validate='many_to_one').merge(area,on='cell',validate='many_to_one')
island=te.strict_single.fillna(False)&te.area_supported&te.landmass_area_km2.lt(7e6)
test=te[island].copy()
test['island_fold']=test.landmass_id.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)
rows=[];folds=[]
for lid,q in test.groupby('landmass_id'):
    touched=set(parts.loc[parts.landmass_id.eq(lid),'cell'])
    training=tr[tr.cell.isin(touched)]
    rows.append(dict(landmass_id=lid,fold=int(q.island_fold.iloc[0]),test_rows=len(q),test_cells=q.cell.nunique(),test_species=q.species.nunique(),
        train_rows_on_intersecting_cells=len(training),train_cells=training.cell.nunique(),train_spatial_blocks=training.spatial_block.nunique(),
        test_spatial_blocks=q.spatial_block.nunique()))
for f in range(5):
    q=test[test.island_fold.eq(f)];held=set(q.landmass_id)
    touched=set(parts.loc[parts.landmass_id.isin(held),'cell'])
    train=tr[~tr.cell.isin(touched)]
    assert not set(train.cell)&touched
    folds.append(dict(fold=f,test_islands=len(held),test_rows=len(q),test_species=q.species.nunique(),
        train_rows=len(train),train_species=train.species.nunique(),excluded_train_rows=len(tr)-len(train)))
pd.DataFrame(rows).to_csv(O/'whole_island_holdout_support.csv',index=False)
pd.DataFrame(folds).to_csv(O/'whole_island_holdout_fold_support.csv',index=False)
test[['speciesKey','cell','landmass_id','island_fold']].to_parquet(W/'whole_island_holdout_test_keys.parquet',index=False)
audit=dict(status='HOLDOUT_SUPPORT_ONLY_NOT_MODEL_RESULT',primary_candidates=225,primary_train_species=tr.species.nunique(),
    secondary_test_species=test.species.nunique(),test_islands=test.landmass_id.nunique(),test_rows=len(test),
    test_rows_GB_Ireland=int(test.landmass_id.isin(['15','26']).sum()),nonempty_folds=sum(x['test_rows']>0 for x in folds),
    limits='Secondary island subset only. Exclude all training cells with any area intersection with held islands; historical predictor records remain available up to origin. Not an unrecorded-island forecast or island-independent occurrence sampling.')
(O/'whole_island_holdout_support_audit.json').write_text(json.dumps(audit,indent=2))
print(pd.DataFrame(folds).to_string(index=False));print(json.dumps(audit))

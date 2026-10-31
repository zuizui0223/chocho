from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers';F=R/'work/temporal/forecast'
area=pd.read_csv(O/'cell_land_overlap.csv',dtype={'centre_landmass_id':str,'dominant_landmass_id':str})
assert len(area)==5282 and area.cell.is_unique and area.area_supported.isin([True,False]).all()
parts=pd.read_parquet(W/'cell_land_overlap_parts.parquet')
totals=parts.groupby('cell').intersection_km2.sum()
np.testing.assert_allclose(area.cell.map(totals).fillna(0)/100,area.land_fraction,atol=1e-8)
geo=pd.read_parquet(W/'cell_geography.parquet')[['cell','landmass_id','landmass_area_km2']]
strict=pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']]
rows=[]
coverage=pd.read_csv(O/'forecast_all_225_final_coverage.csv')[['species','taxon_key','stage_disposition']]
assert len(coverage)==225 and coverage.taxon_key.is_unique
for mode in ['temporal','spatial','temporal_effort20']:
    assert json.loads((O/f'forecast_joint_{mode}_results.json').read_text())['status']=='COMPLETE'
    p=pd.read_parquet(F/f'predictions_joint_{mode}.parquet')
    p=p.merge(geo,on='cell',validate='many_to_one').merge(strict,on='cell',how='left',validate='many_to_one').merge(area[['cell','area_supported']],on='cell',validate='many_to_one')
    island=p.strict_single.fillna(False)&p.landmass_area_km2.lt(7e6)
    groups={'strict_original':island,'strict_area_supported':island&p.area_supported,
            'strict_area_supported_without_GB_Ireland':island&p.area_supported&~p.landmass_id.isin(['15','26'])}
    counts=p[groups['strict_area_supported']].groupby('speciesKey').size()
    coverage[mode+'_area_supported_rows']=coverage.taxon_key.map(counts).fillna(0).astype(int)
    assert coverage[mode+'_area_supported_rows'].sum()==int(groups['strict_area_supported'].sum())
    for cohort,mask in groups.items():
        q=p[mask];assert len(q)>0
        for ref in ['M5_host_composition','B5_sea_terrain','T1_port_connectivity']:
            a=((q.J1_natural_transport-q.target)**2).groupby(q.species).mean().mean()
            b=((q[ref]-q.target)**2).groupby(q.species).mean().mean()
            rows.append(dict(mode=mode,cohort=cohort,reference=ref,rows=len(q),cells=q.cell.nunique(),species=q.species.nunique(),islands=q.landmass_id.nunique(),joint_mse=a,reference_mse=b,improvement_pct=100*(1-a/b)))
d=pd.DataFrame(rows);d.to_csv(O/'joint_area_supported_islands.csv',index=False)
coverage.to_csv(O/'joint_area_supported_all225.csv',index=False)
for mode in d['mode'].unique():
    old=pd.read_csv(O/f'joint_{mode}_comparison.csv').query("cohort == 'strict_islands'").set_index('reference')
    new=d[(d['mode']==mode)&(d.cohort=='strict_original')].set_index('reference')
    np.testing.assert_allclose(new.loc[old.index,'improvement_pct'],old.improvement_pct)
print(d.to_string(index=False))
(O/'joint_area_supported_validation.json').write_text(json.dumps(dict(status='PASS_AREA_SUBSET_RECONCILIATION',
    checks=['5282 unique cell flags','landmass-part areas reconstruct cell land fractions','original strict island scores reproduce previous results'],
    limits='Post-main-result fixed-threshold geographic sensitivity; subset scores without model refitting or independent occupancy verification; does not remove observation bias.'),indent=2))

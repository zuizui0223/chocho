"""Audit port-feature coverage on actual temporal evaluation rows, not the cell pool."""
from pathlib import Path
import json
import pandas as pd

R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
p=pd.read_parquet(R/'work/temporal/forecast/predictions_temporal.parquet')
f=pd.read_parquet(W/'transport_forecast_features.parquet')
f=f[f.period.eq(1)].copy()
q=p.merge(f,on='cell',validate='many_to_one')
assert len(q)==len(p) and q.indexed_port_support_unknown.notna().all()
s=pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']]
q=q.merge(s,on='cell',how='left',validate='many_to_one')
island=q.strict_single.fillna(False)&q.landmass_area_km2.lt(7e6)
cohorts={'all':pd.Series(True,index=q.index),'strict_islands':island,
         'strict_islands_without_GB_Ireland':island&~q.landmass_id.isin(['15','26'])}
rows=[]
for name,mask in cohorts.items():
    a=q[mask];b=a[a.indexed_port_support_unknown.eq(0)]
    rows.append(dict(cohort=name,rows=len(a),supported_rows=len(b),
        cells=a.cell.nunique(),supported_cells=b.cell.nunique(),
        landmasses=a.landmass_id.nunique(),supported_landmasses=b.landmass_id.nunique(),
        species=a.species.nunique(),supported_species=b.species.nunique()))
summary=pd.DataFrame(rows);summary.to_csv(O/'transport_evaluation_support.csv',index=False)
by=q[island].groupby('landmass_id').agg(rows=('cell','size'),cells=('cell','nunique'),
    species=('species','nunique'),unknown_rows=('indexed_port_support_unknown','sum'))
by['supported_rows']=by.rows-by.unknown_rows
assert by.rows.sum()==int(island.sum())
assert by.supported_rows.sum()==int((island&q.indexed_port_support_unknown.eq(0)).sum())
by.to_csv(O/'transport_evaluation_island_support.csv')
(O/'transport_evaluation_support_audit.json').write_text(json.dumps(dict(
    status='PASS_EVALUATION_SUPPORT_AUDIT',period=1,index_years=[2015,2018],
    checks=['many-to-one feature join preserves every temporal test row','no missing support flags',
            'island totals agree with row-level masks'],
    limits='Coverage only. Missing index is not absence of ports or shipping. Not independent island replication.'),indent=2))
print(summary.to_string(index=False));print(by.to_string())

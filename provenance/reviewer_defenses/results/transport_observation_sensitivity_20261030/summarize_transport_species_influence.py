"""Descriptive paired species influence; no independent-species uncertainty claim."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[2]
out=root/'outputs'
mode=sys.argv[1] if len(sys.argv)>1 else 'temporal'
rows=[]
for path in sorted(out.glob(f'transport_{mode}_*_species_errors.csv')):
    cohort=path.name[len(f'transport_{mode}_'):-len('_species_errors.csv')]
    d=pd.read_csv(path)
    assert d.species.is_unique
    for baseline in ['M5_host_composition','T0_port_location']:
        a=d[baseline+'_sq'].to_numpy()
        b=d.T1_port_connectivity_sq.to_numpy()
        assert np.isfinite(a).all() and np.isfinite(b).all()
        n=len(a)
        if n<2: continue
        loo=100*(1-(b.sum()-b)/(a.sum()-a))
        actual=100*(1-b.mean()/a.mean())
        # Independent explicit deletion check of every species.
        direct=np.array([100*(1-np.delete(b,i).mean()/np.delete(a,i).mean()) for i in range(n)])
        np.testing.assert_allclose(loo,direct,atol=1e-10)
        rows.append(dict(cohort=cohort,baseline=baseline,species=n,
                         improved_species=int((b<a).sum()),worsened_species=int((b>a).sum()),
                         improvement_pct=actual,leave_one_species_out_min_pct=loo.min(),
                         leave_one_species_out_max_pct=loo.max(),
                         omitted_at_min=d.species.iloc[loo.argmin()],omitted_at_max=d.species.iloc[loo.argmax()]))
result=pd.DataFrame(rows)
assert len(result)>0
result.to_csv(out/f'transport_{mode}_species_influence.csv',index=False)
print(result.to_string(index=False))

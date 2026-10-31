from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
for name,digest in json.loads((R/'snapshot_manifest.json').read_text()).items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
c=pd.read_csv(R/'regional_host_background_species_comparisons.csv');s=pd.read_csv(R/'regional_host_background_comparison_summary.csv')
assert np.allclose(c.host_rho-c.background_rho,c.difference)
for r in s.itertuples():
    q=c[(c.cohort==r.cohort)&(c.feature==r.feature)&c.regions.ge(r.minimum_regions)]
    assert len(q)==r.species and np.isclose(q.difference.mean(),r.mean_difference)
    assert np.isclose(q.host_rho.mean(),r.host_mean_rho) and np.isclose(q.background_rho.mean(),r.background_mean_rho)
    assert q.partial_rho.notna().sum()==r.partial_species and np.isclose(q.partial_rho.mean(),r.mean_partial_rho)
    rng=np.random.default_rng(20261030);v=q.difference.to_numpy();ix=rng.integers(len(q),size=(2000,len(q)));lo,hi=np.quantile(v[ix].mean(axis=1),[.025,.975])
    assert np.isclose(lo,r.difference_ci_low) and np.isclose(hi,r.difference_ci_high)
raw=pd.read_csv(R/'original_regional_host_distance_correlations.csv');summary=pd.read_csv(R/'original_regional_host_correlation_summary.csv')
for r in summary.itertuples():
    q=raw[(raw.cohort==r.cohort)&(raw.feature==r.feature)&raw.regions.ge(r.minimum_regions)]
    assert len(q)==r.species and np.isclose(q.rho.mean(),r.mean_rho)
print('PASS snapshot hashes and all raw/paired summary means; paired intervals and partial-support counts')

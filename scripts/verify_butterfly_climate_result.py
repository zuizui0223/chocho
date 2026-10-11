#!/usr/bin/env python3
"""Reconstruct resource counts independently of the geographic-null engine.

Inputs are the downloaded exact climate-null artifact and frozen submission
source bundle. Requires numpy, pandas and scipy in the verification environment.
This checks baseline arithmetic and stored traces, not latent MCMC state space.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--result-root', type=Path, required=True)
    parser.add_argument('--frozen-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(next(args.result_root.rglob('climate_geonull_result.json')).read_text())
    hosts = pd.read_csv(next(args.result_root.rglob('insect_host_accepted.csv')), dtype=str)
    hosts = hosts.drop_duplicates(['insect_species', 'accepted_plant_name_id'])
    columns = ['accepted_plant_name_id', 'area_code_l3']
    native = pd.read_csv(next(args.result_root.rglob('native_extant_nondoubtful_wgsrpd3.csv')), dtype=str)[columns].drop_duplicates()
    contemporary = pd.read_csv(next(args.result_root.rglob('extant_nondoubtful_wgsrpd3.csv')), dtype=str)[columns].drop_duplicates()
    frozen = pd.read_csv(args.frozen_root / 'data/frozen/figure_sources/anthropogenic_species_metrics.csv').set_index('species')
    species = sorted(frozen.index)
    focal = hosts[hosts.insect_species.isin(species)].dropna(subset=['family'])
    plants = sorted(focal.accepted_plant_name_id.unique())
    regions = sorted(contemporary.area_code_l3.unique())
    si = {s:i for i,s in enumerate(species)}
    hi = {h:i for i,h in enumerate(plants)}
    ri = {r:i for i,r in enumerate(regions)}
    incidence = csr_matrix((np.ones(len(focal), dtype=np.int32),
        ([si[s] for s in focal.insect_species], [hi[h] for h in focal.accepted_plant_name_id])),
        shape=(len(species),len(plants)))

    def spatial_matrix(frame):
        frame = frame[frame.accepted_plant_name_id.isin(plants)]
        return csr_matrix((np.ones(len(frame), dtype=np.int32),
            ([hi[h] for h in frame.accepted_plant_name_id], [ri[r] for r in frame.area_code_l3])),
            shape=(len(plants),len(regions)))

    N = spatial_matrix(native)
    I = spatial_matrix(contemporary) - N
    assert (I.data >= 0).all(), 'Native geography must be a subset of contemporary geography'
    Nc, Ic = (incidence @ N).toarray(), (incidence @ I).toarray()
    B = (Ic * (Nc > 0)).sum(1)
    U = ((Ic > 0) & (Nc == 0)).sum(1)
    D = (np.maximum(Ic-1, 0) * (Nc == 0)).sum(1)
    singles = ((Ic == 1) & (Nc == 0)).sum(1)
    T = Ic.sum(1)
    assert np.array_equal(T, B + D + U)
    adequate = [s for s in species if frozen.loc[s,'introduced_added_units'] > 0 and frozen.loc[s,'host_taxonomy_lower_bound_adequate']]
    rich = [s for s in adequate if frozen.loc[s,'resolved_host_species'] >= 10]
    panels = {'all_239':species, 'adequate_191':adequate, 'rich_hosts_80':rich}
    verified = {}
    for tag, members in panels.items():
        ix = [si[s] for s in members]
        values = {
            'introduced_host_region_incidences':int(T[ix].sum()),
            'introduced_incidences_already_covered_by_native_hosts':int(B[ix].sum()),
            'duplicated_introduced_incidences_outside_native_resource':int(D[ix].sum()),
            'added_units':int(U[ix].sum()),
            'single_host_added_units':int(singles[ix].sum()),
        }
        for key,value in values.items():
            assert value == report['baseline'][tag][key], (tag,key,value)
        verified[tag] = values
    checked = 0
    for name, model in report['models'].items():
        traces = json.loads(next(args.result_root.rglob(name + '_trace.json')).read_text())
        for tag, rows in traces.items():
            for key, comparison in model['comparisons'][tag].items():
                values = np.asarray([r[key] for r in rows])
                observed = comparison['observed']
                assert np.isclose(values.mean(), comparison['null_mean'], rtol=0, atol=1e-9)
                assert np.isclose((1+(values >= observed).sum())/(len(values)+1), comparison['p_upper'])
                assert np.isclose((1+(values <= observed).sum())/(len(values)+1), comparison['p_lower'])
                assert np.allclose(np.quantile(values,[.025,.975]),
                    [comparison['null_q025'],comparison['null_q975']], rtol=0, atol=1e-9)
            frame = pd.DataFrame(rows)
            assert ((frame['added_units'] + frame['introduced_incidences_already_covered_by_native_hosts'] +
                     frame['duplicated_introduced_incidences_outside_native_resource']) == frame['introduced_host_region_incidences']).all()
            checked += len(rows)
    result = {'status':'PASS_INDEPENDENT_SPARSE_MATRIX_AND_TRACE_AUDIT',
              'method':'SciPy integer sparse matrix products, not the geographic null engine',
              'panels':verified, 'panel_replicates_checked':checked}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()

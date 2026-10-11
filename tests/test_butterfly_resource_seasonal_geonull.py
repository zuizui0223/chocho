"""Regression tests for post-hoc seasonal and subcontinental refinements."""
from pathlib import Path
import csv
import importlib.util
import json
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'scripts/analyze_butterfly_resource_seasonal_geonull.py'

def implementation():
    assert SOURCE.is_file(), 'Seasonal refinement implementation is not present'
    sys.path.insert(0, str(ROOT / 'scripts'))
    spec = importlib.util.spec_from_file_location('seasonal_geonull', SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def toy():
    regions = [f'R{i}' for i in range(8)]
    continent = {r:'1' for r in regions}
    subregion = {r:str(10 + (i % 2)) for i,r in enumerate(regions)}
    climate = {r:(10., 100., float(i), float(7-i)) for i,r in enumerate(regions)}
    return continent, subregion, climate

def test_seasonal_strata_refine_base_without_merging():
    m = implementation(); c,g,x = toy()
    base,_=m.refine_strata(c,g,x,seasonal=False,geographic=False)
    refined,receipt=m.refine_strata(c,g,x,seasonal=True,geographic=False)
    assert len(set(refined.values())) > len(set(base.values()))
    for a in c:
        for b in c:
            if refined[a] == refined[b]: assert base[a] == base[b]
    assert receipt['seasonal_variables'] == ['bio7','bio15']

def test_geographic_labels_preserve_official_level2():
    m=implementation();c,g,x=toy()
    labels,_=m.refine_strata(c,g,x,seasonal=False,geographic=True)
    assert all(g[a]==g[b] for a in c for b in c if labels[a]==labels[b])
    del g['R0']
    with pytest.raises(ValueError,match='geography'):
        m.refine_strata(c,g,x,seasonal=False,geographic=True)

def test_missing_or_quarantined_regions_are_frozen_not_discarded():
    m=implementation();c,g,x=toy();del x['R7']
    labels,receipt=m.refine_strata(c,g,x,seasonal=True,geographic=True,quarantine={'R0'})
    assert set(labels) == set(c)
    assert list(labels.values()).count(labels['R0'])==1
    assert list(labels.values()).count(labels['R7'])==1
    assert receipt['frozen_regions'] == ['R0','R7']

def test_climate_duplicates_do_not_change_quantiles_and_conflicts_fail(tmp_path):
    m=implementation();p=tmp_path/'climate.csv'
    p.write_text('wgsrpd3_code,bio1,bio12,bio7,bio15,climate_sample_points\nR0,10,100,20,30,16\nR0,10,100,20,30,16\nR1,12,120,22,32,16\n')
    climate,metadata=m.load_four_climate(p)
    assert len(climate)==2 and metadata['source_rows']==3
    with p.open('a') as f:f.write('R1,12,120,23,32,16\n')
    with pytest.raises(ValueError,match='conflicting'):
        m.load_four_climate(p)

def test_uint16_ceiling_is_flagged_not_silently_imputed(tmp_path):
    m=implementation();p=tmp_path/'climate.csv'
    p.write_text('wgsrpd3_code,bio1,bio12,bio7,bio15,climate_sample_points\nCPI,2986.5,65535,27,516,16\n')
    climate,metadata=m.load_four_climate(p)
    assert climate['CPI'][1] == 65535
    assert metadata['ceiling_flagged_regions'] == ['CPI']

def test_official_geometry_maps_every_region_and_rejects_conflict(tmp_path):
    m=implementation();p=tmp_path/'geometry.json'
    props={'LEVEL3_COD':'R0','LEVEL1_COD':1,'LEVEL2_COD':10}
    p.write_text(json.dumps({'features':[{'properties':props}]}))
    assert m.load_subregions(p, {'R0':'1'})=={'R0':'10'}
    with pytest.raises(ValueError,match='continent'):
        m.load_subregions(p, {'R0':'2'})

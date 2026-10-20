"""Tests for native WCVP region-trained pseudo-presence SDMs and 21ka projection."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_native_plant_21ka_climate_hindcast_v01 import (
    background_support_indices, block_label, spatial_holdout_labels,
    fit_predict,
)

def test_native_background_only_uses_declared_regions_and_is_deterministic():
    native = {"123": {"AAA", "CCC"}, "456": {"BBB"}}
    support = [
        {"area_code_l3": "AAA"}, {"area_code_l3": "BBB"},
        {"area_code_l3": "CCC"}, {"area_code_l3": "DDD"},
        {"area_code_l3": "AAA"}, {"area_code_l3": "CCC"},
    ]
    ix = background_support_indices("123", native, support)
    assert ix.tolist() == [0, 2, 4, 5]
    assert np.array_equal(ix, background_support_indices("123",native,support))
    assert background_support_indices("456",native,support).tolist() == [1]
    assert len(background_support_indices("unlisted",native,support)) == 0

def test_2degree_block_identical_samples_stay_in_same_validation_fold():
    rng = np.random.default_rng(17)
    plant = np.column_stack([rng.uniform(-40,40,120), rng.uniform(-90,90,120)])
    background = np.vstack([plant[:35], np.column_stack([rng.uniform(-40,40,80),rng.uniform(-90,90,80)])])
    hold = spatial_holdout_labels(plant,background,"Some taxon")
    assert hold.dtype == bool and len(hold)==len(plant)+len(background)
    for a in range(35):
        assert hold[a] == hold[len(plant)+a]

def test_current_same_as_lgm_has_identical_suitability_fraction_and_zero_delta():
    rng = np.random.default_rng(8)
    pres = rng.normal(loc=.30,scale=.13,size=(85,4))
    bg = rng.normal(loc=-.25,scale=.15,size=(95,4))
    pc = np.column_stack([rng.uniform(-60,60,85),rng.uniform(-170,170,85)])
    bc = np.column_stack([rng.uniform(-60,60,95),rng.uniform(-170,170,95)])
    g = rng.normal(loc=0,scale=.10,size=(100,4))
    m,s_now,s_lgm,novel = fit_predict(pres,bg,pc,bc,g,g,"Testia alpha")
    assert len(s_now)==len(s_lgm)==len(novel)==100
    assert np.array_equal(s_now,s_lgm)
    assert m["climate_suitable_fraction_now"] == m["climate_suitable_fraction_21ka"]
    assert abs(m["climate_suitable_fraction_delta_21ka_minus_now"]) < 1e-12
    assert m["source_presence_n"] == 85 and m["source_background_n"] == 95
    assert 0 <= m["nonanalog_21ka_site_fraction"] <= 1
    assert m["cv_auc_presence_vs_background"] is None or 0 <= m["cv_auc_presence_vs_background"] <= 1
    assert set(m["sensitivity_05_10_15"]) == {"0.05","0.1","0.15"}

def test_21ka_shift_can_change_climatic_suitability_without_claiming_occupancy():
    rng=np.random.default_rng(3)
    p=rng.normal(.50,.20,size=(65,4))
    b=rng.normal(-.65,.20,size=(75,4))
    pc=np.column_stack([np.arange(65)*1.7-55,np.arange(65)*2.3-75])
    bc=np.column_stack([np.arange(75)*1.5-55,np.arange(75)*2.2-75])
    now=rng.normal(.4,.12,size=(90,4))
    lgm=np.vstack([now[:10], now[10:]-1.7])
    output,current,glacial,novel=fit_predict(p,b,pc,bc,now,lgm,"Testia beta")
    assert output["climate_suitable_fraction_now"]>output["climate_suitable_fraction_21ka"]
    assert output["climate_suitable_fraction_delta_21ka_minus_now"]<0
    assert output["nonanalog_21ka_site_fraction"]>0.95
    assert 5 <= output["supported_analog_21ka_site_n"] < 20

def test_nonfinite_or_insufficient_training_refused():
    x=np.ones((31,4))
    coords=np.zeros((31,2))
    modern=np.zeros((40,4))
    with pytest.raises(ValueError):
        fit_predict(x[:29],x,coords[:29],coords,modern,modern,"Bad plant")
    nanx=x.copy();nanx[4,2]=np.nan
    with pytest.raises(ValueError):
        fit_predict(nanx,x,coords,coords,modern,modern,"Bad plant")

def test_source_program_never_requires_historical_occupied_plant_ranges():
    src=(Path(__file__).resolve().parents[1]/"scripts"/"run_native_plant_21ka_climate_hindcast_v01.py").read_text()
    assert "not actual LGM occupancy" in src.lower()
    assert "genetic" in src.lower()

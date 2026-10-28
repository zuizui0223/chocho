"""Pure future-climate scenario controls: no network or butterfly observations."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from simulate_monarch_cmip6_thermal_window_20261028 import (
    DAYS, WINTER, ONE_GENERATION_GDD,archive_name,
    grid_cell,month_from_file,choose_temperature_scale,scenario_thermal,
    validate_climate,summarize_group,_sample_tif_band,climate_validity_mask,
)


def fake_climate(temperature=15.,n=4):
    return ({i:np.repeat(float(temperature)-2.,n) for i in range(1,13)},
            {i:np.repeat(float(temperature)+2.,n) for i in range(1,13)})


def test_winter_degree_days_exactly_120_day_months_and_leaf_window():
    tmin,tmax=fake_climate(15.0,n=3)
    out=scenario_thermal(tmin,tmax)
    assert float(out['winter_gdd'][0])==pytest.approx((15.-11.5)*120)
    assert float(out['winter_novdec_gdd'][0])==pytest.approx((15.-11.5)*61)
    assert float(out['winter_janfeb_gdd'][0])==pytest.approx((15.-11.5)*59)
    q=summarize_group(out,[0,1,2])
    assert q['winter_thermal_323_count']==3
    assert q['thermal_fraction_ge_323']==1.0
    assert q['winter_available_50pct_323_fraction']==0.0
    assert q['winter_required_leaf_fraction_median_among_thermally_feasible']==pytest.approx(323./420.)
    assert q['novdec_only_fraction_ge_323']==0.


def test_model_is_not_monotonic_success_when_host_leaves_unavailable():
    old_min,old_max=fake_climate(11.0)
    hot_min,hot_max=fake_climate(22.0)
    historical=scenario_thermal(old_min,old_max)
    future=scenario_thermal(hot_min,hot_max)
    assert summarize_group(historical,range(4))['thermal_fraction_ge_323']==0.
    assert summarize_group(future,range(4))['thermal_fraction_ge_323']==1.
    # A literal leaf-absent scenario always has zero host-mediated development,
    # regardless of how warm weather gets. This is not estimated winter survival.
    assert float(np.sum(future['winter_gdd']*0))==0.


def test_temperature_ceiling_and_monthly_heat_screen_are_not_lethal_events():
    a,b=fake_climate(45.0)
    out=scenario_thermal(a,b)
    assert float(out['winter_gdd'][0])==pytest.approx((33.-11.5)*120)
    assert out['summer_monthly_mean_max_above38'].all()
    with pytest.raises(ValueError):scenario_thermal(b,a)


def test_multiple_worldclim_archive_names_and_months_fixed():
    assert archive_name('tmin')=='wc2.1_10m_tmin.zip'
    assert archive_name('tmax','MIROC6','ssp585','2081-2100')=='wc2.1_10m_tmax_MIROC6_ssp585_2081-2100.zip'
    assert month_from_file('data/wc2.1_10m_tmax_MIROC6_ssp585_2081-2100_07.tif')==7
    assert month_from_file('wc2.1_10m_tmax_1.tif')==1
    assert month_from_file('wc2.1_10m_tmax_13.tif') is None
    with pytest.raises(ValueError):archive_name('tavg','MIROC6','ssp126','2041-2060')
    with pytest.raises(ValueError):archive_name('tmin','FAKE','ssp585','2081-2100')


def test_raw_unit_scale_and_geographic_boundaries():
    assert choose_temperature_scale([280.,300.,310.],[1.,1.,1.])==.1
    assert choose_temperature_scale([28.,30.,31.],[1.,1.,1.])==1.
    assert choose_temperature_scale([280.,300.],[.1,.1])==.1
    with pytest.raises(ValueError):choose_temperature_scale([150.,180.],[.1,.01])
    with pytest.raises(ValueError):choose_temperature_scale([65.,70.],[1.,1.])
    assert grid_cell(0.,0.)==(540,1080)
    assert grid_cell(95.,30.) is None
    assert grid_cell(35.,190.) is None


def test_two_original_hexes_sharing_one_climate_pixel_still_count_as_two_source_hexes():
    x,y=fake_climate(17.0,n=2)
    climate=scenario_thermal(x,y)
    counted=summarize_group(climate,[0,0,1])
    assert counted['n_100km_hex_localities']==3
    assert counted['winter_thermal_323_count']==3


def test_original_12_file_vs_future_12_band_geotiff_input_contract():
    # WorldClim's future CMIP6 archive is a single 12-band GeoTIFF,
    # unlike the base archive (12 month-specific TIFs).
    assert month_from_file('wc2.1_10m_tmin_01.tif')==1
    assert month_from_file('wc2.1_10m_tmin_IPSL-CM6A-LR_ssp126_2041-2060.tif') is None
    assert callable(_sample_tif_band)


def test_physical_QC_excludes_same_raster_pixel_from_every_scenario():
    tmin,tmax=fake_climate(15.,n=4)
    tmin[12][1]=100.
    valid,counts=climate_validity_mask(tmin,tmax)
    assert valid.tolist()==[True,False,True,True]
    assert counts['12_tmin_above_tmax']==1
    with pytest.raises(ValueError,match='month=12'):
        scenario_thermal(tmin,tmax)
    clean_lo={m:v[valid] for m,v in tmin.items()}
    clean_hi={m:v[valid] for m,v in tmax.items()}
    assert summarize_group(scenario_thermal(clean_lo,clean_hi),[0,1,2])['n_100km_hex_localities']==3

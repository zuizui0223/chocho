"""Offline controls for native-host-relative country climate availability, no API."""
import sys
import json
from pathlib import Path

import numpy as np
import pytest
from rasterio.transform import from_origin
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_butterfly_73_country_climate_area_null_20261028 import (
    area_fraction_outside, paired_compare, load_bounds, BOUNDS_SCHEMA,
    EXPECTED_ORIGINAL_JSON_SHA256, EXPECTED_WORLDCLIM_ZIP_SHA256,
)


def fake_climate():
    return {
        "bio1":np.array([[5.,5.,15.,15.]]*4),
        "bio12":np.full((4,4),500.),
        "valid":np.ones((4,4),dtype=bool),
        "transform":from_origin(0,4,1,1),
    }


def rectangles():
    return {"Test butterfly":{"bounds":{
        "area_weighted_5_95":{"bio1":[10.,20.],"bio12":[400.,600.]},
        "area_weighted_10_90":{"bio1":[12.,18.],"bio12":[450.,550.]},
        "equal_pixel_5_95":{"bio1":[10.,20.],"bio12":[400.,600.]},
    }}}


def cases():
    return [
        {"species":"Test butterfly","country":"AA","source_class":"novel",
         "n":10,"outside_area95":8,"outside_area90":9,"outside_equal95":8,
         "strict_n":3,"strict_outside_area95":2,"strict_outside_area90":3,"strict_outside_equal95":2},
        {"species":"Test butterfly","country":"BB","source_class":"analogue",
         "n":10,"outside_area95":2,"outside_area90":3,"outside_equal95":2,
         "strict_n":3,"strict_outside_area95":1,"strict_outside_area90":1,"strict_outside_equal95":1},
    ]


def control():
    return area_fraction_outside(
        {"AA":box(0,0,2,4),"BB":box(2,0,4,4)},
        fake_climate(),cases(),rectangles(),
    )


def test_country_land_is_spatial_opportunity_not_observation_count():
    x=control()
    assert x["Test butterfly","AA"]["country_climate_pixels"]==8
    assert x["Test butterfly","BB"]["country_climate_pixels"]==8
    assert x["Test butterfly","AA"]["status"]=="INSUFFICIENT_10MIN_LAND_PIXEL_CENTERS"
    # Artificial tiny raster below real-country minimum is censored, not zero.
    assert "expected_outside_area95" not in x["Test butterfly","AA"]


def test_geographic_counterfactual_offsets_raw_case_difference():
    # Explicit synthetic 20-pixel per country panels, separated by longitude.
    c={"bio1":np.array([[5.]*5+[15.]*5]*5),"bio12":np.full((5,10),500.),
       "valid":np.ones((5,10),dtype=bool),"transform":from_origin(0,5,1,1)}
    x=area_fraction_outside({"AA":box(0,0,5,5),"BB":box(5,0,10,5)},
                            c,cases(),rectangles())
    assert x["Test butterfly","AA"]["expected_outside_area95"]==1.
    assert x["Test butterfly","BB"]["expected_outside_area95"]==0.
    out=paired_compare(cases(),x)
    assert out["paired_species"]==1
    assert out["mean_observed_country_class_difference"]==pytest.approx(.6)
    assert out["mean_geographic_country_class_difference"]==pytest.approx(1.)
    assert out["mean_observed_minus_country_area_reference"]==pytest.approx(-.4)
    assert out["geography_adjusted_species_directions"]["negative"]==1
    assert out["p_values"] is None
    assert paired_compare(cases(),x,minimum_records=11)["paired_species"]==0
    assert paired_compare(cases(),x,strict=True)["paired_species"]==1


def test_grid_country_nodata_is_not_calculated_as_climate_nonanalog():
    c={"bio1":np.array([[5.]*5+[15.]*5]*5),"bio12":np.full((5,10),500.),
       "valid":np.ones((5,10),dtype=bool),"transform":from_origin(0,5,1,1)}
    c["valid"][:3,:5]=False
    x=area_fraction_outside({"AA":box(0,0,5,5),"BB":box(5,0,10,5)},
                            c,cases(),rectangles())
    assert x["Test butterfly","AA"]["country_climate_pixels"]==10
    assert x["Test butterfly","AA"]["expected_outside_area95"]==1.
    c["valid"][:,:5]=False
    x=area_fraction_outside({"AA":box(0,0,5,5),"BB":box(5,0,10,5)},
                            c,cases(),rectangles())
    assert x["Test butterfly","AA"]["status"]=="INSUFFICIENT_10MIN_LAND_PIXEL_CENTERS"


def test_frozen_bounds_must_match_source_and_all_20_species(tmp_path):
    path=tmp_path/"bounds.json"
    original={"schema":BOUNDS_SCHEMA,
              "origin":{"original_source_json_sha256":EXPECTED_ORIGINAL_JSON_SHA256,
                        "raster_zip_sha256":EXPECTED_WORLDCLIM_ZIP_SHA256},
              "references":rectangles()}
    path.write_text(json.dumps(original))
    with pytest.raises(ValueError,match="20 original species"):
        load_bounds(path,cases())
    original["origin"]["original_source_json_sha256"]="untrusted"
    path.write_text(json.dumps(original))
    with pytest.raises(ValueError,match="source identity"):
        load_bounds(path,cases())


def test_source_verified_real_country_climate_area_null_receipt_is_not_host_effect():
    root = Path(__file__).resolve().parents[1]
    receipt = root / 'provenance/reviewer_defenses/results/butterfly_73_country_land_climate_null_20261028.json'
    r = json.loads(receipt.read_text(encoding='utf8'))
    assert r['schema'] == 'chocho_73_country_land_climate_area_null_20261028_v01'
    assert (r['original_observation_cases'],r['original_GBIF_point_indices']) == (73,6076)
    assert r['country_geography']['source_countries'] == 51
    assert r['main']['cases_used_before_pairing'] == 73
    assert r['main']['paired_species'] == 10
    assert r['main']['p_values'] is None
    assert r['source']['natural_earth_50m_exact_sha256'] == 'b8d421aca6e9e08e8cdf09cc26af111cc3e0deba4fe915611d58ade71e8a4db0'
    assert r['main']['mean_observed_country_class_difference'] == pytest.approx(.42415155141720556)
    assert r['main']['mean_geographic_country_class_difference'] == pytest.approx(.41068989706675174)
    assert r['main']['mean_observed_minus_country_area_reference'] == pytest.approx(.013461654350453805)
    assert r['main']['geography_adjusted_species_directions'] == {'positive':5,'negative':5,'zero':0}
    assert all('lat' not in row and 'lon' not in row and 'observer' not in row
               for row in r['case_summary'])
    assert 'NOT GBIF survey effort' in r['main']['inference_gate'] or \
           'Not GBIF observation-effort matched' in r['main']['inference_gate']

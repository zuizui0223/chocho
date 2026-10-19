"""Synthetic tests for proximity+tethered geographic range rearrangements.

No synthetic null statistic from these tests is a biological result.
"""
from __future__ import annotations

import importlib.util
import json
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SOURCE = ROOT / "scripts/analyze_butterfly_spatially_local_plant_range_null.py"
spec = importlib.util.spec_from_file_location("spatial_host_null", SOURCE)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_antimeridian_is_near_not_world_apart() -> None:
    assert mod.haversine_km((1., 179.5), (1., -179.5)) < 120
    assert mod.haversine_km((0., 0.), (0., 180.)) > 19000
    assert mod.haversine_km((10., 20.), (10., 20.)) == pytest.approx(0)


def test_geometric_codes_and_hierarchy_are_exact(tmp_path: Path) -> None:
    path = tmp_path / "regions.geojson"
    def feature(code: str, area: str, lon: float, lat: float) -> dict:
        return {
            "type":"Feature",
            "properties":{"LEVEL3_COD":code, "LEVEL1_COD":area},
            "geometry":{"type":"Polygon","coordinates":[[
                [lon,lat],[lon+0.1,lat],[lon+0.1,lat+0.1],
                [lon,lat+0.1],[lon,lat],
            ]]},
        }
    data = {"type":"FeatureCollection","features":[
        feature("AAA","1",15.,30.),
        feature("BBB","1",15.5,30.),
    ]}
    path.write_text(json.dumps(data), encoding="utf-8")
    coords, meta = mod.load_geography(path,["AAA","BBB"],{0:"1",1:"1"})
    assert len(coords)==2 and meta["mapped_focal_regions"]==2
    assert len(meta["source_geometry_sha256"])==64
    with pytest.raises(ValueError, match="Level1 hierarchy"):
        mod.load_geography(path,["AAA","BBB"],{0:"1",1:"4"})
    with pytest.raises(ValueError, match="missing Level3"):
        mod.load_geography(path,["AAA","CCC"],{0:"1",1:"1"})
    data["features"].append(feature("AAA","1",15.,30.))
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        mod.load_geography(path,["AAA","BBB"],{0:"1",1:"1"})


def test_local_swap_preserves_plant_degree_region_degree_and_original_tether() -> None:
    # Three neighboring regions (approx 110km each) and two plants.
    coords = {0:(0., 0.), 1:(0., 1.), 2:(0., 2.), 3:(0., 8.)}
    near = mod.neighbour_index(coords, 250)
    assert near[0] == frozenset((0,1,2))
    assert 3 not in near[0]

    native = {"a":{3}, "b":set()}
    original = {"a":{0}, "b":{1}}
    present = {h:set(u) for h,u in original.items()}
    tether = mod.origin_tether(original,near)
    assert 1 in tether["a"]
    assert 3 not in tether["a"]
    assert mod.admissible(("a",0),("b",1),present,native,near,tether)
    assert not mod.admissible(("a",0),("b",3),present,native,near,tether)
    group = [[("a",0),("b",1)]]
    n = mod.local_swap(group,[2],present,native,near,tether,100,random.Random(2026))
    assert n > 0
    assert sorted(len(u) for u in present.values()) == [1,1]
    assert sum(1 in u for u in present.values()) == 1
    assert sum(0 in u for u in present.values()) == 1
    assert all(not present[h]&native[h] for h in present)
    assert all(present[h] <= tether[h] for h in present)
    assert mod.count_mobile_edges(group,present,native,near,tether) == 2


def test_tether_prevents_diffusion_across_multiple_near_swaps() -> None:
    coords = {i:(0.,float(i)) for i in range(5)}
    near = mod.neighbour_index(coords, 130)
    original = {"a":{0}, "b":{1}, "c":{2}, "d":{3}, "e":{4}}
    t = mod.origin_tether(original,near)
    assert 2 not in t["a"]   # a cannot drift 0->1->2 via multiple local swaps
    present = {k:set(v) for k,v in original.items()}
    assert not mod.admissible(("a",1),("c",2),present,{},near,t)
    for h,vals in original.items():
        assert t[h] >= vals


def test_no_geographic_mobile_edges_fails_before_null_sampling() -> None:
    coords = {0:(0.,0.),1:(0.,40.)}
    near = mod.neighbour_index(coords,750)
    source = {"a":{0}, "b":{1}}
    tether = mod.origin_tether(source,near)
    groups = [[("a",0),("b",1)]]
    assert mod.count_mobile_edges(groups,source,{},near,tether) == 0


def test_lag_one_autocorrelation_is_diagnostic_only() -> None:
    assert mod.lag1([1.,1.,1.]) is None
    assert mod.lag1([1.,2.]) is None
    assert -1. <= mod.lag1([1.,2.,3.,2.,1.,2.]) <= 1.


def test_all_caliper_choices_are_named_in_frozen_protocol() -> None:
    assert mod.CALIPERS_KM == (1500,750,3000)
    p = (
        ROOT /
        "provenance/reviewer_defenses/protocols/SPATIALLY_LOCAL_RESOURCE_ABSORPTION_20261015.md"
    )
    txt = p.read_text(encoding="utf-8")
    assert all(f"{r:,} km" in txt for r in mod.CALIPERS_KM)
    assert "NOT_EVALUABLE" in txt

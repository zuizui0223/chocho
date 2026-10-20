"""Geographical containment, slivers and ISO aggregation remain explicit."""
from pathlib import Path
import sys

import pytest
from shapely.geometry import box, mapping

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_country_wgsrpd_crosswalk_20261020 import (
    country_code, dissolve_countries, intersect_crosswalk, equal_area_projector,
)


def feature(code, geometry, *, raw=None, name="unit"):
    return {"type": "Feature", "properties": {
        "ISO_A2": code if raw is None else raw, "ISO_A2_EH": code,
        "NAME": name, "ADMIN": name,
    }, "geometry": mapping(geometry)}


def test_iso_extended_field_retains_territory_and_rejects_real_conflict():
    assert country_code({"ISO_A2": "FR-973", "ISO_A2_EH": "GF"}) == "GF"
    assert country_code({"ISO_A2": "-99", "ISO_A2_EH": "GB"}) == "GB"
    assert country_code({"ISO_A2": "-99", "ISO_A2_EH": "-99"}) is None
    with pytest.raises(ValueError, match="Conflicting"):
        country_code({"ISO_A2": "FR", "ISO_A2_EH": "GF"})


def test_near_antimeridian_source_vertex_does_not_wrap_across_globe():
    projected = equal_area_projector()(box(-180.000026, 60, -179, 61))
    assert projected.bounds[2] < 0
    assert 5e9 < projected.area < 7e9


def test_duplicate_iso_units_dissolve_without_double_area():
    countries, audit = dissolve_countries([
        feature("GB", box(0, 0, 2, 1), name="England"),
        feature("GB", box(1, 0, 3, 1), name="Wales"),
    ], project=lambda geom: geom)
    assert not audit["excluded_map_units"]
    assert countries["GB"]["geometry"].area == 3
    assert countries["GB"]["names"] == ["England", "Wales"]


def test_small_country_inside_large_tdwg_is_not_lost_to_area_threshold():
    countries = {"XX": {"geometry": box(0, 0, 1, 1), "names": ["tiny"]}}
    regions = {"BIG": {"geometry": box(-10, -10, 10, 10), "level1": "1"}}
    rows, summary = intersect_crosswalk(countries, regions)
    assert rows[0]["smaller_polygon_fraction"] == 1
    assert rows[0]["material_link"] == 1
    assert rows[0]["region_fraction"] == .0025
    assert summary[0]["material_coverage"] == 1


def test_union_coverage_does_not_double_count_overlapping_tdwg_regions():
    countries = {"XX": {"geometry": box(0, 0, 10, 10), "names": ["example"]}}
    regions = {
        "AAA": {"geometry": box(0, 0, 9, 10), "level1": "1"},
        # Ten percent of the smaller (country) polygon overlaps.
        "BBB": {"geometry": box(9, 0, 109, 10), "level1": "2"},
        # This overlap duplicates AAA, so coverage must use a union.
        "CCC": {"geometry": box(0, 0, 1, 10), "level1": "1"},
    }
    rows, summary = intersect_crosswalk(countries, regions)
    by_code = {r["wgsrpd3"]: r for r in rows}
    # Country is the smaller polygon, so BBB covers 10%, not a sliver.
    assert by_code["BBB"]["material_link"] == 1
    assert summary[0]["material_coverage"] == 1
    assert summary[0]["material_intersection_sum_fraction"] == 1.1
    assert summary[0]["material_level1_codes"] == ["1", "2"]


def test_guard_only_border_overlap_does_not_become_material_region():
    countries = {"XX": {"geometry": box(0, 0, 10, 10), "names": []}}
    regions = {
        "AAA": {"geometry": box(0, 0, 9.9, 10), "level1": "1"},
        "BBB": {"geometry": box(9.9, 0, 19.9, 10), "level1": "1"},
    }
    rows, summary = intersect_crosswalk(countries, regions)
    by_code = {r["wgsrpd3"]: r for r in rows}
    assert by_code["BBB"]["material_link"] == 0
    assert by_code["BBB"]["guard_link"] == 1
    assert summary[0]["material_regions"] == ["AAA"]
    assert summary[0]["guard_regions"] == ["AAA", "BBB"]
    assert summary[0]["material_coverage"] == pytest.approx(.99)
    assert summary[0]["guard_coverage"] == pytest.approx(1)

"""Exact geographic union areas without serializing source polygons."""
import importlib.util
from pathlib import Path
import sys

import pytest
from shapely.geometry import box


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/build_country_area_atoms_20261020.py"
sys.path.insert(0, str(ROOT / "scripts"))


def module():
    assert SCRIPT.is_file(), "Country atom area builder is not implemented"
    spec = importlib.util.spec_from_file_location("country_area_atoms", SCRIPT)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_overlapping_regions_are_partitioned_and_never_double_counted():
    atoms = module().atom_areas(box(0, 0, 10, 10), {
        "AAA": box(0, 0, 6, 10), "BBB": box(4, 0, 10, 10),
    })
    assert atoms == {("AAA",): 40, ("AAA", "BBB"): 20, ("BBB",): 40}
    assert sum(atoms.values()) == 100
    assert sum(a for members, a in atoms.items() if "AAA" in members) == 60


def test_country_clipping_and_disconnected_parts_keep_exact_membership():
    atoms = module().atom_areas(box(0, 0, 10, 10), {
        "AAA": box(-5, -5, 5, 15), "BBB": box(0, 0, 5, 10),
        "CCC": box(9, 0, 11, 10), "OUT": box(20, 20, 30, 30),
    })
    assert atoms == {("AAA", "BBB"): 50, ("CCC",): 10}


def test_three_way_overlap_preserves_any_selected_mask_union():
    atoms = module().atom_areas(box(0, 0, 10, 10), {
        "AAA": box(0, 0, 6, 10), "BBB": box(4, 0, 10, 10),
        "CCC": box(5, 0, 7, 10),
    })
    assert atoms == {
        ("AAA",): 40, ("AAA", "BBB"): 10,
        ("AAA", "BBB", "CCC"): 10, ("BBB", "CCC"): 10,
        ("BBB",): 30,
    }
    selected = {"AAA", "CCC"}
    assert sum(area for members, area in atoms.items() if selected.intersection(members)) == 70


def test_empty_or_boundary_only_intersections_do_not_make_area():
    api = module()
    country = box(0, 0, 1, 1)
    assert api.atom_areas(country, {}) == {}
    assert api.atom_areas(country, {"TOUCH": box(1, 0, 2, 1)}) == {}


def test_input_order_does_not_change_atom_memberships_or_areas():
    api = module()
    regions = {"BBB": box(.4, 0, 1, 1), "AAA": box(0, 0, .6, 1)}
    a = api.atom_areas(box(0, 0, 1, 1), regions)
    b = api.atom_areas(box(0, 0, 1, 1), dict(reversed(list(regions.items()))))
    assert a == b
    assert sum(a.values()) == pytest.approx(1)

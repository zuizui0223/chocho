"""Independent geometric cases for island isolation and port-location exposure."""
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("geopandas")
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import measure_island_port_geography_20261029 as m


def test_dateline_points_use_short_arc():
    # A planar longitude subtraction would turn a 22 km gap into ~40,000 km.
    a = m.unit_xyz(np.array([[179.9, 0.0]]))
    b = m.unit_xyz(np.array([[-179.9, 0.0]]))
    d = m.chord_to_km(np.linalg.norm(a - b, axis=1))[0]
    assert 22.23 < d < 22.25


def test_long_source_segments_cover_coast_interior_and_bound_gap():
    # True closest coasts are the parallel equatorial ends, one degree apart;
    # retaining vertices only would miss an interior-to-interior nearest pair.
    mainland, _ = m.shoreline_xyz(box(-2, -1, 0, 1), 1.0)
    island, _ = m.shoreline_xyz(box(1, -0.2, 1.2, 0.2), 1.0)
    lower, upper = m.coast_gap_interval(island, cKDTree(mainland), 1.0)
    assert 110.0 < lower <= upper < 112.0
    assert upper - lower <= 1.0000001
    assert len(mainland) > 800


def test_dateline_cut_is_not_a_coast():
    # The middle of the vertical 180-degree edge is a map cut, not sea coast.
    geometry = Polygon([(179, -1), (180, -1), (180, 1), (179, 1)])
    points, qc = m.shoreline_xyz(geometry, 1.0)
    cut_midpoint = m.unit_xyz(np.array([[180.0, 0.0]]))
    dist = m.chord_to_km(cKDTree(points).query(cut_midpoint)[0])[0]
    assert 110 < dist < 112
    assert qc["removed_dateline_cut_segments"] == 1


def test_port_threshold_counts_keep_sampling_ambiguity():
    # At the equator these point distances are ~49.9, 50.1 and 55.6 km.
    coast = m.unit_xyz(np.array([[0.0, 0.0]]))
    ports = m.unit_xyz(np.array([[0.449, 0.0], [0.451, 0.0], [0.5, 0.0]]))
    result = m.port_exposure(coast, ports, 1.0)
    assert result["wpi_entries_within_50km_certain"] == 1
    assert result["wpi_entries_within_50km_possible"] == 2
    assert result["wpi_entries_within_100km_certain"] == 3


def test_nonfinite_coordinates_fail_before_tree_build():
    with pytest.raises(ValueError, match="coordinates"):
        m.unit_xyz(np.array([[np.nan, 4.0]]))


def test_southward_bowing_arc_keeps_antarctic_cap_bound_conservative():
    # The southern great-circle edge bows to ~50.4313 S between the vertices.
    coast, _ = m.shoreline_xyz(box(-10, -50, 10, -49), 1.0)
    bound = m.antarctic_cap_distance_lower(coast, 1.0)
    assert 1063 < bound < 1065
    assert bound < 1080  # A mainland gap of 1080 km must remain ambiguous.


def test_valid_coordinate_mutation_cannot_pass_locked_port_pages():
    import hashlib
    original_page = b'{"latitude":1,"longitude":2}'
    mutated_page = b'{"latitude":3,"longitude":2}'
    original_sha = hashlib.sha256(original_page).hexdigest()
    locked = "8ee3fbfde1791a0184c9000dfd5e88823ac99b269d5d2b5fc6afb546dd995cfa"
    assert m.verify_port_pages([original_sha], locked) == locked
    with pytest.raises(ValueError, match="coordinates.*changed"):
        m.verify_port_pages([hashlib.sha256(mutated_page).hexdigest()], locked)


def test_zerozero_placeholder_does_not_become_an_ocean_port():
    rows = [{"attributes": {"longitude": 0, "latitude": 0}},
            {"attributes": {"longitude": 6, "latitude": 0}}]
    points, excluded = m.usable_port_points(rows)
    assert len(points) == 1
    assert excluded == 1
    assert np.allclose(points, m.unit_xyz(np.array([[6, 0]])))

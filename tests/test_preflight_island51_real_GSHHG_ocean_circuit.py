"""Ocean topology must not invent marine routes or separate collocated anchors."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
for package in ("scipy", "geopandas"):
    pytest.importorskip(package)


SCRIPT = (Path(__file__).resolve().parents[1] / "scripts" /
          "preflight_island51_real_GSHHG_ocean_circuit.py")
spec = importlib.util.spec_from_file_location("ocean_preflight", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.mark.parametrize("side_water", [(False, False), (False, True), (True, False)])
def test_diagonal_cannot_cross_either_adjacent_land_corner(side_water):
    # A diagonal through a blocked cell corner would invent a navigable edge.
    sea = np.zeros((170, 360), dtype=bool)
    sea[80, 100] = sea[81, 101] = True
    sea[80, 101], sea[81, 100] = side_water

    ids, graph, _, _ = m.connect_ocean(sea, neighbors=8)

    assert graph[ids[80, 100], ids[81, 101]] == 0


def test_diagonal_corner_rule_also_applies_across_longitude_seam():
    sea = np.zeros((170, 360), dtype=bool)
    sea[80, 359] = sea[81, 0] = sea[80, 0] = True

    ids, graph, _, _ = m.connect_ocean(sea, neighbors=8)

    assert graph[ids[80, 359], ids[81, 0]] == 0
    assert graph[ids[80, 359], ids[80, 0]] == 1


def test_open_water_keeps_diagonal_and_cardinal_longitude_wrap():
    sea = np.zeros((170, 360), dtype=bool)
    sea[80:82, [0, 359]] = True

    ids, graph, components, _ = m.connect_ocean(sea, neighbors=8)

    assert graph[ids[80, 359], ids[81, 0]] == 1
    assert graph[ids[80, 0], ids[81, 359]] == 1
    assert graph[ids[80, 359], ids[80, 0]] == 1
    assert len(set(components)) == 1


def source_points(*, duplicate):
    # Exactly 51 declared sources; positions are hand-defined grid centers.
    points = {
        str(i + 1): {"longitude": -179.5 + i, "latitude": 0.5}
        for i in range(51)
    }
    if duplicate:
        points["2"] = {"longitude": -179.5, "latitude": 0.5}
    return {
        "status": "ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES",
        "verified_source_51_mean_port_points": points,
    }


def test_duplicate_nearest_anchors_are_reported_without_greedy_displacement(monkeypatch):
    # Only external coastline IO is replaced; graph, snapping and gate are real.
    sea = np.zeros((170, 360), dtype=bool)
    sea[84, :] = True
    monkeypatch.setattr(m, "land_mask_from_GSHHG", lambda _: (~sea, {}))

    result = m.run(Path("unused-source.zip"), source_points(duplicate=True), neighbors=8)
    first = result["island_points"]["1"]
    second = result["island_points"]["2"]
    audit = result["WPI_centroid_marine_snapping"]

    assert first["ocean_grid_node"] == second["ocean_grid_node"]
    assert first["ocean_grid_longitude"] == second["ocean_grid_longitude"] == -179.5
    assert len(result["island_points"]) == 51
    assert audit["islands_snapped"] == 51
    assert audit["unique_marine_nodes"] == 50
    assert audit["duplicate_marine_node_groups"] == {"0": ["1", "2"]}
    assert not audit["all_51_distinct_marine_nodes"]
    assert result["status"] == "INCOMPLETE_1DEG_OCEAN_CONNECTIVITY_DO_NOT_FIT_ALL_51"


def test_distinct_independent_nearest_anchors_pass_all_51_gate(monkeypatch):
    sea = np.zeros((170, 360), dtype=bool)
    sea[84, :] = True
    monkeypatch.setattr(m, "land_mask_from_GSHHG", lambda _: (~sea, {}))

    result = m.run(Path("unused-source.zip"), source_points(duplicate=False), neighbors=4)

    assert result["status"] == "VALID_1DEG_OCEAN_DISCRETE_CIRCUIT_SOURCE_TOPOLOGY"
    audit = result["WPI_centroid_marine_snapping"]
    assert audit["all_51_distinct_marine_nodes"]
    assert audit["unique_marine_nodes"] == 51
    assert audit["duplicate_marine_node_groups"] == {}

"""Aggregation preserves resolved four-neighbor sea connectivity."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
pytest.importorskip("scipy")
from shapely.geometry import box


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_island51_coastal_ocean_graph.py"
spec = importlib.util.spec_from_file_location("coastal_graph", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_two_disconnected_sea_pieces_in_one_block_keep_separate_nodes():
    # Collapsing every wet cell of a block into one node would cross land.
    sea = np.array([[True, False], [False, True]])

    arrays, audit = m.build_coastal_graph(sea, block_size_cells=2, wrap_longitude=False)

    grid = arrays["fine_cell_to_node"]
    assert grid[0, 0] != grid[1, 1]
    assert np.all(grid[~sea] == -1)
    assert arrays["node_block_row"].tolist() == [0, 0]
    assert arrays["node_block_col"].tolist() == [0, 0]
    assert audit["graph_nodes"] == 2
    assert audit["number_of_components"] == 2
    assert len(arrays["edge_i"]) == 0


def test_one_fine_cell_wide_channel_remains_connected_after_aggregation():
    # A majority-water block rule would erase this narrow resolved corridor.
    sea = np.array([
        [False, False, False, False, False, False, False, False],
        [True, True, True, True, True, True, True, True],
        [True, True, False, False, False, False, True, True],
        [False, False, False, False, False, False, False, False],
    ])

    arrays, audit = m.build_coastal_graph(sea, block_size_cells=4, wrap_longitude=False)

    assert audit["graph_nodes"] == 2
    assert audit["number_of_components"] == 1
    assert len(arrays["edge_i"]) == 1
    assert arrays["boundary_crossing_count"].tolist() == [1]
    assert np.count_nonzero(arrays["fine_cell_to_node"] >= 0) == 12


def test_diagonally_touching_block_corners_do_not_form_a_sea_route():
    sea = np.zeros((4, 4), dtype=bool)
    sea[1, 1] = sea[2, 2] = True

    arrays, audit = m.build_coastal_graph(sea, block_size_cells=2, wrap_longitude=False)

    assert audit["graph_nodes"] == 2
    assert audit["number_of_components"] == 2
    assert len(arrays["edge_i"]) == 0


def test_antimeridian_adjacent_water_cells_retain_their_physical_connection():
    sea = np.zeros((2, 8), dtype=bool)
    sea[0, 0] = sea[0, -1] = True

    wrapped, wrapped_audit = m.build_coastal_graph(sea, block_size_cells=2, wrap_longitude=True)
    _, unwrapped_audit = m.build_coastal_graph(sea, block_size_cells=2, wrap_longitude=False)

    assert wrapped_audit["number_of_components"] == 1
    assert unwrapped_audit["number_of_components"] == 2
    assert len(wrapped["edge_i"]) == 1
    assert wrapped["boundary_crossing_count"].tolist() == [1]


def test_parallel_fine_contacts_are_counted_once_each_on_a_shared_boundary():
    sea = np.ones((2, 4), dtype=bool)

    arrays, audit = m.build_coastal_graph(sea, block_size_cells=2, wrap_longitude=False)

    assert audit["graph_nodes"] == 2
    assert len(arrays["edge_i"]) == 1
    assert arrays["boundary_crossing_count"].tolist() == [2]


def test_dry_blocks_create_neither_nodes_nor_edges():
    arrays, audit = m.build_coastal_graph(np.zeros((4, 4), dtype=bool),
                                         block_size_cells=2, wrap_longitude=True)

    assert np.all(arrays["fine_cell_to_node"] == -1)
    assert audit["graph_nodes"] == 0
    assert audit["number_of_components"] == 0
    assert len(arrays["edge_i"]) == 0


def test_chunked_source_polygon_mask_preserves_north_to_south_row_order():
    # Reversing rows would put this northern land polygon south of the equator.
    polygons = np.array([box(-2, 1, 0, 3)], dtype=object)

    land = m.rasterize_land(polygons, np.array([1.5, 0.5, -0.5]),
                            np.array([-1.5, -0.5, 0.5, 1.5]), chunk_rows=1)

    assert land.tolist() == [[True, True, False, False],
                             [False, False, False, False],
                             [False, False, False, False]]


def test_source_polygon_boundaries_are_land_rather_than_zero_width_sea_routes():
    polygons = np.array([box(0, 0, 1, 1)], dtype=object)

    land = m.rasterize_land(polygons, np.array([0.5]), np.array([0.0, 0.5, 1.0]))

    assert land.tolist() == [[True, True, True]]


def test_changed_coastline_archive_is_rejected_before_geometry_is_used(tmp_path):
    changed = tmp_path / "gshhg-shp-2.3.7.zip"
    changed.write_bytes(b"not the pinned original coastline")

    with pytest.raises(ValueError, match="checksum"):
        m.load_source_polygons(changed)


def test_zero_grid_resolution_is_rejected_before_grid_construction(tmp_path):
    with pytest.raises(ValueError):
        m.run(tmp_path / "not-opened.zip", tmp_path, fine_degrees=0)

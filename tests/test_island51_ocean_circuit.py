"""Real numerical contracts for source-aligned marine circuit resistance."""
from pathlib import Path
import importlib.util
import sys
import warnings

import numpy as np
import pytest
pytest.importorskip("scipy")
import scipy.sparse as ss


def module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "island51_ocean_circuit.py"
    assert path.exists(), "The source-aligned marine circuit implementation is missing"
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("ocean_circuit", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_sparse_solver_reproduces_series_resistors_and_duplicate_terminals():
    m = module()
    lap = ss.csr_matrix([[1., -1., 0.], [-1., 2., -1.], [0., -1., 1.]])
    r, meta = m.effective_resistance(lap, np.array([0, 2, 0]))
    assert np.allclose(r, [[0., 2., 0.], [2., 0., 2.], [0., 2., 0.]], atol=1e-11)
    assert meta["unique_terminal_nodes"] == 2
    assert meta["max_linear_system_relative_residual"] < 1e-10


def test_sparse_solver_reproduces_parallel_paths_without_unrelated_components():
    m = module()
    lap = ss.csr_matrix([[2., -1., 0., -1., 0.], [-1., 2., -1., 0., 0.],
                        [0., -1., 2., -1., 0.], [-1., 0., -1., 2., 0.],
                        [0., 0., 0., 0., 0.]])
    r, meta = m.effective_resistance(lap, np.array([0, 1, 2]))
    assert np.allclose(r, [[0., .75, 1.], [.75, 0., .75], [1., .75, 0.]], atol=1e-11)
    assert meta["solved_component_nodes"] == 4
    with pytest.raises(ValueError, match="disconnected"):
        m.effective_resistance(lap, np.array([0, 4]))


def test_distributed_ports_do_not_short_circuit_their_island_nodes():
    m = module()
    lap = ss.csr_matrix([[1., -1., 0.], [-1., 2., -1.], [0., -1., 1.]])
    # Injecting .5 at each of nodes0,1 and withdrawing1 at node2 dissipates
    # (.5**2)*1 + (1**2)*1 = 1.25. Short-circuiting ports would give1.
    weights = ss.csr_matrix([[.5, .5, 0.], [0., 0., 1.]])
    r, meta = m.distributed_effective_resistance(lap, weights)
    assert r[0, 1] == pytest.approx(1.25, abs=1e-11)
    assert meta["island_terminal_distributions"] == 2
    with pytest.raises(ValueError, match="sum to one"):
        m.distributed_effective_resistance(lap, weights * 2)


def test_numerical_overflow_cannot_pass_the_resistance_accuracy_gate():
    m = module()
    lap = ss.csr_matrix(np.array([[1., -1.], [-1., 1.]]) * 1e-310)
    with warnings.catch_warnings(record=True) as emitted:
        warnings.simplefilter("always")
        with pytest.raises(ValueError):
            m.effective_resistance(lap, np.array([0, 1]))
    assert not emitted, "Nonfinite voltages propagated into resistance arithmetic"


def test_edge_resistance_uses_geodesic_distance_and_mean_cell_resistance():
    m = module()
    # A one-degree equatorial edge has length pi*6371.0088/180 km.
    # Endpoint conductances 2 and 6 have harmonic mean 3, not 4.
    distances, conductance = m.edge_conductance(
        np.array([0., 1.]), np.array([0., 0.]), np.array([0]), np.array([1]),
        np.array([2., 6.]))
    assert distances[0] == pytest.approx(111.1950802335329, rel=1e-10)
    assert conductance[0] == pytest.approx(3. / 111.1950802335329, rel=1e-10)
    with pytest.raises(ValueError, match="positive"):
        m.edge_conductance(np.array([0., 1.]), np.array([0., 0.]),
                           np.array([0]), np.array([1]), np.array([0., 6.]))


def test_ais_grid_rejects_offset_raster_even_when_array_shape_matches(tmp_path):
    m = module()
    p = tmp_path / "offset.npz"
    np.savez(p, AIS_mean=np.ones((170, 360)), source_valid_mask=np.ones((170, 360), bool),
             raster_transform=np.array([1.000083, 0., -180., 0., -1., 85.]))
    with pytest.raises(ValueError, match="aligned"):
        m.load_ais(p)


def test_invalid_ais_cells_remain_outside_all_comparison_models(tmp_path):
    m = module()
    p = tmp_path / "valid.npz"
    arr = np.zeros((170, 360))
    valid = np.ones_like(arr, bool)
    valid[80, 10] = False
    arr[80, 10] = 2147483647.
    arr[80, 11] = 3.
    np.savez(p, AIS_mean=arr, source_valid_mask=valid,
             raster_transform=np.array([1., 0., -180., 0., -1., 85.]))
    values, mask = m.load_ais(p)
    sea = np.ones_like(mask)
    sea[80, 12] = False
    domain, fields = m.conductance_fields(values, mask, np.full_like(values, 2.), sea)
    assert not domain[80, 10] and not domain[80, 12]
    assert domain[80, 9]  # A measured zero remains in the observed domain.
    assert fields["ais"][80, 9] == pytest.approx(1.)
    assert fields["ais"][80, 11] == pytest.approx(1. + np.log(4.))
    assert all(np.isnan(field[80, 10]) for field in fields.values())


def test_resistance_predictor_is_invariant_to_arbitrary_electrical_units():
    m = module()
    r = np.array([[0., 1., 4.], [1., 0., 2.], [4., 2., 0.]])
    a = m.resistance_feature(r, np.array([0, 0, 1]), np.array([1, 2, 2]))
    b = m.resistance_feature(r * 1000., np.array([0, 0, 1]), np.array([1, 2, 2]))
    assert np.allclose(a, b)
    assert np.allclose(a, np.log1p(np.array([.5, 2., 1.])))


def test_port_weights_keep_all_original_ports_without_equalizing_unique_cells():
    m = module()
    assert hasattr(m, "port_weights"), "Original-port distributed terminal builder is missing"
    weights = m.port_weights(["a", "b"], {"a": ["1", "2", "3"], "b": ["4"]},
                             {"1": 0, "2": 0, "3": 1, "4": 2}, 3)
    assert np.allclose(weights.toarray(), [[2/3, 1/3, 0], [0, 0, 1]])
    with pytest.raises(ValueError, match="missing"):
        m.port_weights(["a"], {"a": ["1", "5"]}, {"1": 0}, 2)


def test_actual_port_snap_preserves_nearest_node_and_rejects_distant_water():
    m = module()
    assert hasattr(m, "nearest_sea_cell"), "Source-port marine attachment is missing"
    mapping = np.array([[-1, 0], [-1, 1]], dtype=np.int32)
    # Grid covers lon0..0.25, lat0..0.25. Point is just west of sea cell0.
    result = m.nearest_sea_cell(.12, .1875, mapping, .125,
                               west=0., north=.25, distance_limit_km=10.)
    assert result["node"] == 0 and result["raster_row"] == 0 and result["raster_col"] == 1
    assert 7.4 < result["displacement_km"] < 7.6
    with pytest.raises(ValueError, match="marine"):
        m.nearest_sea_cell(.01, .1875, mapping, .125,
                           west=0., north=.25, distance_limit_km=5.)

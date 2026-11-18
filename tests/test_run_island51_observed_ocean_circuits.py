"""Observed-water geometry, source recovery and fixed island-out comparisons."""
from pathlib import Path
import importlib.util
import json
import hashlib
import sys

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
sparse = pytest.importorskip("scipy.sparse")


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_island51_observed_ocean_circuits.py"


def runner():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("observed_ocean_runner", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def topology(rows, cols, latitudes, longitudes, crossings=4):
    return {
        "node_block_row": np.array(rows),
        "node_block_col": np.array(cols),
        "block_latitude": np.array(latitudes),
        "block_longitude": np.array(longitudes),
        "fine_latitude": np.array([84.9375, 84.8125]),
        "fine_longitude": np.array([-179.9375, -179.8125]),
        "edge_i": np.array([0]),
        "edge_j": np.array([1]),
        "boundary_crossing_count": np.array([crossings]),
    }


def test_finite_volume_east_west_edge_uses_north_south_face_width():
    m = runner()
    graph = topology([0, 0], [0, 1], [0., -.5], [-.25, .25])
    values = m.finite_volume_edges(graph, np.array([1., 1.]))
    assert values["distance_km"][0] == pytest.approx(55.59754011676645, rel=1e-10)
    assert values["face_width_km"][0] == pytest.approx(55.59754011676645, rel=1e-10)
    assert values["conductance"][0] == pytest.approx(1., rel=1e-10)


def test_finite_volume_north_south_face_scales_at_shared_latitude():
    m = runner()
    graph = topology([0, 1], [0, 0], [60.25, 59.75], [0., .5])
    values = m.finite_volume_edges(graph, np.array([1., 1.]))
    # Four .125-degree shared faces at60N have half their equatorial width.
    assert values["face_width_km"][0] == pytest.approx(27.79877005838323, rel=1e-10)
    assert values["conductance"][0] == pytest.approx(.5, rel=1e-10)


def test_dateline_is_adjacent_and_narrow_water_face_reduces_conductance():
    m = runner()
    graph = topology([0, 0], [0, 719], [0., -.5], np.arange(-179.75, 180., .5), crossings=1)
    values = m.finite_volume_edges(graph, np.array([1., 1.]))
    assert values["distance_km"][0] == pytest.approx(55.59754011676645, rel=1e-10)
    assert values["conductance"][0] == pytest.approx(.25, rel=1e-10)


def test_finite_volume_edges_reject_diagonal_or_absent_source_contact():
    m = runner()
    diagonal = topology([0, 1], [0, 1], [.25, -.25], [-.25, .25])
    with pytest.raises(ValueError, match="orthogonal"):
        m.finite_volume_edges(diagonal, np.ones(2))
    closed = topology([0, 0], [0, 1], [0., -.5], [-.25, .25], crossings=0)
    with pytest.raises(ValueError, match="positive.*contact"):
        m.finite_volume_edges(closed, np.ones(2))


def test_recovered_raw_port_centroids_must_reproduce_archived_geodesics():
    m = runner()
    groups = {"1": ["a", "b"], "2": ["c"]}
    coordinates = {"a": (179., 0.), "b": (-179., 0.), "c": (90., 0.)}
    pairs = pd.DataFrame([{"island_a": "1", "island_b": "2", "port_distance_km": 10007.557221017962}])
    receipt = m.verify_centroid_distances(pairs, groups, coordinates)
    assert receipt["max_original_pair_distance_error_km"] < 1e-8
    wrong = pairs.copy()
    wrong["port_distance_km"] += 1.
    with pytest.raises(ValueError, match="archived.*distance"):
        m.verify_centroid_distances(wrong, groups, coordinates)


def test_marine_port_filter_is_explicit_and_cannot_remove_an_island():
    m = runner()
    groups = {"1": ["a", "b"], "2": ["c"]}
    attached = {"a": 0, "c": 1}
    with pytest.raises(ValueError, match="original ports"):
        m.eligible_port_groups(groups, attached, marine_ports_only=False)
    kept = m.eligible_port_groups(groups, attached, marine_ports_only=True)
    assert kept == {"1": ["a"], "2": ["c"]}
    with pytest.raises(ValueError, match="island"):
        m.eligible_port_groups(groups, {"a": 0}, marine_ports_only=True)


def test_original_port_count_control_is_symmetric_in_pair_order():
    m = runner()
    q = pd.DataFrame({"i": [0, 1], "j": [1, 0], "port_distance_km": [30., 30.],
                      "earlier_butterfly_jaccard": [.2, .2]})
    features = m.make_base_features(q, ["1", "2"], {"1": 1, "2": 3}, "new_butterfly_jaccard")
    assert np.allclose(features["sum_log1p_original_port_count"], [np.log(8.), np.log(8.)])
    assert np.allclose(features["abs_difference_log1p_original_port_count"], [np.log(2.), np.log(2.)])
    assert features["earlier_butterfly_Jaccard"].tolist() == [.2, .2]
    earlier = m.make_base_features(q, ["1", "2"], {"1": 1, "2": 3}, "earlier_butterfly_jaccard")
    assert "earlier_butterfly_Jaccard" not in earlier


def test_observed_transport_comparison_keeps_uniform_geometry_and_all_51_holdouts():
    m = runner()
    ids = [str(i) for i in range(51)]
    rows = []
    uniform = np.zeros((51, 51))
    ais = np.zeros((51, 51))
    wind = np.zeros((51, 51))
    for i in range(51):
        for j in range(i+1, 51):
            in_group = i < 7 and j < 7
            distance = float(abs(i-j)+1)
            old = float((i+j) % 11) / 20.
            uniform[i, j] = uniform[j, i] = distance
            wind[i, j] = wind[j, i] = distance ** .8
            ais[i, j] = ais[j, i] = 3. if in_group else 1.
            rows.append({"island_a": str(i), "island_b": str(j), "i": i, "j": j,
                         "port_distance_km": distance, "host_jaccard": .1 + .2*in_group,
                         "earlier_butterfly_jaccard": old,
                         "new_butterfly_jaccard": .15 + .2*in_group + .2*old})
    q = pd.DataFrame(rows)
    result = m.compare_models(q, ids, {key: 1 for key in ids},
                              {"uniform": uniform, "ais": ais, "wind": wind})
    assert set(result) == {"base_vs_uniform", "uniform_vs_ais", "uniform_vs_wind", "uniform_vs_ais_and_wind"}
    for comparison in result.values():
        for model in comparison["targets"].values():
            assert len(model["per_island_LOIO"]) == 51
            assert model["n_unique_observed_pairs"] == 1275
            assert model["n_predictions_heldout_island_pairs"] == 2550
            assert all(row["train_pairs"] == 1225 and row["test_pairs"] == 50 for row in model["per_island_LOIO"])
    for target in ["host_jaccard", "new_butterfly_jaccard"]:
        measured = result["uniform_vs_ais"]["targets"][target]
        assert "log1p_normalized_uniform_circuit" in measured["baseline_features"]
        assert "log1p_normalized_uniform_circuit" in measured["augmented_features"]
        assert measured["MSE_baseline"] > 1e-5
        assert measured["MSE_augmented"] < 1e-20


def test_raw_WPI_checksum_is_verified_before_reusing_recovered_groups(tmp_path):
    m = runner()
    ports = tmp_path / "ports.geojson"
    ports.write_text('{"type":"FeatureCollection","features":[]}')
    pairfile = tmp_path / "pairs.csv"
    pairfile.write_text("island_a,island_b,port_distance_km\n1,2,10\n")
    islandfile = tmp_path / "islands.csv"
    islandfile.write_text("landmass_id,ports,eligible\n1,1,True\n2,1,True\n")
    receipt = tmp_path / "recovery.json"
    receipt.write_text(json.dumps({
        "status": "ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES",
        "source_WPI_port_checksum": "0"*64,
        "source_pair_sha256": hashlib.sha256(pairfile.read_bytes()).hexdigest(),
    }))
    with pytest.raises(ValueError, match="WPI.*checksum"):
        m.validate_source_files(ports, receipt, pairfile, islandfile)


def test_wind_repeats_parent_values_without_shifting_or_inventing_resolution(tmp_path):
    m = runner()
    wind = tmp_path / "wind.npz"
    values = np.arange(170*360, dtype=float).reshape(170, 360) + 1.
    np.savez(wind, mean_monthly_wind_vector_speed_10m=values,
             lat_centers=np.arange(84.5, -85., -1.), lon_centers=np.arange(-179.5, 180., 1.))
    expanded = m.load_wind(wind, (340, 720))
    assert expanded.shape == (340, 720)
    assert np.all(expanded[:2, :2] == 1.)
    assert np.all(expanded[2:4, 2:4] == 362.)
    np.savez(wind, mean_monthly_wind_vector_speed_10m=values,
             lat_centers=np.arange(84.6, -84.9, -1.), lon_centers=np.arange(-179.5, 180., 1.))
    with pytest.raises(ValueError, match="aligned"):
        m.load_wind(wind, (340, 720))


def test_network_export_renumbers_component_once_and_preserves_distributed_weights(tmp_path):
    m = runner()
    graph = {"edge_i": np.array([1, 2]), "edge_j": np.array([2, 3])}
    values = {"uniform": {"conductance": np.array([.5, 1.])}}
    weights = sparse.csr_matrix([[0., 1., 0., 0.], [0., 0., .5, .5]])
    receipt = m.export_networks(tmp_path, graph, values, weights, ["A", "B"], np.array([1, 2, 3]))
    observed = np.loadtxt(tmp_path / "uniform_conductances.txt")
    assert np.array_equal(observed, [[1., 2., .5], [2., 3., 1.]])
    ws = pd.read_csv(tmp_path / "island_node_weights.csv")
    assert ws.loc[ws.landmass_id == "A", "network_node_id"].tolist() == [1]
    assert ws.loc[ws.landmass_id == "B", "network_node_id"].tolist() == [2, 3]
    assert ws.groupby("landmass_id").weight.sum().to_dict() == {"A": 1., "B": 1.}
    assert receipt["nodes"] == 3 and receipt["undirected_edges"] == 2


def test_frozen_manifest_checks_derived_bytes_and_prefit_protocol(tmp_path):
    m = runner()
    protocol = tmp_path / "protocol.md"
    protocol.write_text("Fixed geometry and model comparison before outcomes.\n")
    datafile = tmp_path / "observations.npz"
    np.savez(datafile, observed=np.array([1., 2.]))
    manifest_path = tmp_path / "source_bundle_manifest.json"
    manifest = {
        "protocol_sha256_before_outcome_fit": hashlib.sha256(protocol.read_bytes()).hexdigest(),
        "protocol_frozen_at_utc": "2026-11-06T08:29:53Z",
        "files": [{"path": datafile.name, "bytes": datafile.stat().st_size,
                   "sha256": hashlib.sha256(datafile.read_bytes()).hexdigest()}],
    }
    manifest_path.write_text(json.dumps(manifest))
    receipt = m.validate_source_bundle(tmp_path, protocol)
    assert receipt["manifest_sha256"] == hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert receipt["files_verified_count"] == 1
    protocol.write_text("A changed protocol must not be fitted against this frozen bundle.\n")
    with pytest.raises(ValueError, match="protocol.*frozen"):
        m.validate_source_bundle(tmp_path, protocol)
    protocol.write_text("Fixed geometry and model comparison before outcomes.\n")
    np.savez(datafile, observed=np.array([1., 3.]))
    with pytest.raises(ValueError, match="observations.npz.*checksum"):
        m.validate_source_bundle(tmp_path, protocol)


def test_manifest_cannot_verify_a_file_outside_its_source_bundle(tmp_path):
    m = runner()
    protocol = tmp_path / "protocol.md"
    protocol.write_text("Fixed protocol.\n")
    (tmp_path / "source_bundle_manifest.json").write_text(json.dumps({
        "protocol_sha256_before_outcome_fit": hashlib.sha256(protocol.read_bytes()).hexdigest(),
        "files": [{"path": "../unrelated.npz", "bytes": 0, "sha256": "0"*64}],
    }))
    with pytest.raises(ValueError, match="outside.*bundle"):
        m.validate_source_bundle(tmp_path, protocol)

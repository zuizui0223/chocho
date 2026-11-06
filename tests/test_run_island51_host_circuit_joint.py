"""Fixed joint features, source alignment, island exclusion and loss estimands."""
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("scipy")

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/run_island51_host_circuit_joint.py"


@lru_cache(maxsize=1)
def runner():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("host_circuit_joint", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fixed_joint_models_preserve_host_zero_and_both_signed_interactions():
    m = runner()
    q = pd.DataFrame({"i": [0, 0, 1], "j": [1, 2, 2], "port_distance_km": [10., 30., 20.],
                      "host_jaccard": [0., .25, .75], "earlier_butterfly_jaccard": [.1, .2, .3],
                      "new_butterfly_jaccard": [.4, .5, .6], "shared_new_species": [300, 400, 500],
                      "union_new_species": [600, 700, 800], "new_butterfly_species": [900, 1000, 1100]})
    matrices = {"uniform": np.array([[0., 1., 3.], [1., 0., 2.], [3., 2., 0.]]),
                "ais": np.array([[0., 2., 6.], [2., 0., 4.], [6., 4., 0.]]),
                "wind": np.array([[0., 1., 8.], [1., 0., 3.], [8., 3., 0.]])}
    later = m.make_joint_features(q, ["1", "2", "3"], {"1": 1, "2": 2, "3": 3}, matrices,
                                  "new_butterfly_jaccard")
    assert list(later) == ["G", "H", "C", "HC", "I"]
    assert set(later["H"])-set(later["G"]) == {"host_jaccard"}
    assert set(later["C"])-set(later["G"]) == {
        "log1p_normalized_ais_circuit", "log1p_normalized_wind_circuit"}
    assert set(later["HC"]) == set(later["H"]) | set(later["C"])
    assert set(later["I"])-set(later["HC"]) == {
        "host_times_AIS_connectivity", "host_times_wind_connectivity"}
    np.testing.assert_array_equal(later["H"]["host_jaccard"], [0., .25, .75])
    np.testing.assert_allclose(later["I"]["host_times_AIS_connectivity"],
                               [0., -.25*np.log(2.5), -.75*np.log(2.)])
    np.testing.assert_allclose(later["I"]["host_times_wind_connectivity"],
                               [0., -.25*np.log(11./3.), -.75*np.log(2.)])
    assert "earlier_butterfly_Jaccard" in later["G"]
    earlier = m.make_joint_features(q, ["1", "2", "3"], {"1": 1, "2": 2, "3": 3}, matrices,
                                    "earlier_butterfly_jaccard")
    for target, features in [("new_butterfly_jaccard", later), ("earlier_butterfly_jaccard", earlier)]:
        for model in features.values():
            assert target not in model
            assert not {"shared_new_species", "union_new_species", "new_butterfly_species"}.intersection(model)
    assert "earlier_butterfly_Jaccard" not in earlier["G"]
    with pytest.raises(ValueError, match="butterfly targets"):
        m.make_joint_features(q, ["1", "2", "3"], {"1": 1, "2": 2, "3": 3}, matrices, "host_jaccard")


def test_matrix_island_order_is_not_silently_relabelled(tmp_path):
    m = runner()
    path = tmp_path / "matrix.npz"
    matrix = np.array([[0., 1., 4.], [1., 0., 2.], [4., 2., 0.]])
    np.savez(path, island_ids=["9", "15", "20"], resistance=matrix)
    np.testing.assert_array_equal(m.load_resistance_matrix(path, ["9", "15", "20"]), matrix)
    np.savez(path, island_ids=["15", "9", "20"], resistance=matrix)
    with pytest.raises(ValueError, match="island.*order"):
        m.load_resistance_matrix(path, ["9", "15", "20"])


@pytest.mark.parametrize("corruption", ["nonfinite", "asymmetric", "nonzero_diagonal", "zero_pair"])
def test_matrix_must_be_finite_symmetric_and_distinct(tmp_path, corruption):
    m = runner()
    matrix = np.array([[0., 1., 4.], [1., 0., 2.], [4., 2., 0.]])
    if corruption == "nonfinite":
        matrix[0, 1] = np.nan
    elif corruption == "asymmetric":
        matrix[0, 1] = 2.
    elif corruption == "nonzero_diagonal":
        matrix[0, 0] = .1
    else:
        matrix[0, 1] = matrix[1, 0] = 0.
    path = tmp_path / "matrix.npz"
    np.savez(path, island_ids=["9", "15", "20"], resistance=matrix)
    with pytest.raises(ValueError, match="matrix"):
        m.load_resistance_matrix(path, ["9", "15", "20"])


def locked_fixture(tmp_path):
    original_lock = REPO / "results/island51_host_circuit_joint_20261106/input_lock.json"
    lock = json.loads(original_lock.read_text())
    for index, entry in enumerate(lock["files"]):
        path = tmp_path / entry["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"Independent source fixture {index}\n")
        entry["bytes"] = path.stat().st_size
        entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    protocol = tmp_path / lock["protocol_path"]
    protocol.parent.mkdir(parents=True, exist_ok=True)
    protocol.write_text("Five fixed models, no selection by new outcomes.\n")
    lock["protocol_sha256"] = hashlib.sha256(protocol.read_bytes()).hexdigest()
    implementation = protocol.with_name("IMPLEMENTATION_PLAN.md")
    implementation.write_text("Source gate, fixed models, independent reconstruction.\n")
    lock["implementation_plan_sha256"] = hashlib.sha256(implementation.read_bytes()).hexdigest()
    path = tmp_path / "input_lock.json"
    path.write_text(json.dumps(lock))
    return path, lock


@pytest.mark.parametrize("changed", ["input", "protocol", "implementation_plan"])
def test_changed_frozen_bytes_fail_the_prefit_gate(tmp_path, changed):
    m = runner()
    path, lock = locked_fixture(tmp_path)
    receipt = m.validate_input_lock(path, repo=tmp_path)
    assert receipt["locked_files_verified"] == 12
    if changed == "input":
        altered = tmp_path / lock["files"][0]["path"]
    elif changed == "protocol":
        altered = tmp_path / lock["protocol_path"]
    else:
        altered = (tmp_path / lock["protocol_path"]).with_name("IMPLEMENTATION_PLAN.md")
    altered.write_text("Changed after freeze\n")
    with pytest.raises(ValueError, match="checksum"):
        m.validate_input_lock(path, repo=tmp_path)


def test_port_controls_use_original194_counts_not_marine_counts(tmp_path):
    m = runner()
    ids = [str(index+1) for index in range(51)]
    counts = [4]*41 + [3]*10  #194 original ports.
    path = tmp_path / "islands.csv"
    pd.DataFrame({"landmass_id": ids+["999"], "ports": counts+[20], "eligible": [True]*51+[False]}).to_csv(path, index=False)
    receipt = {"original_ports_audited": 194, "marine_ports_included": 182, "islands_retained": 51,
               "by_island": {key: {"original_ports": count, "marine_ports": count-(index < 12)}
                             for index, (key, count) in enumerate(zip(ids, counts))}}
    observed = m.load_original_port_counts(path, ids, receipt)
    assert sum(observed.values()) == 194
    assert observed["1"] == 4 and receipt["by_island"]["1"]["marine_ports"] == 3
    receipt["by_island"]["1"]["original_ports"] = 3
    with pytest.raises(ValueError, match="original port"):
        m.load_original_port_counts(path, ids, receipt)


def synthetic_source():
    ids = [str(index+1) for index in range(51)]
    matrices = {key: np.zeros((51, 51)) for key in ("uniform", "ais", "wind")}
    rows = []
    for i in range(51):
        for j in range(i+1, 51):
            h = ((i*17+j*7) % 19)/20.
            old = ((i*3+j*11) % 23)/30.
            distance = (j-i+1)*35.+(i*j) % 7
            for key, value in [("uniform", distance/30.),
                               ("ais", distance/20.+((i+j) % 5)/3.),
                               ("wind", np.sqrt(distance)+((i*2+j) % 7))]:
                matrices[key][i, j] = matrices[key][j, i] = value
            rows.append({"island_a": ids[i], "island_b": ids[j], "i": i, "j": j,
                         "port_distance_km": distance, "host_jaccard": h,
                         "earlier_butterfly_jaccard": old,
                         "new_butterfly_jaccard": .07+.25*h+.15*old+.025*((i+j) % 7)})
    return pd.DataFrame(rows), ids, {key: 1+index % 4 for index, key in enumerate(ids)}, matrices


@pytest.fixture(scope="module")
def evaluated_synthetic():
    m = runner()
    q, ids, counts, matrices = synthetic_source()
    predictions, folds, features = m.fit_joint_models(q, ids, counts, matrices)
    return q, ids, counts, matrices, predictions, folds, features


def test_all_predictions_retain_two_endpoint_holdouts_and_individual_losses(evaluated_synthetic):
    m = runner()
    q, ids, _, _, predictions, folds, _ = evaluated_synthetic
    assert len(predictions) == 25500 and len(folds) == 510
    assert predictions.groupby(["target", "model", "source_pair_index"]).size().eq(2).all()
    assert predictions.groupby(["target", "model", "heldout_island_id"]).size().eq(50).all()
    assert ((predictions.heldout_island_id == predictions.island_a) |
            (predictions.heldout_island_id == predictions.island_b)).all()
    assert folds.train_pairs.eq(1225).all() and folds.test_pairs.eq(50).all()
    model_summary, summary, comparison_folds = m.summarize_results(folds)
    assert len(model_summary) == 10 and len(summary) == 10 and len(comparison_folds) == 510
    selected = predictions[(predictions.target == "new_butterfly_jaccard") & (predictions.model == "G")]
    direct_mse = np.mean((selected.observed.to_numpy()-selected.predicted.to_numpy())**2)
    observed = model_summary[(model_summary.target == "new_butterfly_jaccard") & (model_summary.model == "G")].iloc[0]
    assert observed.MSE == pytest.approx(direct_mse, abs=1e-14)
    # Averaging the two endpoint predictions first is a different estimand.
    ensemble = selected.groupby("source_pair_index")[["observed", "predicted"]].mean()
    ensemble_mse = np.mean((ensemble.observed-ensemble.predicted)**2)
    assert abs(direct_mse-ensemble_mse) > 1e-8
    compared = summary[(summary.target == "new_butterfly_jaccard") & (summary.comparison == "H_to_HC")].iloc[0]
    assert compared.relative_MSE_reduction_percent == pytest.approx(
        100*(compared.MSE_reference-compared.MSE_augmented)/compared.MSE_reference, abs=1e-12)
    assert compared.islands_improved+compared.islands_worse+compared.islands_tied == 51


@pytest.mark.parametrize("target", ["new_butterfly_jaccard", "earlier_butterfly_jaccard"])
def test_poisoning_withheld_outcomes_cannot_change_their_predictions(evaluated_synthetic, target):
    m = runner()
    q, ids, counts, matrices, original, _, _ = evaluated_synthetic
    heldout = "17"
    changed = q.copy()
    mask = changed.island_a.eq(heldout) | changed.island_b.eq(heldout)
    changed.loc[mask, target] = 1.-changed.loc[mask, target]
    poisoned, _, _ = m.fit_joint_models(changed, ids, counts, matrices)
    select = lambda p: p[(p.target == target) & (p.heldout_island_id == heldout)].sort_values(["model", "source_pair_index"])
    np.testing.assert_array_equal(select(original).predicted, select(poisoned).predicted)


def test_prior_circuit_gate_checks_every_fold_not_only_global_MSE(evaluated_synthetic):
    m = runner()
    _, _, _, _, _, folds, _ = evaluated_synthetic
    _, summary, compared = m.summarize_results(folds)
    reference = compared[compared.comparison.eq("G_to_C")].copy()
    prior_folds = reference.rename(columns={"reference_mse": "baseline_mse", "augmented_mse": "circuit_augmented_mse",
        "delta_mse_augmented_minus_reference": "delta_mse_circuit_minus_baseline",
        "reference_mae": "MAE_base", "augmented_mae": "MAE_circuit", "improved": "circuit_improves_this_island"})
    prior_folds["comparison"] = "uniform_vs_ais_and_wind"
    prior_summary = summary[summary.comparison.eq("G_to_C")].copy()
    prior_summary["comparison"] = "uniform_vs_ais_and_wind"
    gate = m.verify_prior_circuit_results(folds, summary, prior_folds, prior_summary)
    assert gate["prior_comparison_fold_rows_checked"] == 102
    # Equal and opposite corruption leaves the aggregate unchanged; the gate
    # must still reject the per-island discrepancy.
    prior_folds.iloc[0, prior_folds.columns.get_loc("baseline_mse")] += .001
    prior_folds.iloc[1, prior_folds.columns.get_loc("baseline_mse")] -= .001
    with pytest.raises(ValueError, match="prior.*fold"):
        m.verify_prior_circuit_results(folds, summary, prior_folds, prior_summary)

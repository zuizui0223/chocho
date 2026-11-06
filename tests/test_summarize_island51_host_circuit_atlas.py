"""Synthetic gain accounting and geography; no frozen results are required."""
from copy import deepcopy
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
import shutil

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/summarize_island51_host_circuit_atlas.py"


@lru_cache(maxsize=1)
def module():
    assert SCRIPT.is_file(), "The atlas summary implementation is not present yet"
    spec = importlib.util.spec_from_file_location("synthetic_island_atlas", SCRIPT)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def fold(island, reference, augmented, *, n=2, comparison="C_to_HC", mae=None):
    references = {"C_to_HC": "C", "H_to_HC": "H"}
    mr, ma = mae if mae is not None else (reference, augmented)
    return {"target": "new_butterfly_jaccard", "comparison": comparison,
            "reference_model": references[comparison], "augmented_model": "HC",
            "heldout_island_id": str(island), "train_pairs": (n-1)*(n-2)//2,
            "test_pairs": n-1, "reference_mse": reference, "augmented_mse": augmented,
            "delta_mse_augmented_minus_reference": augmented-reference,
            "reference_mae": mr, "augmented_mae": ma,
            "improved": augmented < reference, "worse": augmented > reference,
            "tied": augmented == reference}


def test_reference_weighted_contributions_are_fourteen_percent_not_mean_island_percent():
    rows, summaries = module().decompose_losses([
        fold("2", .01, .005, mae=(.10, .05)),
        fold("10", .09, .081, mae=(.30, .27))], ["2", "10"])
    assert [r["mse_gain"] for r in rows] == pytest.approx([.005, .009])
    assert [r["mse_relative_gain_percent"] for r in rows] == pytest.approx([50., 10.])
    assert [r["mse_global_gain_contribution_pp"] for r in rows] == pytest.approx([5., 9.])
    assert [r["mse_mean_gain_contribution"] for r in rows] == pytest.approx([.0025, .0045])
    assert summaries[0]["mse_net_global_pp"] == pytest.approx(14.)
    assert summaries[0]["mse_net_mean_gain"] == pytest.approx(.007)
    assert summaries[0]["mae_net_global_pp"] == pytest.approx(20.)
    assert [r["mse_sign"] for r in rows] == ["positive", "positive"]


def test_cancellation_uses_gross_positive_and_adverse_mass_for_concentration():
    rows, summary = module().decompose_losses([
        fold("10", 1., .4, n=4), fold("2", 1., .6, n=4),
        fold("30", .1, 1., n=4), fold("40", .1, .1, n=4)], ["2", "10", "30", "40"])
    s = summary[0]
    assert s["mse_positive_mass"] == pytest.approx(1.)
    assert s["mse_adverse_mass"] == pytest.approx(.9)
    assert s["mse_net_mass"] == pytest.approx(.1)
    assert s["mse_net_mean_gain"] == pytest.approx(.025)
    assert s["positive_top1_share"] == pytest.approx(.6)
    assert s["positive_top3_share"] == pytest.approx(1.)
    assert s["adverse_top5_share"] == pytest.approx(1.)
    assert s["positive_half_mass_count"] == s["adverse_half_mass_count"] == 1
    assert s["ranked_positive_island_ids"] == ["10", "2"]
    assert s["ranked_adverse_island_ids"] == ["30"]
    assert (s["positive_islands"], s["adverse_islands"], s["tied_islands"]) == (2, 1, 1)
    assert next(r for r in rows if r["island_id"] == "30")["mse_sign"] == "negative"


@pytest.mark.parametrize("augmented", [[0., 0.], [.1, .2]])
def test_zero_reference_and_zero_positive_mass_remain_explicitly_undefined(augmented):
    rows, summaries = module().decompose_losses([
        fold("1", 0., augmented[0]), fold("2", 0., augmented[1])], ["1", "2"])
    s = summaries[0]
    assert all(r["mse_relative_gain_percent"] is None for r in rows)
    assert all(r["mse_global_gain_contribution_pp"] is None for r in rows)
    assert s["mse_positive_mass"] == s["positive_islands"] == 0
    assert s["positive_top1_share"] is None
    assert s["positive_top3_share"] is None
    assert s["positive_top5_share"] is None
    assert s["positive_half_mass_count"] is None
    assert s["mse_net_global_pp"] is None
    assert s["ranked_positive_island_ids"] == []
    if augmented == [0., 0.]:
        assert s["adverse_top1_share"] is None
        assert s["adverse_half_mass_count"] is None
        assert s["tied_islands"] == 2
    else:
        assert s["adverse_top1_share"] == pytest.approx(2/3)
        assert s["adverse_half_mass_count"] == 1


def test_concentration_ties_use_numeric_island_ids_and_half_mass_reaches_threshold():
    _, summaries = module().decompose_losses([
        fold("10", 1., .5, n=4), fold("2", 1., .5, n=4),
        fold("30", 1., .5, n=4), fold("40", 1., .5, n=4)], ["2", "10", "30", "40"])
    s = summaries[0]
    assert s["ranked_positive_island_ids"] == ["2", "10", "30", "40"]
    assert s["positive_half_mass_count"] == 2
    assert s["positive_top3_share"] == .75


@pytest.mark.parametrize("damage", ["missing", "duplicate", "unknown", "delta", "sign", "nonfinite", "counts"])
def test_fold_support_or_saved_loss_corruption_is_rejected(damage):
    rows = [fold("1", .1, .05), fold("2", .2, .15)]
    if damage == "missing":
        rows.pop()
    elif damage == "duplicate":
        rows[1] = deepcopy(rows[0])
    elif damage == "unknown":
        rows[1]["heldout_island_id"] = "3"
    elif damage == "delta":
        rows[0]["delta_mse_augmented_minus_reference"] = .05
    elif damage == "sign":
        rows[0]["improved"] = False
    elif damage == "nonfinite":
        rows[0]["reference_mse"] = float("nan")
    else:
        rows[0]["test_pairs"] = 99
    with pytest.raises(ValueError):
        module().decompose_losses(rows, ["1", "2"])


def test_parent_summary_identity_is_checked_with_independent_literal_values():
    m = module()
    _, actual = m.decompose_losses([fold("2", .01, .005), fold("10", .09, .081)], ["2", "10"])
    expected = [{"target": "new_butterfly_jaccard", "comparison": "C_to_HC",
                 "comparison_role": "primary", "reference_model": "C", "augmented_model": "HC",
                 "MSE_reference": .05, "MSE_augmented": .043, "absolute_MSE_gain": .007,
                 "relative_MSE_reduction_percent": 14., "MAE_reference": .05,
                 "MAE_augmented": .043, "islands_improved": 2, "islands_worse": 0, "islands_tied": 0}]
    assert m.verify_parent_summaries(actual, expected)["comparisons_verified"] == 1
    expected[0]["relative_MSE_reduction_percent"] = 30.
    with pytest.raises(ValueError):
        m.verify_parent_summaries(actual, expected)


def test_conditional_categories_keep_all_signs_and_mae_disagreement():
    ids = ["1", "2", "3", "4", "5"]
    pairs = [(.6, .7), (.6, .3), (.3, .6), (.3, .2), (.5, .8)]
    folds = []
    for island, (host_ref, circuit_ref) in zip(ids, pairs):
        folds.extend([fold(island, host_ref, .5, n=5, mae=(.3, .4)),
                      fold(island, circuit_ref, .5, n=5, comparison="H_to_HC", mae=(.6, .4))])
    decomposition, _ = module().decompose_losses(folds, ids)
    status = module().conditional_status(decomposition, ids)
    assert [r["mse_category"] for r in status] == [
        "both_positive", "host_positive_circuit_negative", "host_negative_circuit_positive",
        "both_negative", "at_least_one_tie"]
    assert status[0]["host_mae_sign"] == "negative"
    assert status[0]["circuit_mae_sign"] == "positive"
    assert not any(r["both_mae_positive"] for r in status)
    assert status[0]["joint_mae"] == .4
    damaged = deepcopy(decomposition)
    next(r for r in damaged if r["island_id"] == "1" and r["comparison"] == "H_to_HC")["augmented_mse"] = .4
    with pytest.raises(ValueError, match="HC"):
        module().conditional_status(damaged, ids)
    with pytest.raises(ValueError):
        module().conditional_status(decomposition[:-1], ids)


def geography_fixture():
    features = []
    for port, lon, lat, name, country in [(10, 179., 0., "EAST", "AA"),
                                        (2, -179., 0., "WEST", "BB"),
                                        (30, 20., 10., "OTHER", "AA")]:
        features.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [lon, lat]},
                         "properties": {"OBJECTID": port, "PORT_NAME": name, "COUNTRY": country,
                                        "HARBORSIZE": "S", "HARBORTYPE": "CN"}})
    recovery = {"source_51_port_groups": {"1": ["10", "2"], "20": ["30"]},
                "verified_source_51_mean_port_points": {
                    "1": {"longitude": 180., "latitude": 0., "n_complete_single_landmass_ports": 2},
                    "20": {"longitude": 20., "latitude": 10., "n_complete_single_landmass_ports": 1}}}
    ledger = []
    for feature in features:
        p = feature["properties"]; port = str(p["OBJECTID"])
        lon, lat = feature["geometry"]["coordinates"]; included = port != "2"
        ledger.append({"landmass_id": "20" if port == "30" else "1", "port_id": port,
                       "longitude": lon, "latitude": lat, "source_harbor_type": "CN",
                       "source_port_name": p["PORT_NAME"], "included": included,
                       "node": int(port) if included else "", "raster_row": 1 if included else "",
                       "raster_col": 2 if included else "", "marine_longitude": lon if included else "",
                       "marine_latitude": lat if included else "", "displacement_km": 1. if included else "",
                       "status": "included" if included else "excluded", "reason": "" if included else "outside bound"})
    prep = {"original_ports_audited": 3, "marine_ports_included": 2, "islands_retained": 2,
            "ports_excluded": 1, "by_island": {
                "1": {"original_ports": 2, "marine_ports": 1, "unique_terminal_nodes": 1},
                "20": {"original_ports": 1, "marine_ports": 1, "unique_terminal_nodes": 1}}}
    roster = [{"landmass_id": "1", "ports": 2, "eligible": True},
              {"landmass_id": "20", "ports": 1, "eligible": True}]
    return recovery, {"type": "FeatureCollection", "features": features}, ledger, prep, roster


def test_geography_uses_numeric_objectids_all_original_ports_and_spherical_dateline_mean():
    recovery, raw, ledger, prep, roster = geography_fixture()
    islands, ports, checks = module().build_geography(recovery, raw, ledger, prep, roster, ["1", "20"])
    first = islands[0]
    assert first["representative_longitude"] == 180.
    assert first["representative_latitude"] == 0.
    assert first["example_port_id"] == "2"
    assert first["example_port_name"] == "WEST"
    assert first["source_port_ids"] == ["2", "10"]
    assert first["source_port_names"] == ["WEST", "EAST"]
    assert first["sampled_country_codes"] == ["AA", "BB"]
    assert first["display_label"] == "ID 1 | example port: WEST"
    assert (first["original_ports"], first["marine_ports"], first["excluded_marine_ports"]) == (2, 1, 1)
    excluded = next(r for r in ports if r["port_id"] == "2")
    assert excluded["original_longitude"] == -179.
    assert not excluded["included"]
    for key in ["marine_node", "raster_row", "raster_col", "marine_longitude", "marine_latitude", "displacement_km"]:
        assert excluded[key] is None
    assert checks["all_original_port_groups_spherical_means_verified"]


@pytest.mark.parametrize("damage", ["missing_port", "duplicate_port", "wrong_group", "name", "coordinate", "marine_count", "centroid", "excluded_zero"])
def test_geographic_join_or_original_marine_scope_mismatch_is_rejected(damage):
    recovery, raw, ledger, prep, roster = geography_fixture()
    if damage == "missing_port":
        ledger.pop()
    elif damage == "duplicate_port":
        ledger.append(deepcopy(ledger[0]))
    elif damage == "wrong_group":
        ledger[0]["landmass_id"] = "20"
    elif damage == "name":
        ledger[0]["source_port_name"] = "UNRELATED"
    elif damage == "coordinate":
        ledger[0]["longitude"] = 10.
    elif damage == "marine_count":
        prep["by_island"]["1"]["marine_ports"] = 2
    elif damage == "centroid":
        recovery["verified_source_51_mean_port_points"]["1"]["longitude"] = 0.
    else:
        ledger[1]["marine_longitude"] = 0.
    with pytest.raises(ValueError):
        module().build_geography(recovery, raw, ledger, prep, roster, ["1", "20"])


def test_import_needs_no_result_files_and_does_not_import_the_parent_auditor(tmp_path):
    m = module()
    isolated = tmp_path / "scripts" / SCRIPT.name
    isolated.parent.mkdir()
    shutil.copy2(SCRIPT, isolated)
    spec = importlib.util.spec_from_file_location("isolated_atlas_without_results", isolated)
    copied = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(copied)
    assert copied.decompose_losses([fold("1", .1, .05), fold("2", .2, .1)], ["1", "2"])[1][0]["mse_net_global_pp"] == pytest.approx(50.)
    assert not (tmp_path / "results").exists()


def test_output_guard_rejects_parent_source_and_symlink_output_locations(tmp_path):
    m = module()
    source = tmp_path / "results" / "parent"
    source.mkdir(parents=True)
    locked = source / "losses.csv"
    locked.write_text("source bytes\n")
    output = tmp_path / "results" / "atlas"
    assert m.validate_output_directory(output, [source, tmp_path / "data"]) == output.resolve()
    for rejected in [source, source / "child", source.parent]:
        with pytest.raises(ValueError):
            m.validate_output_directory(rejected, [source, tmp_path / "data"])
    alias = tmp_path / "alias"
    alias.symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError):
        m.validate_output_directory(alias, [source])
    assert locked.read_text() == "source bytes\n"


def test_csv_writer_keeps_nulls_blank_and_arrays_machine_readable(tmp_path):
    path = tmp_path / "example.csv"
    module().write_csv(path, [{"island_id": "2", "undefined": None, "ports": ["2", "10"], "included": False}])
    import csv
    saved = next(csv.DictReader(path.open()))
    assert saved["undefined"] == ""
    assert json.loads(saved["ports"]) == ["2", "10"]
    assert saved["included"] == "False"

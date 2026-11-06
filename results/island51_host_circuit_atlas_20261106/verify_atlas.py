#!/usr/bin/env python3
"""Independent source/geographic accounting audit; never imports atlas code or fits models."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[2]
SPEC = Path(__file__).resolve().parent
PARENT = Path("results/island51_host_circuit_joint_20261106")
CIRCUIT = Path("results/island51_observed_ocean_circuit_20261106")
BUNDLE = Path("data/exploratory/island51_observed_ocean_20261106")
ORIGINAL = Path("provenance/reviewer_defenses/results/island_host_network_20261030")
TARGETS = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
COMPARISONS = {
    "G_to_H": ("G", "H", "contextual"), "G_to_C": ("G", "C", "contextual"),
    "C_to_HC": ("C", "HC", "primary"), "H_to_HC": ("H", "HC", "primary"),
    "HC_to_I": ("HC", "I", "secondary_interaction"),
}
CATEGORIES = ("both_positive", "host_positive_circuit_negative", "host_negative_circuit_positive",
              "both_negative", "at_least_one_tie")
KS = (1, 3, 5)
ATOL = RTOL = 1e-12


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(value, label):
    try:
        answer = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric {label}: {value!r}") from exc
    require(math.isfinite(answer), f"Nonfinite {label}")
    return answer


def integer(value, label):
    answer = number(value, label)
    require(answer.is_integer(), f"Noninteger {label}")
    return int(answer)


def boolean(value, label):
    if isinstance(value, bool):
        return value
    require(str(value).lower() in {"true", "false"}, f"Invalid boolean {label}")
    return str(value).lower() == "true"


def close(actual, expected, label, *, atol=ATOL, rtol=RTOL):
    value = number(actual, label)
    require(math.isclose(value, expected, abs_tol=atol, rel_tol=rtol),
            f"Mismatch {label}: saved={value!r}; independent={expected!r}")
    return abs(value - expected)


def sign(value):
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def indexed(rows, keys, label):
    answer = {}
    for row in rows:
        key = tuple(str(row[column]) for column in keys)
        require(key not in answer, f"Duplicate {label} identity: {key}")
        answer[key] = row
    return answer


def summarize_group(parent_rows):
    """Reconstruct one contrast/target using independent stable summation."""
    require(bool(parent_rows), "Empty comparison group")
    first = parent_rows[0]
    target, comparison = str(first["target"]), str(first["comparison"])
    require(target in TARGETS and comparison in COMPARISONS, "Unexpected target/comparison")
    ref, aug, role = COMPARISONS[comparison]
    n = len(parent_rows)
    rows, seen = [], set()
    for original in parent_rows:
        require(str(original["target"]) == target and str(original["comparison"]) == comparison,
                "Mixed target/comparison group")
        require((original["reference_model"], original["augmented_model"]) == (ref, aug),
                "Parent reference/augmented model direction changed")
        row = dict(original)
        island = str(integer(row["heldout_island_id"], "island ID"))
        require(island not in seen, "Duplicate island in comparison group")
        seen.add(island)
        row.update(island_id=island, comparison_role=role)
        for col in ("train_pairs", "test_pairs"):
            row[col] = integer(row[col], col)
        for metric in ("mse", "mae"):
            reference = number(row[f"reference_{metric}"], f"reference_{metric}")
            augmented = number(row[f"augmented_{metric}"], f"augmented_{metric}")
            require(0 <= reference <= 1 and 0 <= augmented <= 1, "Loss outside the Jaccard range")
            gain = reference - augmented
            row[f"reference_{metric}"], row[f"augmented_{metric}"] = reference, augmented
            row[f"{metric}_gain"] = gain
            row[f"{metric}_relative_gain_percent"] = 100 * gain / reference if reference else None
            row[f"{metric}_mean_gain_contribution"] = gain / n
            row[f"{metric}_sign"] = sign(gain)
        delta = number(row["delta_mse_augmented_minus_reference"], "parent delta")
        close(delta, -row["mse_gain"], "parent delta sign")
        row["delta_mse_augmented_minus_reference"] = delta
        for col, expected in (("improved", row["mse_gain"] > 0),
                              ("worse", row["mse_gain"] < 0), ("tied", row["mse_gain"] == 0)):
            require(boolean(row[col], col) == expected, f"Parent strict {col} flag disagrees")
            row[col] = expected
        rows.append(row)
    summary = dict(target=target, comparison=comparison, comparison_role=role,
                   reference_model=ref, augmented_model=aug, n_islands=n)
    for metric in ("mse", "mae"):
        reference = math.fsum(row[f"reference_{metric}"] for row in rows)
        augmented = math.fsum(row[f"augmented_{metric}"] for row in rows)
        gains = [row[f"{metric}_gain"] for row in rows]
        positive = math.fsum(max(g, 0.) for g in gains)
        adverse = math.fsum(max(-g, 0.) for g in gains)
        net = math.fsum(gains)
        close(net, positive - adverse, f"{metric} mass identity")
        close(net, reference - augmented, f"{metric} parent loss identity")
        for row in rows:
            row[f"{metric}_global_gain_contribution_pp"] = (
                100 * row[f"{metric}_gain"] / reference if reference else None)
        values = dict(reference_sum=reference, augmented_sum=augmented,
                      reference_mean=reference / n, augmented_mean=augmented / n,
                      positive_mass=positive, adverse_mass=adverse, net_mass=net,
                      positive_mean_gain=positive / n, adverse_mean_loss=adverse / n,
                      net_mean_gain=net / n,
                      positive_global_pp=100 * positive / reference if reference else None,
                      adverse_global_pp=100 * adverse / reference if reference else None,
                      net_global_pp=100 * net / reference if reference else None)
        summary.update({f"{metric}_{key}": value for key, value in values.items()})
    positive_rows = sorted((r for r in rows if r["mse_gain"] > 0),
                           key=lambda r: (-r["mse_gain"], int(r["island_id"])))
    adverse_rows = sorted((r for r in rows if r["mse_gain"] < 0),
                          key=lambda r: (r["mse_gain"], int(r["island_id"])))
    summary.update(positive_islands=len(positive_rows), adverse_islands=len(adverse_rows),
                   tied_islands=n-len(positive_rows)-len(adverse_rows))
    for label, ranked in (("positive", positive_rows), ("adverse", adverse_rows)):
        masses = [abs(row["mse_gain"]) for row in ranked]
        total = math.fsum(masses)
        summary[f"ranked_{label}_island_ids"] = [row["island_id"] for row in ranked]
        for k in KS:
            summary[f"{label}_top{k}_share"] = math.fsum(masses[:k]) / total if total else None
        summary[f"{label}_half_mass_count"] = (
            next(k for k in range(1, len(masses)+1) if math.fsum(masses[:k]) >= total/2)
            if total else None)
    return rows, summary


def conditional_status(gains):
    selected = [row for row in gains if row["comparison"] in ("C_to_HC", "H_to_HC")]
    rows = indexed(selected, ("target", "comparison", "island_id"), "conditional source")
    endpoints = {(str(row["target"]), str(row["island_id"])) for row in selected}
    answer = []
    for target, island in sorted(endpoints, key=lambda key: (TARGETS.index(key[0]), int(key[1]))):
        require((target, "C_to_HC", island) in rows and (target, "H_to_HC", island) in rows,
                "Incomplete conditional comparison coverage")
        host, circuit = rows[target, "C_to_HC", island], rows[target, "H_to_HC", island]
        out = dict(target=target, island_id=island)
        for metric in ("mse", "mae"):
            close(circuit[f"augmented_{metric}"], host[f"augmented_{metric}"], f"shared HC {metric}")
            out[f"host_reference_{metric}"] = host[f"reference_{metric}"]
            out[f"circuit_reference_{metric}"] = circuit[f"reference_{metric}"]
            out[f"joint_{metric}"] = host[f"augmented_{metric}"]
            for label, row in (("host", host), ("circuit", circuit)):
                for suffix in ("gain", "sign", "relative_gain_percent", "global_gain_contribution_pp"):
                    out[f"{label}_{metric}_{suffix}"] = row[f"{metric}_{suffix}"]
        h, c = host["mse_sign"], circuit["mse_sign"]
        if "zero" in (h, c):
            category = "at_least_one_tie"
        elif h == c:
            category = "both_positive" if h == "positive" else "both_negative"
        else:
            category = "host_positive_circuit_negative" if h == "positive" else "host_negative_circuit_positive"
        out["mse_category"] = category
        out["both_mae_positive"] = host["mae_gain"] > 0 and circuit["mae_gain"] > 0
        answer.append(out)
    return answer


def spherical_mean(coordinates):
    xyz = []
    for longitude, latitude in coordinates:
        lon, lat = math.radians(longitude), math.radians(latitude)
        xyz.append((math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)))
    require(bool(xyz), "Empty original port group")
    x, y, z = (math.fsum(row[axis] for row in xyz)/len(xyz) for axis in range(3))
    require(math.sqrt(x*x+y*y+z*z) > 1e-12, "Undefined antipodal spherical mean")
    return math.degrees(math.atan2(y, x)), math.degrees(math.atan2(z, math.hypot(x, y)))


def reconstruct_geography(recovery, wpi, ledger, roster, preparation):
    groups = {str(key): sorted(map(str, value), key=int)
              for key, value in recovery["source_51_port_groups"].items()}
    ids = sorted(groups, key=int)
    points = recovery["verified_source_51_mean_port_points"]
    eligible = {str(row["landmass_id"]): row for row in roster if boolean(row["eligible"], "eligible")}
    require(set(eligible) == set(groups) == set(points) == set(preparation["by_island"]),
            "Original geography island coverage differs")
    all_ports = [port for ports in groups.values() for port in ports]
    require(len(all_ports) == len(set(all_ports)), "Original port belongs to multiple groups")
    features = {}
    for feature in wpi["features"]:
        key = str(integer(feature["properties"]["OBJECTID"], "OBJECTID"))
        require(key not in features, "Duplicate source WPI OBJECTID")
        features[key] = feature
    attached = indexed(ledger, ("port_id",), "original port ledger")
    require({key[0] for key in attached} == set(all_ports), "Original port ledger coverage differs")
    require(set(all_ports) <= set(features), "An original port is absent from source WPI")
    islands, ports = [], []
    maximum_centroid_error = maximum_ledger_coordinate_error = 0.
    marine_fields = {"node": "marine_node", "raster_row": "raster_row", "raster_col": "raster_col",
                     "marine_longitude": "marine_longitude", "marine_latitude": "marine_latitude",
                     "displacement_km": "displacement_km"}
    for island in ids:
        current, coords = [], []
        for port in groups[island]:
            feature, source = features[port], attached[port,]
            properties, geometry = feature["properties"], feature["geometry"]
            require(geometry["type"] == "Point", "Original WPI geometry is not a point")
            lon, lat = (number(x, "WPI coordinate") for x in geometry["coordinates"])
            require(-180 <= lon <= 180 and -90 <= lat <= 90, "Original WPI coordinate outside range")
            require(str(source["landmass_id"]) == island, "Ledger island/port identity mismatch")
            require(source["source_port_name"] == properties["PORT_NAME"], "Source port name disagrees with ledger")
            require(source["source_harbor_type"] == properties["HARBORTYPE"], "Source harbor type disagrees")
            for key, value in (("longitude", lon), ("latitude", lat)):
                error = close(source[key], value, f"source port {key}", atol=1e-10, rtol=0)
                maximum_ledger_coordinate_error = max(maximum_ledger_coordinate_error, error)
            included = boolean(source["included"], "marine included")
            out = dict(island_id=island, landmass_id=island, port_id=port,
                       port_name=properties["PORT_NAME"], country_code=properties["COUNTRY"],
                       harbor_size=properties["HARBORSIZE"], harbor_type=properties["HARBORTYPE"],
                       original_longitude=lon, original_latitude=lat, included=included,
                       status=source["status"], reason=source["reason"])
            for original, derived in marine_fields.items():
                if not included:
                    require(source[original] in ("", None), "An excluded marine field must remain blank")
                    out[derived] = None
                elif original in ("node", "raster_row", "raster_col"):
                    out[derived] = integer(source[original], original)
                else:
                    out[derived] = number(source[original], original)
            current.append(out)
            ports.append(out)
            coords.append((lon, lat))
        actual_lon, actual_lat = spherical_mean(coords)
        point = points[island]
        frozen_lon, frozen_lat = number(point["longitude"], "representative longitude"), number(point["latitude"], "representative latitude")
        lon_error = abs((actual_lon-frozen_lon+180) % 360 - 180)
        lat_error = abs(actual_lat-frozen_lat)
        require(max(lon_error, lat_error) <= 1e-10, "Original spherical mean differs from frozen geography")
        maximum_centroid_error = max(maximum_centroid_error, lon_error, lat_error)
        marine = [row for row in current if row["included"]]
        terminals = len({row["marine_node"] for row in marine})
        known = preparation["by_island"][island]
        require(len(current) == integer(eligible[island]["ports"], "original ports") ==
                integer(point["n_complete_single_landmass_ports"], "source mean port count") == known["original_ports"],
                "Original island port counts disagree")
        require(len(marine) == known["marine_ports"] and terminals == known["unique_terminal_nodes"],
                "Marine port/terminal counts disagree")
        example = current[0]
        islands.append(dict(island_id=island, landmass_id=island,
                            representative_longitude=frozen_lon, representative_latitude=frozen_lat,
                            original_ports=len(current), marine_ports=len(marine),
                            excluded_marine_ports=len(current)-len(marine), unique_terminal_nodes=terminals,
                            example_port_id=example["port_id"], example_port_name=example["port_name"],
                            example_port_country_code=example["country_code"],
                            source_port_ids=[row["port_id"] for row in current],
                            source_port_names=[row["port_name"] for row in current],
                            sampled_country_codes=sorted({row["country_code"] for row in current}),
                            display_label=f"ID {island} | example port: {example['port_name']}"))
    return islands, ports, {"islands": len(islands), "original_ports": len(ports),
                            "marine_ports": sum(row["included"] for row in ports),
                            "excluded_ports": sum(not row["included"] for row in ports),
                            "maximum_centroid_error_degrees": maximum_centroid_error,
                            "maximum_ledger_coordinate_error_degrees": maximum_ledger_coordinate_error,
                            "excluded_marine_values_are_missing_not_zero": True}


def compare_table(expected, saved, keys, label):
    reference = indexed(expected, keys, label + " expected")
    actual = indexed(saved, keys, label + " saved")
    require(reference.keys() == actual.keys(), f"{label} row coverage differs")
    maxima = {}
    for key, row in reference.items():
        for field, want in row.items():
            require(field in actual[key], f"Missing {label}.{field}")
            value = actual[key][field]
            if want is None:
                require(value in ("", None), f"{label}.{field} must remain undefined/blank at {key}")
            elif isinstance(want, bool):
                require(boolean(value, field) == want, f"Mismatch {label}.{field} at {key}")
            elif isinstance(want, int):
                require(integer(value, field) == want, f"Mismatch {label}.{field} at {key}")
            elif isinstance(want, float):
                error = close(value, want, f"{label}.{field} at {key}")
                maxima[field] = max(maxima.get(field, 0.), error)
            elif isinstance(want, list):
                values = json.loads(value) if isinstance(value, str) else value
                require(values == want, f"Mismatch ordered {label}.{field} at {key}")
            else:
                require(str(value) == str(want), f"Mismatch {label}.{field} at {key}")
    return {"rows_checked": len(reference), "all_fields_match": True,
            "maximum_absolute_numeric_errors": maxima}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing or duplicate CSV headers: {path}")
        rows = list(reader)
    require(all(None not in row and None not in row.values() for row in rows), f"Malformed CSV: {path}")
    return rows


def verify_manifest(root, entries):
    checked = {}
    for entry in entries:
        path = (root / entry["path"]).resolve()
        require(path.is_relative_to(root.resolve()), "Frozen input path leaves the repository")
        require(entry["path"] not in checked, "Duplicate frozen input path")
        sha = digest(path)
        require(path.stat().st_size == entry["bytes"] and sha == entry["sha256"],
                f"Frozen input size/hash mismatch: {entry['path']}")
        checked[entry["path"]] = sha
    return checked


def verify_lock():
    path = SPEC / "input_lock.json"
    lock = json.loads(path.read_text())
    require(lock["schema"] == "chocho_island51_host_circuit_atlas_input_lock_v1" and
            lock["status"] == "ACCOUNTING_RULES_FIXED_BEFORE_ATLAS_COMPUTATION" and
            lock["new_atlas_computation_already_run"] is False, "Atlas source lock status changed")
    require((lock["fixed_islands"], lock["fixed_pairs"], lock["original_ports"], lock["marine_ports"])
            == (51, 1275, 194, 182), "Frozen atlas support changed")
    require(lock["primary_contrasts"] == ["C_to_HC", "H_to_HC"] and lock["concentration_k"] == list(KS),
            "Frozen primary contrasts or concentration rules changed")
    required = {str(PARENT / name) for name in (
        "heldout_predictions.csv.gz", "model_fold_losses.csv", "comparison_fold_losses.csv",
        "model_summary.csv", "comparison_summary.csv", "joint_run_receipt.json", "source_validation.json",
        "prior_circuit_reproduction.json", "input_lock.json", "independent_prediction_audit.json",
        "PROTOCOL.md", "IMPLEMENTATION_PLAN.md", "verify_saved_predictions.py")}
    required.update(str(p) for p in (
        BUNDLE / "coordinates/island51_original_WPI_ports_coords_unique_source.json",
        BUNDLE / "wpi/ports_2019_source.geojson", CIRCUIT / "original194_port_attachment_audit.csv",
        CIRCUIT / "source_preparation.json", ORIGINAL / "island_host_network_islands.csv",
        ORIGINAL / "island_host_network_pairs.csv", BUNDLE / "coast_0p0625/island51_coastal_ocean_graph.npz"))
    require(len(lock["files"]) == 20 and {entry["path"] for entry in lock["files"]} == required,
            "All 20 fixed atlas inputs must be present in the lock")
    checked = verify_manifest(REPO, lock["files"])
    require(digest(REPO / lock["protocol_path"]) == lock["protocol_sha256"], "Frozen atlas protocol hash mismatch")
    require(digest(SPEC / "IMPLEMENTATION_PLAN.md") == lock["implementation_plan_sha256"],
            "Frozen atlas implementation-plan hash mismatch")
    return dict(input_lock_sha256=digest(path), protocol_sha256=lock["protocol_sha256"],
                implementation_plan_sha256=lock["implementation_plan_sha256"],
                frozen_at_utc=lock["frozen_at_utc"], files_sha256=checked)


def verify_lineage(results_dir, frozen, category_counts, mae_counts, counts):
    receipt = json.loads((results_dir / "atlas_receipt.json").read_text())
    parent = json.loads((results_dir / "parent_prediction_audit.json").read_text())
    require(receipt["schema"] == "chocho_island51_host_circuit_atlas_v1" and
            receipt["status"] == "COMPLETE_POSTHOC_SAVED_PREDICTION_GEOGRAPHIC_ACCOUNTING" and
            receipt["no_models_refitted"] is True, "Atlas completion/no-refit declaration changed")
    for key, value in frozen.items():
        require(receipt["source_validation"][key] == value, f"Atlas recorded source mismatch: {key}")
    require(receipt["row_counts"] == counts, "Atlas receipt row counts disagree")
    require(receipt["conditional_MSE_category_counts"] == category_counts,
            "Atlas receipt joint-sign category counts disagree")
    require(receipt["both_conditional_MAE_positive_counts"] == mae_counts,
            "Atlas receipt MAE sign counts disagree")
    names = set(counts) | {"parent_prediction_audit.json"}
    require(set(receipt["outputs_sha256"]) == names, "Atlas output hash coverage differs")
    for name, sha in receipt["outputs_sha256"].items():
        require(digest(results_dir / name) == sha, f"Atlas output hash mismatch: {name}")
    code_paths = {"scripts/summarize_island51_host_circuit_atlas.py", str(PARENT / "verify_saved_predictions.py")}
    require(set(receipt["implementation_sha256"]) == code_paths, "Atlas implementation hash coverage differs")
    for name, sha in receipt["implementation_sha256"].items():
        require(digest(REPO / name) == sha, f"Atlas implementation hash mismatch: {name}")
    parent_files = receipt["parent_artifacts_verified_unchanged_sha256"]
    required_parent = {name for name in frozen["files_sha256"] if name.startswith(str(PARENT) + "/")}
    require(required_parent <= parent_files.keys(), "Original parent preservation record is incomplete")
    for name, sha in parent_files.items():
        require((REPO / name).resolve().is_relative_to((REPO / PARENT).resolve()),
                "Parent preservation path is outside parent results")
        require(digest(REPO / name) == sha, f"Parent scientific artifact changed: {name}")
    require(parent["status"] == "PASS_SAVED_PREDICTIONS_AND_LOCKED_COMPARISONS_INDEPENDENTLY_RECONSTRUCTED" and
            parent["numerical_audit"]["verdict"] == "PASS", "Fresh parent prediction audit did not pass")
    require(parent["auditor_sha256"] == frozen["files_sha256"][str(PARENT / "verify_saved_predictions.py")],
            "Parent auditor implementation identity changed")
    require(parent["prior_parent_audit_sha256"] == frozen["files_sha256"][str(PARENT / "independent_prediction_audit.json")],
            "Parent audit evidence identity changed")
    require(parent["reconstructed_table_row_counts"] == {
        "model_fold_losses": 510, "comparison_fold_losses": 510, "model_summaries": 10, "comparison_summaries": 10},
        "Fresh parent audit table coverage changed")
    require(parent["numerical_audit"]["prediction_checks"]["rows"] == 25500,
            "Fresh parent audit prediction coverage changed")
    for name, sha in parent["audited_output_sha256"].items():
        require(digest(REPO / PARENT / name) == sha, f"Fresh parent audit file differs: {name}")
    timestamps = [datetime.fromisoformat(text) for text in (
        frozen["frozen_at_utc"], receipt["started_at_utc"],
        receipt["parent_audit_completed_before_accounting_at_utc"], receipt["completed_at_utc"])]
    require(all(value.tzinfo is not None for value in timestamps) and timestamps == sorted(timestamps),
            "Atlas freeze/parent-audit/accounting timestamps disagree")
    return {"all_recorded_output_and_implementation_hashes_match": True,
            "parent_scientific_files_preserved": len(parent_files),
            "fresh_parent_audit_rows": 25500,
            "fresh_parent_audit_and_frozen_auditor_identity_verified": True,
            "declared_parent_audit_precedes_accounting_completion": True,
            "outputs_sha256": {name: digest(results_dir / name) for name in names},
            "atlas_receipt_sha256": digest(results_dir / "atlas_receipt.json"),
            "implementation_sha256": receipt["implementation_sha256"]}


def audit(results_dir):
    frozen = verify_lock()
    islands, ports, geography = reconstruct_geography(
        json.loads((REPO / BUNDLE / "coordinates/island51_original_WPI_ports_coords_unique_source.json").read_text()),
        json.loads((REPO / BUNDLE / "wpi/ports_2019_source.geojson").read_text()),
        read_csv(REPO / CIRCUIT / "original194_port_attachment_audit.csv"),
        read_csv(REPO / ORIGINAL / "island_host_network_islands.csv"),
        json.loads((REPO / CIRCUIT / "source_preparation.json").read_text()))
    require((len(islands), len(ports), geography["marine_ports"]) == (51, 194, 182), "Fixed geographic support changed")
    island_ids = {row["island_id"] for row in islands}
    parent_folds = read_csv(REPO / PARENT / "comparison_fold_losses.csv")
    require(len(parent_folds) == 510, "All 510 parent comparison folds are required")
    groups = defaultdict(list)
    for row in parent_folds:
        require((integer(row["train_pairs"], "train_pairs"), integer(row["test_pairs"], "test_pairs")) == (1225, 50),
                "Parent fold support changed")
        groups[row["target"], row["comparison"]].append(row)
    require(set(groups) == {(target, comparison) for target in TARGETS for comparison in COMPARISONS},
            "Parent comparison/target coverage differs")
    gains, summaries = [], []
    for target in TARGETS:
        for comparison in COMPARISONS:
            rows = groups[target, comparison]
            require(len(rows) == 51 and {str(row["heldout_island_id"]) for row in rows} == island_ids,
                    "A parent comparison lacks complete 51-island support")
            current, summary = summarize_group(rows)
            gains.extend(current)
            summaries.append(summary)
    conditional = conditional_status(gains)
    require(len(conditional) == 102, "Both targets must have all 51 conditional statuses")
    expected = {
        "island_geography.csv": (islands, ("island_id",)),
        "original_port_geography.csv": (ports, ("port_id",)),
        "island_gain_decomposition.csv": (gains, ("target", "comparison", "island_id")),
        "island_conditional_status.csv": (conditional, ("target", "island_id")),
        "gain_concentration_summary.csv": (summaries, ("target", "comparison")),
    }
    table_checks = {name: compare_table(rows, read_csv(results_dir / name), keys, name)
                    for name, (rows, keys) in expected.items()}
    parent_expected = []
    for row in summaries:
        parent_expected.append(dict(
            target=row["target"], comparison=row["comparison"], comparison_role=row["comparison_role"],
            reference_model=row["reference_model"], augmented_model=row["augmented_model"],
            MSE_reference=row["mse_reference_mean"], MSE_augmented=row["mse_augmented_mean"],
            relative_MSE_reduction_percent=row["mse_net_global_pp"], absolute_MSE_gain=row["mse_net_mean_gain"],
            MAE_reference=row["mae_reference_mean"], MAE_augmented=row["mae_augmented_mean"],
            islands_improved=row["positive_islands"], islands_worse=row["adverse_islands"], islands_tied=row["tied_islands"]))
    parent_summary = compare_table(parent_expected, read_csv(REPO / PARENT / "comparison_summary.csv"),
                                   ("target", "comparison"), "parent_summary")
    category_counts = {target: {category: sum(row["target"] == target and row["mse_category"] == category
                                            for row in conditional) for category in CATEGORIES} for target in TARGETS}
    mae_counts = {target: sum(row["target"] == target and row["both_mae_positive"] for row in conditional) for target in TARGETS}
    counts = {name: len(rows) for name, (rows, _) in expected.items()}
    lineage = verify_lineage(results_dir, frozen, category_counts, mae_counts, counts)
    # Source bytes must remain unchanged even while the verification is running.
    require(verify_lock() == frozen, "Atlas frozen sources changed during independent verification")
    return {
        "schema": "chocho_island51_independent_atlas_audit_v1",
        "status": "PASS_INDEPENDENT_SAVED_LOSS_ACCOUNTING_AND_SOURCE_GEOGRAPHY",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(), "auditor_sha256": digest(Path(__file__)),
        "specification_verdict": "PASS_FOR_AUDITED_TABLES_GEOGRAPHY_AND_FROZEN_INPUTS",
        "numerical_quality_verdict": "PASS_NO_DISCREPANCY_BEYOND_FLOATING_POINT_TOLERANCE",
        "absolute_tolerance": ATOL, "relative_tolerance": RTOL,
        "row_counts": counts, "table_checks": table_checks, "parent_global_summary_checks": parent_summary,
        "geography_checks": geography, "source_checks": frozen, "lineage_checks": lineage,
        "conditional_MSE_category_counts": category_counts,
        "both_conditional_MAE_positive_counts": mae_counts,
        "independently_reconstructed_concentration_summaries": summaries,
        "independence": "Stdlib-only independent CSV/JSON, math.fsum and spherical trigonometry; expected values derive from frozen parent folds and original WPI/ledger/port-group sources. No atlas-runner imports, model fitting or predictor selection.",
        "scope_limits": [
            "Audits five numerical/source tables and their provenance, not rendered map placement or label legibility.",
            "The fresh parent prediction audit is hash-linked evidence; this verifier does not rerun a biological model.",
            "Island contributions describe 50 held-out pairs and share observations/HC predictions; they are not deletion refits, independent replications or causal effects.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=SPEC)
    args = parser.parse_args()
    results_dir = args.results_dir.resolve()
    protected = (REPO / PARENT, REPO / CIRCUIT, REPO / "data", REPO / "provenance")
    if any(results_dir.is_relative_to(path.resolve()) for path in protected):
        print("Refusing to write an audit receipt inside parent results or source inputs", file=sys.stderr)
        return 1
    destination = results_dir / "independent_atlas_audit.json"
    try:
        receipt = audit(results_dir)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        receipt = dict(schema="chocho_island51_independent_atlas_audit_v1", status="FAIL",
                       audited_at_utc=datetime.now(timezone.utc).isoformat(), error_type=type(exc).__name__, error=str(exc))
        if results_dir.is_dir():
            destination.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        print(json.dumps(receipt), file=sys.stderr)
        return 1
    destination.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": receipt["status"], "receipt": str(destination),
                      "row_counts": receipt["row_counts"], "geography_checks": receipt["geography_checks"],
                      "conditional_MSE_category_counts": receipt["conditional_MSE_category_counts"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

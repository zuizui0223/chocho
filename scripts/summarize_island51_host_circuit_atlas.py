#!/usr/bin/env python3
"""Describe geographic gains in frozen island predictions without fitting models.

Calculation helpers use only the standard library and synthetic data. The
independent parent auditor is loaded only by the production entry point.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import types


REPO = Path(__file__).resolve().parents[1]
ATLAS = "results/island51_host_circuit_atlas_20261106"
PARENT = "results/island51_host_circuit_joint_20261106"
CIRCUIT = "results/island51_observed_ocean_circuit_20261106"
SOURCE = "provenance/reviewer_defenses/results/island_host_network_20261030"
DATA = "data/exploratory/island51_observed_ocean_20261106"
DEFAULT_LOCK = REPO / ATLAS / "input_lock.json"
TARGETS = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
COMPARISONS = {
    "G_to_H": ("G", "H", "contextual"),
    "G_to_C": ("G", "C", "contextual"),
    "C_to_HC": ("C", "HC", "primary"),
    "H_to_HC": ("H", "HC", "primary"),
    "HC_to_I": ("HC", "I", "secondary_interaction"),
}
CATEGORIES = ("both_positive", "host_positive_circuit_negative", "host_negative_circuit_positive",
              "both_negative", "at_least_one_tie")
INPUT_PATHS = {
    **{key: PARENT + "/" + name for key, name in {
        "predictions": "heldout_predictions.csv.gz", "model_folds": "model_fold_losses.csv",
        "comparison_folds": "comparison_fold_losses.csv", "model_summary": "model_summary.csv",
        "comparison_summary": "comparison_summary.csv", "run_receipt": "joint_run_receipt.json",
        "source_validation": "source_validation.json", "prior_reproduction": "prior_circuit_reproduction.json",
        "parent_lock": "input_lock.json", "previous_parent_audit": "independent_prediction_audit.json",
        "parent_protocol": "PROTOCOL.md", "parent_plan": "IMPLEMENTATION_PLAN.md",
        "parent_auditor": "verify_saved_predictions.py"}.items()},
    "recovery": DATA + "/coordinates/island51_original_WPI_ports_coords_unique_source.json",
    "ports": DATA + "/wpi/ports_2019_source.geojson",
    "attachment": CIRCUIT + "/original194_port_attachment_audit.csv",
    "preparation": CIRCUIT + "/source_preparation.json",
    "islands": SOURCE + "/island_host_network_islands.csv",
    "pairs": SOURCE + "/island_host_network_pairs.csv",
    "coast": DATA + "/coast_0p0625/island51_coastal_ocean_graph.npz",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing or repeated CSV headers: {path}")
        rows = list(reader)
    require(all(None not in row and None not in row.values() for row in rows), f"Malformed CSV: {path}")
    return rows


def number(value, label):
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid number for {label}: {value!r}") from exc
    require(math.isfinite(result), f"Nonfinite number for {label}")
    return result


def integer(value, label):
    result = number(value, label)
    require(result.is_integer() and result >= 0, f"Invalid nonnegative integer for {label}")
    return int(result)


def boolean(value):
    require(str(value).lower() in {"true", "false"}, f"Invalid boolean: {value!r}")
    return str(value).lower() == "true"


def close(actual, expected, label):
    require(math.isclose(number(actual, label), number(expected, label), rel_tol=1e-12, abs_tol=1e-12),
            f"Identity mismatch for {label}: {actual!r} != {expected!r}")


def canonical_ids(island_ids):
    ids = [str(value) for value in island_ids]
    require(len(ids) > 1 and len(ids) == len(set(ids)), "Missing or duplicate canonical island IDs")
    require(all(value == str(integer(value, "island ID")) for value in ids), "Noncanonical island ID")
    require(ids == sorted(ids, key=int), "Canonical island IDs must be in numeric order")
    return ids


def percentage(gain, reference):
    if reference == 0:
        return None
    return number(100 * (gain / reference), "derived percentage")


def sign(value):
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def decompose_losses(folds, island_ids):
    """Retain each saved fold and account for signed gains with common denominators."""
    ids = canonical_ids(island_ids)
    groups = defaultdict(dict)
    n = len(ids)
    for source in folds:
        target, comparison = str(source["target"]), str(source["comparison"])
        require(target in TARGETS and comparison in COMPARISONS, "Unknown target or comparison")
        reference, augmented, role = COMPARISONS[comparison]
        require((source["reference_model"], source["augmented_model"]) == (reference, augmented),
                "Comparison model labels are misaligned")
        island = str(source["heldout_island_id"])
        require(island in ids and island not in groups[target, comparison], "Unknown or duplicate fold island")
        row = dict(source, island_id=island, comparison_role=role)
        for field, expected in (("train_pairs", (n-1)*(n-2)//2), ("test_pairs", n-1)):
            row[field] = integer(source[field], field)
            require(row[field] == expected, f"Incorrect held-island {field}")
        for metric in ("mse", "mae"):
            for prefix in ("reference", "augmented"):
                field = f"{prefix}_{metric}"
                row[field] = number(source[field], field)
                require(row[field] >= 0, f"Negative loss: {field}")
            gain = row[f"reference_{metric}"] - row[f"augmented_{metric}"]
            row[f"{metric}_gain"] = gain
            row[f"{metric}_relative_gain_percent"] = percentage(gain, row[f"reference_{metric}"])
            row[f"{metric}_mean_gain_contribution"] = gain / n
            row[f"{metric}_sign"] = sign(gain)
        delta = "delta_mse_augmented_minus_reference"
        row[delta] = number(source[delta], delta)
        close(row[delta], -row["mse_gain"], "opposite saved delta sign")
        for field, expected in (("improved", row["mse_gain"] > 0), ("worse", row["mse_gain"] < 0),
                                ("tied", row["mse_gain"] == 0)):
            row[field] = boolean(source[field])
            require(row[field] == expected, f"Saved strict {field} sign disagrees")
        groups[target, comparison][island] = row
    require(groups, "No comparison folds")
    all_rows, summaries = [], []
    for target in TARGETS:
        for comparison, (reference, augmented, role) in COMPARISONS.items():
            if (target, comparison) not in groups:
                continue
            group = groups[target, comparison]
            require(set(group) == set(ids), "Incomplete canonical island support in a comparison")
            rows = [group[key] for key in ids]
            summary = {"target": target, "comparison": comparison, "comparison_role": role,
                       "reference_model": reference, "augmented_model": augmented, "n_islands": n}
            for metric in ("mse", "mae"):
                r = math.fsum(row[f"reference_{metric}"] for row in rows)
                a = math.fsum(row[f"augmented_{metric}"] for row in rows)
                positive = math.fsum(max(row[f"{metric}_gain"], 0.) for row in rows)
                adverse = math.fsum(max(-row[f"{metric}_gain"], 0.) for row in rows)
                net = positive - adverse
                close(net, math.fsum(row[f"{metric}_gain"] for row in rows), metric + " gross-to-net")
                close(net, r-a, metric + " sum-to-net")
                for row in rows:
                    row[f"{metric}_global_gain_contribution_pp"] = percentage(row[f"{metric}_gain"], r)
                summary.update({f"{metric}_reference_sum": r, f"{metric}_augmented_sum": a,
                                f"{metric}_reference_mean": r/n, f"{metric}_augmented_mean": a/n,
                                f"{metric}_positive_mass": positive, f"{metric}_adverse_mass": adverse,
                                f"{metric}_net_mass": net, f"{metric}_positive_mean_gain": positive/n,
                                f"{metric}_adverse_mean_loss": adverse/n, f"{metric}_net_mean_gain": net/n,
                                f"{metric}_positive_global_pp": percentage(positive, r),
                                f"{metric}_adverse_global_pp": percentage(adverse, r),
                                f"{metric}_net_global_pp": percentage(net, r)})
                close(math.fsum(row[f"{metric}_mean_gain_contribution"] for row in rows), net/n,
                      metric + " mean contributions")
                if r:
                    close(math.fsum(row[f"{metric}_global_gain_contribution_pp"] for row in rows),
                          summary[f"{metric}_net_global_pp"], metric + " percentage-point contributions")
            for direction, multiplier in (("positive", 1), ("adverse", -1)):
                ranked = sorted(((multiplier*row["mse_gain"], row["island_id"]) for row in rows
                                 if multiplier*row["mse_gain"] > 0), key=lambda item: (-item[0], int(item[1])))
                mass = summary[f"mse_{direction}_mass"]
                summary[f"{direction}_islands"] = len(ranked)
                summary[f"ranked_{direction}_island_ids"] = [island for _, island in ranked]
                for k in (1, 3, 5):
                    summary[f"{direction}_top{k}_share"] = math.fsum(v for v, _ in ranked[:k])/mass if mass else None
                summary[f"{direction}_half_mass_count"] = next(
                    (k for k in range(1, len(ranked)+1)
                     if math.fsum(v for v, _ in ranked[:k]) >= mass/2), None)
            summary["tied_islands"] = sum(row["mse_gain"] == 0 for row in rows)
            require(summary["positive_islands"] + summary["adverse_islands"] + summary["tied_islands"] == n,
                    "MSE sign counts do not cover every island")
            all_rows.extend(rows)
            summaries.append(summary)
    return all_rows, summaries


def verify_parent_summaries(summaries, parent):
    lookup = {}
    for row in parent:
        key = (row["target"], row["comparison"])
        require(key not in lookup, "Duplicate parent summary")
        lookup[key] = row
    require(len(lookup) == len(summaries), "Parent summary coverage changed")
    maxima = {}
    fields = {"mse_reference_mean": "MSE_reference", "mse_augmented_mean": "MSE_augmented",
              "mse_net_mean_gain": "absolute_MSE_gain", "mse_net_global_pp": "relative_MSE_reduction_percent",
              "mae_reference_mean": "MAE_reference", "mae_augmented_mean": "MAE_augmented",
              "positive_islands": "islands_improved", "adverse_islands": "islands_worse", "tied_islands": "islands_tied"}
    for row in summaries:
        key = (row["target"], row["comparison"])
        require(key in lookup, "Parent comparison missing")
        old = lookup[key]
        for field in ("comparison_role", "reference_model", "augmented_model"):
            require(row[field] == old[field], "Parent comparison metadata changed")
        for actual_key, saved_key in fields.items():
            close(row[actual_key], old[saved_key], saved_key)
            maxima[saved_key] = max(maxima.get(saved_key, 0.), abs(row[actual_key]-float(old[saved_key])))
    return {"comparisons_verified": len(summaries), "absolute_tolerance": 1e-12,
            "relative_tolerance": 1e-12, "maximum_absolute_errors": maxima}


def conditional_status(decomposition, island_ids):
    ids = canonical_ids(island_ids)
    lookup = {}
    targets = {row["target"] for row in decomposition}
    for row in decomposition:
        if row["comparison"] not in {"C_to_HC", "H_to_HC"}:
            continue
        key = (row["target"], row["island_id"], row["comparison"])
        require(key not in lookup, "Duplicate conditional comparison")
        lookup[key] = row
    result = []
    for target in TARGETS:
        if target not in targets:
            continue
        for island in ids:
            keys = [(target, island, comparison) for comparison in ("C_to_HC", "H_to_HC")]
            require(all(key in lookup for key in keys), "Missing conditional island/target join")
            host, circuit = (lookup[key] for key in keys)
            row = {"target": target, "island_id": island}
            for metric in ("mse", "mae"):
                close(host[f"augmented_{metric}"], circuit[f"augmented_{metric}"], "shared HC " + metric)
                row[f"host_reference_{metric}"] = host[f"reference_{metric}"]
                row[f"circuit_reference_{metric}"] = circuit[f"reference_{metric}"]
                row[f"joint_{metric}"] = host[f"augmented_{metric}"]
                for label, source in (("host", host), ("circuit", circuit)):
                    for suffix in ("gain", "sign", "relative_gain_percent", "global_gain_contribution_pp"):
                        row[f"{label}_{metric}_{suffix}"] = source[f"{metric}_{suffix}"]
            hs, cs = row["host_mse_sign"], row["circuit_mse_sign"]
            row["mse_category"] = ("at_least_one_tie" if "zero" in (hs, cs) else
                                   "both_positive" if hs == cs == "positive" else
                                   "both_negative" if hs == cs == "negative" else
                                   "host_positive_circuit_negative" if hs == "positive" else
                                   "host_negative_circuit_positive")
            row["both_mae_positive"] = row["host_mae_sign"] == row["circuit_mae_sign"] == "positive"
            result.append(row)
    require(len(lookup) == 2*len(result), "Unexpected conditional island/target rows")
    return result


def build_geography(recovery, raw_ports, ledger, preparation, roster, island_ids):
    """Join original OBJECTIDs and reproduce means using every original port."""
    ids = canonical_ids(island_ids)
    groups = {str(key): [str(integer(port, "group port ID")) for port in ports]
              for key, ports in recovery["source_51_port_groups"].items()}
    points = recovery["verified_source_51_mean_port_points"]
    require(set(groups) == set(points) == set(ids), "Geography island coverage differs from canonical IDs")
    all_ports = [port for island in ids for port in groups[island]]
    require(len(all_ports) == len(set(all_ports)), "Original port belongs to multiple island groups")
    require(raw_ports.get("type") == "FeatureCollection", "WPI is not a FeatureCollection")
    features = {}
    for feature in raw_ports["features"]:
        port = str(integer(feature["properties"]["OBJECTID"], "WPI OBJECTID"))
        require(port not in features, "Duplicate WPI OBJECTID")
        features[port] = feature
    require(set(all_ports) <= set(features), "Original port missing in WPI")
    attached = {}
    for row in ledger:
        port = str(integer(row["port_id"], "ledger port ID"))
        require(port not in attached, "Duplicate port ledger entry")
        attached[port] = row
    require(set(attached) == set(all_ports), "Port ledger is not an exact one-to-one original port join")
    eligible = {}
    for row in roster:
        if boolean(row["eligible"]):
            key = str(row["landmass_id"])
            require(key not in eligible, "Duplicate eligible island roster row")
            eligible[key] = integer(row["ports"], "roster port count")
    require(set(eligible) == set(preparation["by_island"]) == set(ids), "Roster/preparation island mismatch")
    islands, port_rows, max_point_error = [], [], 0.
    for island in ids:
        port_ids = sorted(groups[island], key=int)
        require(port_ids, "Empty original island port group")
        island_ports, xyz = [], []
        for port in port_ids:
            feature, old = features[port], attached[port]
            require(str(old["landmass_id"]) == island, "Port ledger landmass disagrees with original group")
            geometry, properties = feature["geometry"], feature["properties"]
            require(geometry.get("type") == "Point" and len(geometry["coordinates"]) == 2, "Invalid WPI point")
            lon, lat = [number(v, "WPI coordinate") for v in geometry["coordinates"]]
            require(-180 <= lon <= 180 and -90 <= lat <= 90, "WPI coordinate outside geographic bounds")
            close(old["longitude"], lon, "ledger source longitude")
            close(old["latitude"], lat, "ledger source latitude")
            require(old["source_port_name"] == properties["PORT_NAME"] and
                    old["source_harbor_type"] == properties["HARBORTYPE"], "Source port name/type mismatch")
            included = boolean(old["included"])
            row = {"island_id": island, "landmass_id": island, "port_id": port,
                   "port_name": properties["PORT_NAME"], "country_code": properties["COUNTRY"],
                   "harbor_size": properties["HARBORSIZE"], "harbor_type": properties["HARBORTYPE"],
                   "original_longitude": lon, "original_latitude": lat, "included": included}
            for destination, original in (("marine_node", "node"), ("raster_row", "raster_row"),
                                          ("raster_col", "raster_col"), ("marine_longitude", "marine_longitude"),
                                          ("marine_latitude", "marine_latitude"), ("displacement_km", "displacement_km")):
                value = old[original]
                if included:
                    row[destination] = (integer(value, destination) if destination in {"marine_node", "raster_row", "raster_col"}
                                        else number(value, destination))
                else:
                    require(value is None or value == "", "Excluded port has a fabricated marine coordinate/node/displacement")
                    row[destination] = None
            if included:
                require(-180 <= row["marine_longitude"] <= 180 and -90 <= row["marine_latitude"] <= 90 and
                        0 <= row["displacement_km"] <= 25, "Marine attachment outside frozen geometry bounds")
            row.update(status=old["status"], reason=old["reason"])
            island_ports.append(row)
            ll, aa = math.radians(lon), math.radians(lat)
            xyz.append((math.cos(aa)*math.cos(ll), math.cos(aa)*math.sin(ll), math.sin(aa)))
        x, y, z = [math.fsum(v[j] for v in xyz)/len(xyz) for j in range(3)]
        require(math.sqrt(x*x+y*y+z*z) >= 1e-7, "Original port spherical mean is undefined")
        mean_lon = math.degrees(math.atan2(y, x))
        mean_lat = math.degrees(math.atan2(z, math.hypot(x, y)))
        point = points[island]
        lon, lat = number(point["longitude"], "representative longitude"), number(point["latitude"], "representative latitude")
        require(-180 <= lon <= 180 and -90 <= lat <= 90, "Representative point outside geographic bounds")
        error = max(abs((mean_lon-lon+180) % 360-180), abs(mean_lat-lat))
        require(error <= 1e-10, "Original all-port spherical mean does not reproduce saved representative point")
        max_point_error = max(max_point_error, error)
        require(integer(point["n_complete_single_landmass_ports"], "representative point port count") == len(port_ids),
                "Representative point original port count mismatch")
        marine = sum(row["included"] for row in island_ports)
        unique_nodes = len({row["marine_node"] for row in island_ports if row["included"]})
        declared = preparation["by_island"][island]
        require(eligible[island] == declared["original_ports"] == len(port_ids) and
                marine == declared["marine_ports"] and unique_nodes == declared["unique_terminal_nodes"],
                "Original/marine island port counts differ from source preparation")
        example = island_ports[0]
        islands.append({"island_id": island, "landmass_id": island,
                        "representative_longitude": lon, "representative_latitude": lat,
                        "original_ports": len(port_ids), "marine_ports": marine,
                        "excluded_marine_ports": len(port_ids)-marine, "unique_terminal_nodes": unique_nodes,
                        "example_port_id": example["port_id"], "example_port_name": example["port_name"],
                        "example_port_country_code": example["country_code"], "source_port_ids": port_ids,
                        "source_port_names": [row["port_name"] for row in island_ports],
                        "sampled_country_codes": sorted({row["country_code"] for row in island_ports}),
                        "display_label": f"ID {island} | example port: {example['port_name']}"})
        port_rows.extend(island_ports)
    marine_total = sum(row["included"] for row in port_rows)
    require((len(port_rows), marine_total, len(islands), len(port_rows)-marine_total) ==
            (preparation["original_ports_audited"], preparation["marine_ports_included"],
             preparation["islands_retained"], preparation["ports_excluded"]), "Global original/marine scope mismatch")
    return islands, port_rows, {
        "islands": len(islands), "original_ports": len(port_rows), "marine_ports": marine_total,
        "excluded_marine_ports": len(port_rows)-marine_total,
        "all_original_port_groups_spherical_means_verified": True,
        "maximum_representative_coordinate_error_degrees": max_point_error,
        "exact_OBJECTID_ledger_join_and_source_names_coordinates_verified": True,
        "source_labels": "Canonical landmass ID and smallest-OBJECTID example port; country sets describe sampled WPI ports",
        "representative_point": "Saved spherical mean of all original complete-port coordinates; marine exclusions retained",
    }


def validate_output_directory(out_dir, protected_roots):
    out = Path(out_dir).resolve()
    for root in protected_roots:
        protected = Path(root).resolve()
        require(not out.is_relative_to(protected) and not protected.is_relative_to(out),
                f"Output would overlap a protected parent/source location: {out}")
    require(not out.exists() or out.is_dir(), "Output location is not a directory")
    return out


def write_csv(path, rows):
    require(rows, "Cannot write an empty atlas table")
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, separators=(",", ":"))
                             if isinstance(value, (list, dict)) else value for key, value in row.items()})


def validate_input_lock(lock_path, *, repo=REPO):
    root, lock_path = Path(repo).resolve(), Path(lock_path).resolve()
    lock = json.loads(lock_path.read_text())
    require(lock.get("schema") == "chocho_island51_host_circuit_atlas_input_lock_v1" and
            lock.get("status") == "ACCOUNTING_RULES_FIXED_BEFORE_ATLAS_COMPUTATION" and
            lock.get("existing_fold_outputs_and_aggregate_results_already_available") is True and
            lock.get("new_atlas_computation_already_run") is False, "Atlas post-hoc accounting lock is invalid")
    require((lock.get("fixed_islands"), lock.get("fixed_pairs"), lock.get("original_ports"), lock.get("marine_ports")) ==
            (51, 1275, 194, 182), "Frozen atlas support changed")
    require(lock.get("primary_contrasts") == ["C_to_HC", "H_to_HC"] and lock.get("concentration_k") == [1, 3, 5],
            "Fixed contrasts or concentration rules changed")
    files = lock.get("files", [])
    require(len(files) == 20 and {entry["path"] for entry in files} == set(INPUT_PATHS.values()),
            "Exactly all 20 fixed atlas inputs must be pinned")
    verified = {}
    for entry in files:
        path = (root / entry["path"]).resolve()
        require(path.is_relative_to(root) and path.is_file(), "Locked atlas input is missing or outside repository")
        require(path.stat().st_size == entry["bytes"] and digest(path) == entry["sha256"],
                "Locked atlas input checksum/size changed: " + entry["path"])
        verified[entry["path"]] = entry["sha256"]
    require(lock.get("protocol_path") == ATLAS + "/PROTOCOL.md", "Atlas specification path changed")
    protocol = root / lock["protocol_path"]
    require(digest(protocol) == lock["protocol_sha256"] and
            digest(protocol.with_name("IMPLEMENTATION_PLAN.md")) == lock["implementation_plan_sha256"],
            "Frozen atlas protocol/plan checksum changed")
    frozen = datetime.fromisoformat(lock["frozen_at_utc"])
    require(frozen.tzinfo is not None and frozen <= datetime.now(timezone.utc), "Invalid atlas freeze timestamp")
    return {"status": "ALL_20_ATLAS_INPUTS_AND_FROZEN_ACCOUNTING_RULES_VERIFIED",
            "input_lock_sha256": digest(lock_path), "protocol_sha256": lock["protocol_sha256"],
            "implementation_plan_sha256": lock["implementation_plan_sha256"],
            "frozen_at_utc": lock["frozen_at_utc"], "parent_commit": lock["parent_commit"],
            "files_sha256": verified, "existing_joint_results_already_known": True,
            "accounting_rules_fixed_before_atlas_computation": True}


def run(out_dir, input_lock=DEFAULT_LOCK):
    started = datetime.now(timezone.utc).isoformat()
    validation = validate_input_lock(input_lock)
    out = validate_output_directory(out_dir, [REPO / path for path in (PARENT, CIRCUIT, SOURCE, DATA, "scripts", "tests")])
    paths = {key: REPO / value for key, value in INPUT_PATHS.items()}
    parent_hashes = {str(path.relative_to(REPO)): digest(path) for path in (REPO/PARENT).rglob("*")
                     if path.is_file() and "__pycache__" not in path.parts}
    # Compiling the pinned auditor directly also avoids creating a parent __pycache__.
    auditor = types.ModuleType("frozen_parent_prediction_auditor")
    auditor.__file__ = str(paths["parent_auditor"])
    exec(compile(paths["parent_auditor"].read_bytes(), auditor.__file__, "exec"), auditor.__dict__)
    parent_audit = auditor.audit(REPO / PARENT)
    require(parent_audit["status"] == "PASS_SAVED_PREDICTIONS_AND_LOCKED_COMPARISONS_INDEPENDENTLY_RECONSTRUCTED",
            "Independent parent prediction audit did not pass")
    audited_at = datetime.now(timezone.utc).isoformat()
    roster = read_csv(paths["islands"])
    ids = sorted([row["landmass_id"] for row in roster if boolean(row["eligible"])], key=int)
    require(len(ids) == 51, "Production atlas requires exactly 51 islands")
    recovery = json.loads(paths["recovery"].read_text())
    require(recovery["status"] == "ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES" and
            recovery["source_WPI_port_checksum"] == digest(paths["ports"]) and
            recovery["source_pair_sha256"] == digest(paths["pairs"]), "Original coordinate recovery provenance mismatch")
    geography, ports, geography_checks = build_geography(
        recovery, json.loads(paths["ports"].read_text()), read_csv(paths["attachment"]),
        json.loads(paths["preparation"].read_text()), roster, ids)
    require((len(ports), sum(row["included"] for row in ports)) == (194, 182), "Production port scope changed")
    folds = read_csv(paths["comparison_folds"])
    decomposition, concentration = decompose_losses(folds, ids)
    require(len(decomposition) == 510 and len(concentration) == 10, "All five contrasts/two targets must be retained")
    summary_checks = verify_parent_summaries(concentration, read_csv(paths["comparison_summary"]))
    conditional = conditional_status(decomposition, ids)
    require(len(conditional) == 102, "Both conditional targets must retain every island")
    category_counts = {target: {category: sum(row["target"] == target and row["mse_category"] == category
                                              for row in conditional) for category in CATEGORIES} for target in TARGETS}
    fresh = validate_input_lock(input_lock)
    require(fresh == validation, "Atlas inputs changed during accounting")
    require(all(digest(REPO/path) == sha for path, sha in parent_hashes.items()), "Parent artifacts were modified")
    compact_audit = {key: value for key, value in parent_audit.items() if not key.startswith("reconstructed_")}
    compact_audit["reconstructed_table_row_counts"] = {
        key.removeprefix("reconstructed_"): len(value) for key, value in parent_audit.items() if key.startswith("reconstructed_")}
    compact_audit["prior_parent_audit_sha256"] = digest(paths["previous_parent_audit"])
    tables = {"island_geography.csv": geography, "original_port_geography.csv": ports,
              "island_gain_decomposition.csv": decomposition, "island_conditional_status.csv": conditional,
              "gain_concentration_summary.csv": concentration}
    output_names = [*tables, "parent_prediction_audit.json", "atlas_receipt.json"]
    locked_paths = {path.resolve() for path in paths.values()} | {Path(input_lock).resolve(),
                    (REPO/ATLAS/"PROTOCOL.md").resolve(), (REPO/ATLAS/"IMPLEMENTATION_PLAN.md").resolve()}
    for name in output_names:
        require(not (out/name).is_symlink() and (out/name).resolve() not in locked_paths,
                "Output file overlaps a frozen input or is a symlink")
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        write_csv(out/name, rows)
    (out/"parent_prediction_audit.json").write_text(json.dumps(compact_audit, indent=2, allow_nan=False) + "\n")
    require(validate_input_lock(input_lock) == validation and
            all(digest(REPO/path) == sha for path, sha in parent_hashes.items()), "Source/parent changed during export")
    receipt = {
        "schema": "chocho_island51_host_circuit_atlas_v1",
        "status": "COMPLETE_POSTHOC_SAVED_PREDICTION_GEOGRAPHIC_ACCOUNTING",
        "started_at_utc": started, "parent_audit_completed_before_accounting_at_utc": audited_at,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(), "source_validation": validation,
        "geography_checks": geography_checks, "parent_summary_checks": summary_checks,
        "row_counts": {name: len(rows) for name, rows in tables.items()},
        "conditional_MSE_category_counts": category_counts,
        "both_conditional_MAE_positive_counts": {target: sum(row["target"] == target and row["both_mae_positive"]
                                                             for row in conditional) for target in TARGETS},
        "parent_artifacts_verified_unchanged_sha256": parent_hashes,
        "outputs_sha256": {name: digest(out/name) for name in output_names if name != "atlas_receipt.json"},
        "implementation_sha256": {str(Path(__file__).resolve().relative_to(REPO)): digest(__file__),
                                  INPUT_PATHS["parent_auditor"]: digest(paths["parent_auditor"])},
        "python_version": platform.python_version(), "no_models_refitted": True,
        "undefined_values": "Blank CSV cells / JSON null; no denominator threshold or invented zero",
        "concentration": "MSE only; k=1,3,5; gross positive/adverse masses separate; numeric island ID breaks exact ties",
        "interpretation": "Post-hoc geography of 50 endpoint-held-out pair predictions per island. Shared HC and overlapping folds; no independent replication, causal transport, abundance or establishment inference.",
    }
    (out/"atlas_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--input-lock", type=Path, default=DEFAULT_LOCK)
    args = parser.parse_args()
    result = run(args.out_dir, args.input_lock)
    print(json.dumps({"status": result["status"], "row_counts": result["row_counts"],
                      "conditional_MSE_category_counts": result["conditional_MSE_category_counts"]}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the frozen observed-water circuit comparison on the original51 islands.

Raw source rasters and exported large network files stay outside the repository.
The full194-port ledger is retained even in the explicit182-port marine analysis.
No biology is used to select coast geometry, marine ports or conductance fields.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy import sparse

import island51_ocean_circuit as circuit
from exploratory_island51_circuit_network_holdout import source, make_features, holdout
from reconstruct_island51_verified_port_coordinates import spherical_mean, arc_distance


REPO = Path(__file__).resolve().parents[1]
SOURCE_DIRECTORY = REPO / "provenance/reviewer_defenses/results/island_host_network_20261030"
PROTOCOL = REPO / "results/island51_observed_ocean_circuit_20261106/PROTOCOL.md"
AIS_SOURCE_SHA = "0cfdce41d9934982bbc8f16bd92ce7545cc088839327cbaa9b0ed1469dda31e8"
GSHHG_SOURCE_SHA = "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
RECOVERY_STATUS = "ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4*1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_source_bundle(input_root, protocol_path=PROTOCOL):
    """Verify frozen derived inputs and protocol before reading any outcomes."""
    root = Path(input_root).resolve()
    manifest_path = root / "source_bundle_manifest.json"
    if not manifest_path.exists():
        return {"manifest_present": False}
    manifest = json.loads(manifest_path.read_text())
    protocol_hash = sha256(protocol_path)
    if protocol_hash != manifest.get("protocol_sha256_before_outcome_fit"):
        raise ValueError("Current protocol differs from the frozen prefit source-bundle protocol")
    entries = manifest.get("files", [])
    if not entries:
        raise ValueError("Frozen source bundle has no file checksums")
    seen = set()
    for entry in entries:
        path = (root / entry["path"]).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Manifest file is outside its source bundle: " + entry["path"])
        if path in seen:
            raise ValueError("Duplicate frozen source-bundle file: " + entry["path"])
        seen.add(path)
        if (not path.is_file() or path.stat().st_size != entry["bytes"] or
                sha256(path) != entry["sha256"]):
            raise ValueError(entry["path"] + " differs from its frozen size or checksum")
    return {"manifest_present": True, "manifest_sha256": sha256(manifest_path),
            "files_verified_count": len(entries),
            "protocol_sha256_before_outcome_fit": protocol_hash,
            "protocol_frozen_at_utc": manifest.get("protocol_frozen_at_utc")}


def verify_centroid_distances(pairs, groups, coordinates, tolerance_km=1e-6):
    centers = {key: spherical_mean([coordinates[port] for port in ports])
               for key, ports in groups.items()}
    errors = []
    for row in pairs.itertuples(index=False):
        a = centers[str(row.island_a)]["unit_sphere_xyz"]
        b = centers[str(row.island_b)]["unit_sphere_xyz"]
        errors.append(abs(arc_distance(a, b)-float(row.port_distance_km)))
    maximum = max(errors, default=float("inf"))
    if not np.isfinite(maximum) or maximum > tolerance_km:
        raise ValueError(f"Original WPI groups do not reproduce archived pair distances: {maximum} km")
    return {"source_pairs_recomputed": len(errors),
            "max_original_pair_distance_error_km": maximum,
            "original_port_centroids_recomputed_from_raw_WPI": True}


def validate_source_files(ports_path, recovery_path, pairs_path, islands_path):
    recovery = json.loads(Path(recovery_path).read_text())
    if recovery.get("status") != RECOVERY_STATUS:
        raise ValueError("Original51-island recovery did not pass its source gate")
    if sha256(ports_path) != recovery.get("source_WPI_port_checksum"):
        raise ValueError("Original WPI source checksum differs from the verified recovery")
    if sha256(pairs_path) != recovery.get("source_pair_sha256"):
        raise ValueError("Original biological/geographic pair checksum differs from recovery")
    q, ids, _ = source(Path(pairs_path), Path(islands_path))
    groups = {str(k): [str(p) for p in v] for k, v in recovery["source_51_port_groups"].items()}
    all_ports = [p for ports in groups.values() for p in ports]
    if (set(groups) != set(ids) or len(groups) != 51 or len(all_ports) != 194 or
            len(set(all_ports)) != 194 or recovery.get("total_51_original_port_assignments") != 194):
        raise ValueError("The original51 groups and194 unique port assignments must be preserved")
    geojson = json.loads(Path(ports_path).read_text())
    if geojson.get("type") != "FeatureCollection":
        raise ValueError("Original WPI coordinates are not a FeatureCollection")
    coordinates, properties = {}, {}
    for feature in geojson["features"]:
        key = str(int(feature["properties"]["OBJECTID"]))
        if key in coordinates or feature["geometry"]["type"] != "Point":
            raise ValueError("Duplicated original WPI ID or nonpoint geometry")
        xy = feature["geometry"]["coordinates"]
        if len(xy) != 2 or not np.isfinite(xy).all():
            raise ValueError("Original WPI port has invalid coordinates")
        coordinates[key] = tuple(map(float, xy))
        properties[key] = feature["properties"]
    if set(all_ports)-set(coordinates):
        raise ValueError("Original recovered port IDs are missing from WPI coordinates")
    roster = pd.read_csv(islands_path, dtype={"landmass_id": str})
    counts = {str(row.landmass_id): int(row.ports) for row in roster.itertuples()
              if str(row.landmass_id) in groups}
    if any(counts[key] != len(groups[key]) for key in ids):
        raise ValueError("Original island port counts differ from recovered groups")
    receipt = verify_centroid_distances(q, groups, coordinates)
    receipt.update({"original_islands": len(ids), "original_unique_ports": len(all_ports),
                    "WPI_sha256": sha256(ports_path), "recovery_sha256": sha256(recovery_path),
                    "pairs_sha256": sha256(pairs_path), "islands_sha256": sha256(islands_path)})
    return q, ids, groups, coordinates, properties, counts, receipt


def eligible_port_groups(groups, attached, *, marine_ports_only=False):
    missing = [p for ports in groups.values() for p in ports if p not in attached]
    if missing and not marine_ports_only:
        raise ValueError(f"{len(missing)} original ports lack verified marine access; explicit marine filter required")
    kept = {key: [p for p in ports if p in attached] for key, ports in groups.items()}
    empty = [key for key, ports in kept.items() if not ports]
    if empty:
        raise ValueError("Marine-port filter would remove an original island: " + ",".join(empty))
    return kept


def load_wind(path, target_shape):
    with np.load(path) as data:
        values = np.asarray(data["mean_monthly_wind_vector_speed_10m"], dtype=float)
        lat = np.asarray(data["lat_centers"], dtype=float)
        lon = np.asarray(data["lon_centers"], dtype=float)
    if (values.shape != (170, 360) or lat.shape != (170,) or lon.shape != (360,) or
            not np.allclose(lat, np.arange(84.5, -85., -1.), atol=1e-12, rtol=0.) or
            not np.allclose(lon, np.arange(-179.5, 180., 1.), atol=1e-12, rtol=0.)):
        raise ValueError("NOAA wind values are not aligned with their exact1-degree source grid")
    if tuple(target_shape) != (340, 720):
        raise ValueError("Frozen wind transfer requires the exact0.5-degree output grid")
    return np.repeat(np.repeat(values, 2, axis=0), 2, axis=1)


def finite_volume_edges(topology, conductivity):
    rows = np.asarray(topology["node_block_row"], dtype=int)
    cols = np.asarray(topology["node_block_col"], dtype=int)
    latitude = np.asarray(topology["block_latitude"], dtype=float)[rows]
    longitude = np.asarray(topology["block_longitude"], dtype=float)[cols]
    i = np.asarray(topology["edge_i"], dtype=int)
    j = np.asarray(topology["edge_j"], dtype=int)
    crossing = np.asarray(topology["boundary_crossing_count"], dtype=float)
    if (crossing.shape != i.shape or j.shape != i.shape or
            not np.isfinite(crossing).all() or np.any(crossing <= 0) or
            np.any(crossing != np.rint(crossing))):
        raise ValueError("Every edge requires a positive integer source water-face contact count")
    dr = np.abs(rows[i]-rows[j])
    raw_dc = np.abs(cols[i]-cols[j])
    dc = np.minimum(raw_dc, len(topology["block_longitude"])-raw_dc)
    east_west = (dr == 0) & (dc == 1)
    north_south = (dr == 1) & (dc == 0)
    if not np.all(east_west | north_south):
        raise ValueError("Source edges must join orthogonal adjacent blocks, including longitude wrap")
    fine_step = float(abs(np.diff(topology["fine_longitude"])[0]))
    face = crossing*circuit.EARTH_RADIUS_KM*np.deg2rad(fine_step)
    shared_latitude = (latitude[i]+latitude[j])/2
    face[north_south] *= np.cos(np.deg2rad(shared_latitude[north_south]))
    distance, per_km = circuit.edge_conductance(longitude, latitude, i, j, conductivity)
    conductance = per_km*face
    if not np.isfinite(conductance).all() or np.any(conductance <= 0):
        raise ValueError("Source finite-volume edge conductance is nonpositive or nonfinite")
    return {"distance_km": distance, "face_width_km": face, "conductance": conductance}


def make_base_features(q, ids, original_port_counts, target):
    features = make_features(q, None, response=target)
    counts = np.array([original_port_counts[key] for key in ids], dtype=float)
    if not np.isfinite(counts).all() or np.any(counts <= 0):
        raise ValueError("Original eligible port counts must be finite and positive")
    logs = np.log1p(counts)
    left, right = logs[q.i.to_numpy()], logs[q.j.to_numpy()]
    features["sum_log1p_original_port_count"] = left+right
    features["abs_difference_log1p_original_port_count"] = np.abs(left-right)
    return features


def compare_models(q, ids, original_port_counts, matrices):
    if set(matrices) != {"uniform", "ais", "wind"}:
        raise ValueError("All three frozen marine circuit matrices are required")
    features = {key: circuit.resistance_feature(value, q.i.to_numpy(), q.j.to_numpy())
                for key, value in matrices.items()}
    result = {key: {"targets": {}} for key in
              ("base_vs_uniform", "uniform_vs_ais", "uniform_vs_wind", "uniform_vs_ais_and_wind")}
    for target in ("host_jaccard", "new_butterfly_jaccard", "earlier_butterfly_jaccard"):
        base = make_base_features(q, ids, original_port_counts, target)
        uniform = dict(base, log1p_normalized_uniform_circuit=features["uniform"])
        ais = dict(uniform, log1p_normalized_ais_circuit=features["ais"])
        wind = dict(uniform, log1p_normalized_wind_circuit=features["wind"])
        combined = dict(ais, log1p_normalized_wind_circuit=features["wind"])
        for label, reference, augmented in [
                ("base_vs_uniform", base, uniform),
                ("uniform_vs_ais", uniform, ais),
                ("uniform_vs_wind", uniform, wind),
                ("uniform_vs_ais_and_wind", uniform, combined)]:
            model = holdout(q, reference, augmented, target, ids)
            model["earlier_target_is_temporal_comparator_not_zero_association_null"] = target == "earlier_butterfly_jaccard"
            result[label]["targets"][target] = model
    return result


def export_networks(directory, topology, edge_fields, weights, ids, component_nodes):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    component_nodes = np.asarray(component_nodes, dtype=int)
    remap = np.full(weights.shape[1], -1, dtype=int)
    remap[component_nodes] = np.arange(1, len(component_nodes)+1)
    i, j = topology["edge_i"], topology["edge_j"]
    retain = (remap[i] > 0) & (remap[j] > 0)
    hashes = {}
    for key, value in edge_fields.items():
        path = directory / f"{key}_conductances.txt"
        edges = np.column_stack([remap[i[retain]], remap[j[retain]], value["conductance"][retain]])
        np.savetxt(path, edges, fmt=["%d", "%d", "%.17g"])
        hashes[path.name] = sha256(path)
    rows = []
    for index, island in enumerate(ids):
        row = weights.getrow(index)
        for node, weight in zip(row.indices, row.data):
            if remap[node] < 1:
                raise ValueError("A terminal port is outside the exported common marine component")
            rows.append({"landmass_id": island, "network_node_id": int(remap[node]),
                         "source_graph_node_id": int(node), "weight": float(weight)})
    wpath = directory / "island_node_weights.csv"
    pd.DataFrame(rows).to_csv(wpath, index=False, float_format="%.17g")
    hashes[wpath.name] = sha256(wpath)
    np.savez_compressed(directory / "network_node_mapping.npz", source_graph_nodes=component_nodes)
    return {"nodes": len(component_nodes), "undirected_edges": int(retain.sum()),
            "one_edge_orientation_exported": True, "node_numbering": "compact1-based",
            "files_sha256": hashes}


def prepare_inputs(input_root, coast_dir, out_dir, network_dir, *, marine_ports_only):
    root, coast, out = map(Path, (input_root, coast_dir, out_dir))
    bundle_receipt = validate_source_bundle(root)
    out.mkdir(parents=True, exist_ok=True)
    ports = root / "wpi/ports_2019_source.geojson"
    recovery = root / "coordinates/island51_original_WPI_ports_coords_unique_source.json"
    pairs = SOURCE_DIRECTORY / "island_host_network_pairs.csv"
    islands = SOURCE_DIRECTORY / "island_host_network_islands.csv"
    q, ids, groups, coordinates, properties, counts, source_receipt = validate_source_files(
        ports, recovery, pairs, islands)
    graph_path = coast / "island51_coastal_ocean_graph.npz"
    graph_receipt_path = coast / "island51_coastal_ocean_graph_receipt.json"
    graph_receipt = json.loads(graph_receipt_path.read_text())
    if (sha256(graph_path) != graph_receipt["graph_npz_sha256"] or
            graph_receipt["source"]["GSHHG_original_archive_sha256"] != GSHHG_SOURCE_SHA):
        raise ValueError("Source coastline graph checksum or GSHHG archive identity changed")
    if (graph_receipt["fine_grid_degrees"], graph_receipt["block_grid_degrees"]) != (.0625, .5):
        raise ValueError("Frozen primary requires0.0625-degree coast subcells and0.5-degree blocks")
    with np.load(graph_path) as data:
        topology = {key: np.asarray(data[key]) for key in data.files}
    fine_step = float(graph_receipt["fine_grid_degrees"])
    components = topology["node_component"]
    total_water = np.bincount(components, weights=topology["node_fine_water_cells"])
    marine_component = int(np.argmax(total_water))
    marine_nodes = np.flatnonzero(components == marine_component)
    # Restrict candidate attachments by source coastline alone. Inland-water
    # components are not marine ports and no bridge is created for them.
    fine_nodes = topology["fine_cell_to_node"].copy()
    wet = fine_nodes >= 0
    nonmarine = np.zeros(fine_nodes.shape, dtype=bool)
    nonmarine[wet] = components[fine_nodes[wet]] != marine_component
    fine_nodes[nonmarine] = -1
    attached, attachment_rows = {}, []
    for island in ids:
        for port in groups[island]:
            lon, lat = coordinates[port]
            row = {"landmass_id": island, "port_id": port, "longitude": lon, "latitude": lat,
                   "source_harbor_type": properties[port].get("HARBORTYPE", ""),
                   "source_port_name": properties[port].get("PORT_NAME", "")}
            try:
                snapped = circuit.nearest_sea_cell(lon, lat, fine_nodes, fine_step, distance_limit_km=25.)
                attached[port] = snapped["node"]
                row.update(snapped)
                row.update(included=True, status="VERIFIED_MARINE_ACCESS_WITHIN_25_KM")
            except ValueError as exc:
                row.update(included=False, status="NO_VERIFIED_MARINE_ACCESS_WITHIN_25_KM", reason=str(exc))
            attachment_rows.append(row)
    ledger = pd.DataFrame(attachment_rows)
    ledger.to_csv(out / "original194_port_attachment_audit.csv", index=False)
    kept = eligible_port_groups(groups, attached, marine_ports_only=marine_ports_only)
    if marine_ports_only and len(attached) != 182:
        raise ValueError("Frozen primary marine eligibility changed from182 of194 original ports")
    weights = circuit.port_weights(ids, kept, attached, len(components))
    aliases = []
    for i in range(len(ids)):
        for j in range(i+1, len(ids)):
            if (weights.getrow(i)-weights.getrow(j)).nnz == 0:
                aliases.append([ids[i], ids[j]])
    if aliases:
        raise ValueError("Original islands have indistinguishable marine source distributions: " + str(aliases))
    if len(np.unique(components[weights.nonzero()[1]])) != 1:
        raise ValueError("Qualifying original ports occupy disconnected marine components")
    ais_path = root / "ais/commercial_AIS_mean_0p5deg_2015_2021.npz"
    ais_receipt_path = root / "ais/commercial_AIS_mean_0p5deg_2015_2021_receipt.json"
    ais_receipt = json.loads(ais_receipt_path.read_text())
    if ais_receipt.get("source_original_SHA256") != AIS_SOURCE_SHA:
        raise ValueError("Original World Bank commercial AIS archive identity changed")
    ais, valid = circuit.load_ais(ais_path)
    if ais.shape != (340, 720):
        raise ValueError("Frozen primary AIS grid is not0.5-degree")
    wind_path = root / "wind/NOAA_NCEP_R1_10m_wind2015to2021_1deg.npz"
    wind_receipt_path = root / "wind/NOAA_NCEP_R1_10m_wind2015to2021_1deg_receipt.json"
    wind_receipt = json.loads(wind_receipt_path.read_text())
    if wind_receipt.get("observed_months") != 74 or wind_receipt.get("source_wind_height_m") != 10:
        raise ValueError("Frozen NOAA source period or wind height changed")
    wind = load_wind(wind_path, ais.shape)
    node_row, node_col = topology["node_block_row"], topology["node_block_col"]
    sea = np.zeros(ais.shape, dtype=bool)
    sea[node_row, node_col] = True
    if not np.all(valid[sea]):
        raise ValueError("Original AIS has missing coastal-graph cells; do not silently zero-fill or alter common support")
    _, fields = circuit.conductance_fields(ais, valid, wind, sea)
    edge_fields = {key: finite_volume_edges(topology, value[node_row, node_col])
                   for key, value in fields.items()}
    exported = export_networks(network_dir, topology, edge_fields, weights, ids, marine_nodes)
    included = ledger.loc[ledger.included]
    excluded = ledger.loc[~ledger.included]
    by_island = {key: {"original_ports": len(groups[key]), "marine_ports": len(kept[key]),
                       "unique_terminal_nodes": int(weights.getrow(index).nnz)}
                 for index, key in enumerate(ids)}
    receipt = {
        "schema": "chocho_island51_source_observed_ocean_circuit_preparation_v1",
        "status": "SOURCE_GEOMETRY_AND_CONDUCTANCE_READY_NO_BIOLOGICAL_MODELS_FITTED",
        "protocol_sha256": sha256(PROTOCOL), "original_source": source_receipt,
        "frozen_source_bundle": bundle_receipt,
        "graph": graph_receipt, "marine_component": marine_component,
        "original_ports_audited": len(ledger), "marine_ports_included": len(included),
        "ports_excluded": len(excluded), "islands_retained": len(ids),
        "explicit_marine_port_filter": marine_ports_only,
        "excluded_harbor_type_counts": dict(Counter(excluded.source_harbor_type)),
        "excluded_original_port_ids": sorted(excluded.port_id.tolist(), key=int),
        "maximum_snap_distance_km": float(included.displacement_km.max()),
        "median_snap_distance_km": float(included.displacement_km.median()),
        "by_island": by_island, "island_distributions_distinct": True,
        "original_port_counts_used_as_model_controls": True,
        "all_models_use_identical_marine_domain_and_port_weights": True,
        "AIS_source_receipt": ais_receipt, "wind_source_receipt": wind_receipt,
        "wind_regridding": "repeat each observed1-degree parent over four0.5-degree children; no new observed resolution",
        "inputs_sha256": {"AIS_grid": sha256(ais_path), "AIS_receipt": sha256(ais_receipt_path),
                          "wind_grid": sha256(wind_path), "wind_receipt": sha256(wind_receipt_path),
                          "coast_graph": sha256(graph_path), "coast_receipt": sha256(graph_receipt_path)},
        "exported_network": exported,
        "official_crosscheck_pairs_selected_from_geometry": [["9", "15"], ["715", "1577"]],
    }
    (out / "source_preparation.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print("OBSERVED_OCEAN_PREPARED", json.dumps({"islands": len(ids), "ports": len(included),
          "original_port_denominator": len(ledger), "marine_nodes": len(marine_nodes),
          "network_edges": exported["undirected_edges"]}), flush=True)
    return q, ids, counts, topology, edge_fields, weights, receipt


def run(input_root, coast_dir, out_dir, network_dir, *, marine_ports_only=False, prepare_only=False):
    started = time.monotonic()
    q, ids, counts, topology, edge_fields, weights, preparation = prepare_inputs(
        input_root, coast_dir, out_dir, network_dir, marine_ports_only=marine_ports_only)
    if prepare_only:
        return preparation
    out = Path(out_dir)
    i, j = topology["edge_i"], topology["edge_j"]
    n = weights.shape[1]
    matrices, numerical = {}, {}
    for label in ("uniform", "ais", "wind"):
        began = time.monotonic()
        g = edge_fields[label]["conductance"]
        adjacency = sparse.csr_matrix((np.r_[g, g], (np.r_[i, j], np.r_[j, i])), shape=(n, n))
        laplacian = sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency
        result, diagnostics = circuit.distributed_effective_resistance(laplacian, weights)
        diagnostics["elapsed_seconds"] = time.monotonic()-began
        matrices[label], numerical[label] = result, diagnostics
        np.savez_compressed(out / f"{label}_distributed_port_resistance.npz",
                            island_ids=np.array(ids), resistance=result)
        print("OBSERVED_OCEAN_CIRCUIT_SOLVED", label, json.dumps(diagnostics), flush=True)
    comparisons = compare_models(q, ids, counts, matrices)
    pair_output = q.drop(columns=["i", "j"]).copy()
    for label, values in matrices.items():
        pair_output[f"{label}_distributed_port_resistance"] = values[q.i.to_numpy(), q.j.to_numpy()]
        pair_output[f"log1p_normalized_{label}_circuit"] = circuit.resistance_feature(values, q.i.to_numpy(), q.j.to_numpy())
    pair_output.to_csv(out / "island51_observed_circuit_pairs.csv", index=False, float_format="%.17g")
    fold_rows, summary_rows = [], []
    for label, comparison in comparisons.items():
        for target, model in comparison["targets"].items():
            summary_rows.append({"comparison": label, "target": target,
                                 "MSE_reference": model["MSE_baseline"],
                                 "MSE_augmented": model["MSE_augmented"],
                                 "relative_MSE_reduction_percent": model["relative_MSE_reduction_percent"],
                                 "MAE_reference": model["MAE_baseline"], "MAE_augmented": model["MAE_augmented"],
                                 "islands_improved": model["islands_with_better_augmented_prediction"]})
            for fold in model["per_island_LOIO"]:
                fold_rows.append(dict(comparison=label, target=target, **fold))
    pd.DataFrame(summary_rows).to_csv(out / "model_comparison_summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(fold_rows).to_csv(out / "all_island_holdout_losses.csv", index=False, float_format="%.17g")
    result = {
        "schema": "chocho_island51_observed_ocean_distributed_port_holdout_v1",
        "status": "THREE_OBSERVED_MARINE_CIRCUITS_AND_FIXED_ISLAND_HOLDOUT_COMPLETE",
        "source_preparation": preparation, "numerical_diagnostics": numerical,
        "comparisons": comparisons, "elapsed_seconds": time.monotonic()-started,
        "official_circuitscape_crosscheck_performed_by_this_script": False,
        "official_circuitscape_validation_receipt": "Separate pinned Julia/Circuitscape validation receipt, if executed",
        "interpretation_limits": [
            "AIS is modern position-report density, not bilateral voyages, historical cargo or plant carriage.",
            "Wind is a scalar mean of monthly10m vector speeds, not directional or flight-altitude dispersal.",
            "Marine-access filtering retains182/194 audited original ports and all51 original biological islands.",
            "Observed first records do not establish true butterfly colonization; pooled plants can postdate butterfly records.",
            "All1275 pairs share endpoints and51 fold losses are dependent; no pair-independent p-values or confidence intervals.",
            "Earlier butterfly composition is a temporal comparison, not a null expected to lack geographic association.",
            "Coastline block aggregation preserves connectivity, not exact fine-grid electrical resistance.",
        ],
    }
    (out / "observed_ocean_circuit_holdout.json").write_text(json.dumps(result, indent=2)+"\n")
    print("OBSERVED_OCEAN_HOLDOUT_COMPLETE", json.dumps(summary_rows), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--coast-dir", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--network-dir", type=Path, required=True)
    parser.add_argument("--marine-ports-only", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    run(args.input_root, args.coast_dir or args.input_root / "coast_0p0625", args.out_dir,
        args.network_dir, marine_ports_only=args.marine_ports_only, prepare_only=args.prepare_only)


if __name__ == "__main__":
    main()

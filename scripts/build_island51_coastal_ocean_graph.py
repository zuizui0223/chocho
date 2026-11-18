"""Build a coastline-resolved ocean graph for the fixed island51 analysis.

Each node is one four-connected sea component within an aggregation block.
Edges exist only across observed adjacent fine-grid sea centers. This preserves
the fine-grid topology, including resolved straits and separated water bodies
within a single block; it is not a resistance-equivalent circuit reduction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import shapely
from shapely.strtree import STRtree
from scipy.ndimage import label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


SOURCE_GSHHG_SHA = "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_source_polygons(original_zip):
    checksum = sha256(original_zip)
    if checksum != SOURCE_GSHHG_SHA:
        raise ValueError("Pinned GSHHG original v2.3.7 coastline checksum mismatch")
    import geopandas as gpd

    member = "GSHHS_shp/i/GSHHS_i_L1.shp"
    geodata = gpd.read_file(f"zip://{Path(original_zip).resolve()}!{member}")
    if not geodata.crs:
        raise ValueError("Source coastline has no coordinate reference system")
    geodata = geodata.to_crs(4326)
    polygons = geodata.geometry[geodata.geometry.notna()].to_numpy()
    polygons = polygons[~shapely.is_empty(polygons)]
    if not len(polygons):
        raise ValueError("Source coastline has no nonempty polygons")
    invalid = ~shapely.is_valid(polygons)
    if np.any(invalid):
        polygons[invalid] = shapely.make_valid(polygons[invalid])
    return polygons, {
        "GSHHG_original_archive_sha256": checksum,
        "GSHHG_member": member,
        "source_polygon_features": int(len(polygons)),
        "invalid_source_features_repaired_with_make_valid": int(invalid.sum()),
    }


def rasterize_land(polygons, latitudes, longitudes, chunk_rows=32, progress=False):
    if chunk_rows < 1:
        raise ValueError("Positive chunk row count is required")
    polygons = np.asarray(polygons, dtype=object)
    tree = STRtree(polygons)
    shapely.prepare(polygons)
    land = np.zeros((len(latitudes), len(longitudes)), dtype=bool)
    for start in range(0, len(latitudes), chunk_rows):
        end = min(start + chunk_rows, len(latitudes))
        longitude, latitude = np.meshgrid(longitudes, latitudes[start:end])
        points = shapely.points(longitude.ravel(), latitude.ravel())
        candidates = tree.query(points)
        hit = shapely.intersects(polygons[candidates[1]], points[candidates[0]])
        flat_land = np.zeros(len(points), dtype=bool)
        flat_land[candidates[0, hit]] = True
        land[start:end] = flat_land.reshape(end - start, len(longitudes))
        if progress and (end % 256 == 0 or end == len(latitudes)):
            print("COAST_MASK_PROGRESS", json.dumps({"completed_rows": end,
                                                     "total_rows": len(latitudes)}), flush=True)
    return land


def build_coastal_graph(sea, block_size_cells=4, wrap_longitude=True):
    sea = np.asarray(sea, dtype=bool)
    size = int(block_size_cells)
    if sea.ndim != 2 or size < 1 or size != block_size_cells:
        raise ValueError("A two-dimensional sea mask and positive integer block size are required")
    height, width = sea.shape
    if not height or not width or height % size or width % size:
        raise ValueError("Fine-grid dimensions must be positive multiples of block size")
    ny, nx = height // size, width // size
    blocks = sea.reshape(ny, size, nx, size).transpose(0, 2, 1, 3)
    wet_counts = blocks.sum(axis=(2, 3))
    full_water = wet_counts == size * size
    mixed = (wet_counts > 0) & ~full_water
    local_labels = np.zeros(blocks.shape, dtype=np.int32)
    local_labels[full_water] = 1
    local_counts = full_water.astype(np.int32)
    cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    for by, bx in np.argwhere(mixed):
        local_labels[by, bx], local_counts[by, bx] = label(blocks[by, bx], structure=cross)

    counts = local_counts.ravel()
    offsets = np.cumsum(counts, dtype=np.int64) - counts
    total_nodes = int(counts.sum())
    if total_nodes > np.iinfo(np.int32).max:
        raise ValueError("Graph exceeds supported int32 node index range")
    block_nodes = np.where(local_labels > 0,
                           local_labels + offsets.reshape(ny, nx, 1, 1) - 1,
                           -1).astype(np.int32)
    fine_nodes = block_nodes.transpose(0, 2, 1, 3).reshape(height, width)
    block_ids = np.arange(ny * nx, dtype=np.int32)
    node_rows = np.repeat(block_ids // nx, counts)
    node_cols = np.repeat(block_ids % nx, counts)

    # The slices include only cross-block fine contacts. Internal contacts
    # are already represented by the connected components within each block.
    contacts = [
        (fine_nodes[size - 1:-1:size, :], fine_nodes[size::size, :]),
        (fine_nodes[:, size - 1:-1:size], fine_nodes[:, size::size]),
    ]
    if wrap_longitude:
        contacts.append((fine_nodes[:, -1], fine_nodes[:, 0]))
    codes = []
    for left, right in contacts:
        valid = (left >= 0) & (right >= 0) & (left != right)
        lo = np.minimum(left[valid], right[valid]).astype(np.int64)
        hi = np.maximum(left[valid], right[valid]).astype(np.int64)
        codes.append(lo * total_nodes + hi)
    all_codes = np.concatenate(codes)
    if total_nodes:
        unique_codes, crossings = np.unique(all_codes, return_counts=True)
        edge_i = (unique_codes // total_nodes).astype(np.int32)
        edge_j = (unique_codes % total_nodes).astype(np.int32)
        adjacency = coo_matrix((np.ones(len(edge_i) * 2, dtype=np.uint8),
                                (np.concatenate((edge_i, edge_j)),
                                 np.concatenate((edge_j, edge_i)))),
                               shape=(total_nodes, total_nodes)).tocsr()
        n_components, components = connected_components(adjacency, directed=False)
        node_water_cells = np.bincount(fine_nodes[sea], minlength=total_nodes).astype(np.int32)
        component_nodes = np.bincount(components)
        component_water = np.bincount(components, weights=node_water_cells).astype(np.int64)
    else:
        edge_i = edge_j = crossings = components = node_water_cells = np.array([], dtype=np.int32)
        component_nodes = component_water = np.array([], dtype=np.int64)
        n_components = 0
    arrays = {
        "fine_cell_to_node": fine_nodes,
        "node_block_row": node_rows,
        "node_block_col": node_cols,
        "node_fine_water_cells": node_water_cells,
        "node_component": components.astype(np.int32),
        "edge_i": edge_i,
        "edge_j": edge_j,
        "boundary_crossing_count": crossings.astype(np.int32),
    }
    audit = {
        "fine_grid_shape": [height, width],
        "block_grid_shape": [ny, nx],
        "block_size_cells": size,
        "fine_water_cells": int(sea.sum()),
        "fine_land_cells": int((~sea).sum()),
        "wet_blocks": int(np.count_nonzero(wet_counts)),
        "mixed_land_water_blocks": int(mixed.sum()),
        "blocks_with_separated_water_components": int(np.sum(local_counts > 1)),
        "graph_nodes": total_nodes,
        "undirected_edges": int(len(edge_i)),
        "boundary_fine_contacts": int(crossings.sum()),
        "number_of_components": int(n_components),
        "component_node_counts": component_nodes.tolist(),
        "component_fine_water_cell_counts": component_water.tolist(),
        "fine_neighbor_definition": 4,
        "longitude_wrap": bool(wrap_longitude),
        "aggregation": "separate four-connected sea components inside each block",
        "conductance_scope": "topology is preserved; aggregation is not an exact circuit reduction",
    }
    return arrays, audit


def run(original_zip, out_dir, *, fine_degrees=0.125, block_degrees=0.5, chunk_rows=32):
    started = time.monotonic()
    if fine_degrees <= 0 or block_degrees <= 0 or not np.all(np.isfinite([fine_degrees, block_degrees])):
        raise ValueError("Positive finite fine and block resolutions are required")
    ratios = np.array([170 / fine_degrees, 360 / fine_degrees, block_degrees / fine_degrees])
    if not np.all(np.isclose(ratios, np.rint(ratios), atol=1e-8, rtol=0)):
        raise ValueError("Fine resolution must tile the globe and aggregation blocks exactly")
    height, width, size = np.rint(ratios).astype(int)
    latitude = 85 - (np.arange(height) + 0.5) * fine_degrees
    longitude = -180 + (np.arange(width) + 0.5) * fine_degrees
    polygons, source = load_source_polygons(original_zip)
    land = rasterize_land(polygons, latitude, longitude, chunk_rows=chunk_rows, progress=True)
    arrays, graph = build_coastal_graph(~land, block_size_cells=int(size), wrap_longitude=True)
    arrays.update({
        "fine_latitude": latitude,
        "fine_longitude": longitude,
        "block_latitude": 85 - (np.arange(height // size) + 0.5) * block_degrees,
        "block_longitude": -180 + (np.arange(width // size) + 0.5) * block_degrees,
    })
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    destination = out_dir / "island51_coastal_ocean_graph.npz"
    np.savez_compressed(destination, **arrays)
    receipt = {
        "schema": "chocho_coastline_connected_block_ocean_graph_v1",
        "status": "SOURCE_VERIFIED_COASTLINE_GRAPH_BUILT_ANCHORS_NOT_EVALUATED",
        "source": source,
        "fine_grid_degrees": fine_degrees,
        "block_grid_degrees": block_degrees,
        "extent_degrees": {"west": -180, "east": 180, "north": 85, "south": -85},
        "mask_method": "intermediate GSHHG L1 polygons intersect fine-grid centers; polygon boundaries are land",
        "graph": graph,
        "graph_npz_sha256": sha256(destination),
        "graph_npz_bytes": destination.stat().st_size,
        "elapsed_seconds": time.monotonic() - started,
        "limits": [
            "Only channels resolved by the declared fine grid are retained; no inferred bridges are added.",
            "Fine-grid center adjacency can miss coastline features narrower than the fine resolution.",
            "Topology is preserved on aggregation, but effective resistance is not an exact fine-circuit reduction.",
            "Island terminals and biological comparisons have not been evaluated by this topology builder.",
        ],
    }
    (out_dir / "island51_coastal_ocean_graph_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gshhg-zip", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--fine-degrees", type=float, default=0.125)
    parser.add_argument("--block-degrees", type=float, default=0.5)
    parser.add_argument("--chunk-rows", type=int, default=32)
    args = parser.parse_args()
    receipt = run(args.gshhg_zip, args.out_dir, fine_degrees=args.fine_degrees,
                  block_degrees=args.block_degrees, chunk_rows=args.chunk_rows)
    graph = receipt["graph"]
    print("COASTAL_GRAPH_BUILT", json.dumps({
        "status": receipt["status"], "graph_nodes": graph["graph_nodes"],
        "undirected_edges": graph["undirected_edges"],
        "number_of_components": graph["number_of_components"],
        "largest_components_nodes": sorted(graph["component_node_counts"], reverse=True)[:8],
        "elapsed_seconds": receipt["elapsed_seconds"],
    }), flush=True)


if __name__ == "__main__":
    main()

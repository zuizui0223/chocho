#!/usr/bin/env python3
"""Exact source-equivalent UNIFORM GSHHG finite-volume ocean circuit for WCVP126.

The requested uniform-sea-CIRCUIT control is DIFFERENT from shortest-water
distance. Use the original 0.0625deg GSHHG water-component-aware 0.5deg
finite-volume conductance definitions from the independent island51 original.
Attach original 104 botanical regions by independently coastline-verified
TDWG3 water nodes, NOT original WPI ports and never botanical outcomes.

104×104 effective-resistance matrix solves same original 189550-node
connected marine component Laplacian with SciPy sparse LU, original verified
face-width / geographic edge-length conductive weights. Expensive but exact
at the source grid's coarse aggregation, NOT exact continuum coastline physics.

No ocean current direction, vessel transport, viable seed, historical date or
butterfly establishment is measured by this undirected electrical circuit.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time
from pathlib import Path

import numpy as np
from scipy import sparse
import island51_ocean_circuit as circuit
from run_island51_observed_ocean_circuits import finite_volume_edges

SOURCE_GSHHG_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
EXPECTED_NODES=193516
EXPECTED_EDGES=372096
EXPECTED_GIANT_NODES=189550
EXPECTED_GIANT_EDGES=370996

def SHA(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for v in iter(lambda:f.read(1024*1024),b""):h.update(v)
    return h.hexdigest()

def rows(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-GSHHG-graph",type=Path,required=True)
    p.add_argument("--source-coast-audit",type=Path,required=True)
    p.add_argument("--source-pairs",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    started=time.monotonic()
    if SHA(a.original_GSHHG_graph)!=SOURCE_GSHHG_SHA:
        raise RuntimeError("Original 0.0625deg source GSHHG graph hash changed")
    coast=json.loads(a.source_coast_audit.read_text())
    if (coast["source_original_126_regions"]!=126 or
        coast["source_126_regions_with_at_least_one_verified_coastal_water_anchor"]!=104 or
        coast["source_valid_original_island_pair_geometric_shortest_sea"]!=5356):
        raise RuntimeError("Original 126 island GSHHG source coverage changed")
    original=rows(a.source_pairs)
    if len(original)!=7875:
        raise RuntimeError("Original 7875 accepted WCVP source pair denominator changed")
    anchored=coast["coastal_anchors_by_TDWG_region"]
    island_ids=sorted(k for k,v in anchored.items() if v)
    if len(island_ids)!=104:raise RuntimeError("Original source 104 island terminal IDs not restored")
    topo_keys=("node_component","node_block_row","node_block_col",
      "block_latitude","block_longitude","edge_i","edge_j",
      "boundary_crossing_count","fine_longitude")
    with np.load(a.original_GSHHG_graph) as data:
        topo={k:np.asarray(data[k]) for k in topo_keys}
    node_components=topo["node_component"]
    largest=int(np.argmax(np.bincount(node_components)))
    giant=np.flatnonzero(node_components==largest)
    u=np.asarray(topo["edge_i"],int)
    v=np.asarray(topo["edge_j"],int)
    if (len(node_components)!=EXPECTED_NODES or len(u)!=EXPECTED_EDGES or
        len(giant)!=EXPECTED_GIANT_NODES):
        raise RuntimeError("Original frozen GSHHG graph source nodes/edges changed")
    ret=(node_components[u]==largest)&(node_components[v]==largest)
    if int(ret.sum())!=EXPECTED_GIANT_EDGES:
        raise RuntimeError("Original frozen largest marine component edge population changed")
    # Physically original finite volume boundary face-width / center-distance
    # conductance with a UNIFORM (1) ocean field. Do not substitute shortest
    # paths or spurious 1/distance-only edge weights.
    result=finite_volume_edges(topo,np.ones(len(node_components),dtype=float))
    cond=result["conductance"]
    if np.any(~np.isfinite(cond)) or np.any(cond<=0):
        raise RuntimeError("Source original marine graph source conductances not positive")
    adjacency=sparse.csr_matrix((np.r_[cond,cond],
        (np.r_[u,v],np.r_[v,u])),
        shape=(len(node_components),)*2)
    lap=sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency
    rr=[];cc=[];ww=[]
    for i,code in enumerate(island_ids):
        pts=anchored[code]
        if not (1<=len(pts)<=3):raise RuntimeError("Unexpected TDWG coastline terminal counts")
        for anchor in pts:
            node=int(anchor["node"])
            if node<0 or node>=len(node_components) or node_components[node]!=largest:
                raise RuntimeError(f"Original region coastal terminal outside largest source marine component {code}")
            if not anchor.get("official_polygon_coast_verified"):
                raise RuntimeError("Original TDWG polygon coastal identity was not verified")
            rr.append(i);cc.append(node);ww.append(1./len(pts))
    weights=sparse.csr_matrix((ww,(rr,cc)),shape=(len(island_ids),len(node_components)))
    if not np.allclose(np.asarray(weights.sum(axis=1)).ravel(),1.,atol=1e-12):
        raise RuntimeError("Original region coastal weights don't sum to one")
    uniq=len(np.unique(weights.nonzero()[1]))
    print("WCVP126_ORIGINAL_UNIFORM_CIRCUIT_LU_STARTED",json.dumps({
        "source_regions":len(island_ids),"marine_nodes":len(giant),
        "total_source_nodes":len(node_components),"source_boundary_edges":int(ret.sum()),
        "verified_weight_nodes":uniq
    }),flush=True)
    resistance,qa=circuit.distributed_effective_resistance(lap,weights)
    if resistance.shape!=(104,104) or not np.isfinite(resistance).all() or
       not np.allclose(resistance,resistance.T,atol=1e-8,rtol=1e-8) or
       not np.allclose(np.diag(resistance),0,atol=1e-10):
        raise RuntimeError("Original exact 104-terminal uniform circuit resistance invalid")
    upper=resistance[np.triu_indices(104,k=1)]
    med=float(np.median(upper[upper>0]))
    if med<=0 or not math.isfinite(med):raise RuntimeError("Original ocean resistance median invalid")
    ix={code:i for i,code in enumerate(island_ids)}
    pairs=[]
    for r in original:
        ca,cb=r["source_region"],r["destination_region"]
        if ca not in ix or cb not in ix:continue
        value=float(resistance[ix[ca],ix[cb]])
        if value<=0:raise RuntimeError(f"Source WCVP island coastal terminal alias/source resistance zero {ca}/{cb}")
        pairs.append({
           "source_region":ca,"destination_region":cb,
           "source_verified_original_GSHHG_uniform_finite_volume_effective_resistance":value,
           "log1p_normalized_uniform_circuit":math.log1p(value/med),
           "source_original_TDWG3_coastal_anchor_count":len(anchored[ca]),
           "destination_original_TDWG3_coastal_anchor_count":len(anchored[cb]),
           "is_current_velocity_or_particle_transport":False,
           "is_historic_seed_or_ship_arrival":False
        })
    if len(pairs)!=5356:raise RuntimeError(f"Original 104 botanical source circuit pair count changed: {len(pairs)}")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"WCVP126_original104_uniform_finvol_electrical_resistance.csv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(pairs[0]));w.writeheader();w.writerows(pairs)
    receipt={
       "schema":"chocho_126_WCVP_island_original_GSHHG_i_0p0625_0p5_uniform_circuit_effective_resistance_v1",
       "status":"SOURCE_VERIFIED_ORIGINAL_UNIFORM_MARINE_CIRCUIT_104_REGIONS_NOT_OCEAN_CURRENTS",
       "original_GSHHG_source_graph_SHA256":SOURCE_GSHHG_SHA,
       "source_original_126_WCVP_island_regions":126,
       "source_verified_104_TDWG_coastal_region_terminals":len(island_ids),
       "source_original_104_unordered_pair_resistances":len(pairs),
       "source_component_marine_nodes":len(giant),
       "source_component_marine_edges":int(ret.sum()),
       "source_original_connected_GSHHG_graph_node_count":len(node_components),
       "conductance_method":"original island51 finite_volume_edges with uniform input 1: shared fine-water interface width / spherical center distance",
       "source_terminal_method":"equal-weighted up to 3 independently TDWG polygon GSHHG land-verified coastal-water source nodes, not WPI ports",
       "effective_resistance_normalization":"log1p(value / positive 104×104 median resistance)",
       "median_positive_uniform_resistance":med,
       "min_positive_uniform_resistance":float(np.min(upper[upper>0])),
       "max_uniform_resistance":float(upper.max()),
       "original_source_solver_numerical_gates":qa,
       "execution_elapsed_seconds":time.monotonic()-started,
       "limits":[
          "Original GSHHG 0.0625° fine raster and 0.5° component-aware block graph, so finite-volume weights are approximate discretized ocean, not continuum shoreline physics.",
          "The true uniform electrical effective resistance is distinct from shortest-water path length and from sea currents or directed marine transport.",
          "No island is used as a land/stepping-stone bridge, only coastal injection terminals independently fixed by source geometry.",
          "A region with no verified GSHHG shore anchor is missing source coverage, not a zero natural ocean connection.",
          "No historical current, cargo, shipping OD, viable seed arrival, seed buoyancy or butterfly establishment information enters circuit.",
       ]
    }
    (a.out_dir/"uniform104_WCVP_GSHHG_source_circuit_numerical_audit.json").write_text(
       json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("WCVP126_ORIGINAL_UNIFORM_CIRCUIT_FINISHED",json.dumps({
        "status":receipt["status"],"islands":len(island_ids),
        "pairs":len(pairs),"median_resistance":med,
        "max_relative_residual":qa["max_linear_system_relative_residual"],
        "seconds":receipt["execution_elapsed_seconds"]
    }),flush=True)

if __name__=="__main__":
    main()

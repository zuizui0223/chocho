#!/usr/bin/env python3
"""Map modelled over-water current corridors, NOT biological island relays.

Compute a source-locked, equal-pair-weighted mean of absolute edge currents
across all 1,275 original island pairs. A single grounded LU and 50
island-minus-reference solves suffice by electrical superposition.
"""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

REPO = Path(__file__).resolve().parents[1]


def pair_mean_absolute_edge_current(potential, edge_i, edge_j, conductance, chunk_size=16384):
    """Exact pairwise absolute flow average from sorted island voltage drops.

    Sum_(i<j) abs(x_i-x_j) = Sum_k (2k-n+1) sorted(x)_k.
    """
    p = np.asarray(potential, float)
    u, v = np.asarray(edge_i, int), np.asarray(edge_j, int)
    g = np.asarray(conductance, float)
    if p.ndim != 2 or p.shape[1] < 2 or u.shape != v.shape or u.shape != g.shape:
        raise ValueError("Inconsistent potentials and graph edge shapes")
    if (not np.isfinite(p).all() or not np.isfinite(g).all() or np.any(g <= 0) or
            np.any(u < 0) or np.any(v < 0) or
            np.any(u >= len(p)) or np.any(v >= len(p)) or np.any(u == v)):
        raise ValueError("Invalid source potentials, graph endpoints or conductance")
    n = p.shape[1]
    coeff = 2*np.arange(n,dtype=float)-(n-1)
    result = np.empty(len(g))
    for start in range(0,len(g),chunk_size):
        stop=min(start+chunk_size,len(g))
        drops=p[u[start:stop]]-p[v[start:stop]]
        drops.sort(axis=1)
        result[start:stop]=g[start:stop]*np.maximum(drops@coeff,0.)/(n*(n-1)//2)
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ValueError("Invalid pairwise mean currents")
    return result


def marine_network(topology, conductance, weights, component_label):
    component=np.asarray(topology["node_component"],int)
    marine=np.flatnonzero(component==int(component_label))
    remap=np.full(len(component),-1,dtype=int)
    remap[marine]=np.arange(len(marine))
    src=np.asarray(topology["edge_i"],int)
    dst=np.asarray(topology["edge_j"],int)
    conductance=np.asarray(conductance,float)
    if src.shape!=dst.shape or src.shape!=conductance.shape:
        raise ValueError("Source graph arrays differ in length")
    good=(remap[src]>=0)&(remap[dst]>=0)
    if np.any((remap[src]>=0)!=(remap[dst]>=0)):
        raise ValueError("Source graph bridges different components")
    u,v,g=remap[src[good]],remap[dst[good]],conductance[good]
    w=sparse.csr_matrix(weights[:,marine],dtype=float)
    if (w.shape[0]!=51 or not np.allclose(np.asarray(w.sum(axis=1)).ravel(),1.,atol=1e-12) or
            len(marine)!=189550 or len(g)!=370996 or
            not np.isfinite(g).all() or np.any(g<=0)):
        raise ValueError("Original frozen marine-source graph or source weights changed")
    return marine,u,v,g,w


def solve_island_potentials(node_count, u, v, g, weights, batch_size=4):
    """Ground one marine node, solve for w_i - w_0, and retain all 51 columns."""
    n=int(node_count)
    w=sparse.csr_matrix(weights,dtype=float)
    if n<2 or w.shape[1]!=n or w.shape[0]<2:
        raise ValueError("Invalid node/source count")
    u,v,g=np.asarray(u,int),np.asarray(v,int),np.asarray(g,float)
    adjacency=sparse.csr_matrix((np.r_[g,g],(np.r_[u,v],np.r_[v,u])),shape=(n,n))
    if connected_components(adjacency,directed=False,return_labels=False)!=1:
        raise ValueError("Marine graph disconnected")
    laplacian=sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency
    reduced=laplacian[1:,1:].tocsc()
    ref=sparse.vstack([w[0]]*(w.shape[0]-1),format="csr")
    injections=(w[1:]-ref).T.tocsr()[1:,:]
    factor=splu(reduced,permc_spec="MMD_AT_PLUS_A",diag_pivot_thresh=0.,
                options={"SymmetricMode":True})
    p=np.zeros((n,w.shape[0]))
    max_residual=0.
    for start in range(0,w.shape[0]-1,batch_size):
        stop=min(start+batch_size,w.shape[0]-1)
        rhs=injections[:,start:stop].toarray()
        volts=factor.solve(rhs)
        if not np.isfinite(volts).all():
            raise ValueError("Nonfinite solved node potentials")
        residual=reduced@volts-rhs
        den=np.maximum(np.max(np.abs(rhs),axis=0),1e-30)
        max_residual=max(max_residual,float(np.max(np.max(np.abs(residual),axis=0)/den)))
        p[1:,start+1:stop+1]=volts
    if max_residual>1e-7:
        raise ValueError("Electrical residual exceeds frozen tolerance")
    return p,max_residual


def node_through_current(edge_current,u,v,node_count):
    x=np.asarray(edge_current,float)
    return np.bincount(np.r_[u,v],weights=.5*np.r_[x,x],minlength=node_count)


def verify_against_resistance(p,u,v,g,ids,source):
    with np.load(source) as data:
        original=np.asarray(data["resistance"],float)
        archived_ids=[str(x) for x in data["island_ids"]]
    if archived_ids!=list(map(str,ids)) or original.shape!=(len(ids),len(ids)):
        raise ValueError("Archived resistance matrix/source IDs differ")
    reports=[]
    for i,j in [(0,1),(0,len(ids)-1),(len(ids)//2,len(ids)-1)]:
        volts=p[:,i]-p[:,j]
        energy=float(np.sum(g*(volts[u]-volts[v])**2))
        expected=float(original[i,j])
        relative=abs(energy-expected)/max(abs(expected),1e-12)
        reports.append({"a":str(ids[i]),"b":str(ids[j]),"energy":energy,
                        "archived":expected,"relative_error":relative})
        if not np.isfinite(relative) or relative>2e-6:
            raise ValueError("Current energy did not reproduce archived resistance")
    return reports


def original_island_centers(input_root):
    sys.path.insert(0,str(REPO/"scripts"))
    from reconstruct_island51_verified_port_coordinates import spherical_mean
    from run_island51_observed_ocean_circuits import SOURCE_DIRECTORY,validate_source_files
    root=Path(input_root)
    pairs,ids,groups,coords,_,_,_=validate_source_files(
        root/"wpi/ports_2019_source.geojson",
        root/"coordinates/island51_original_WPI_ports_coords_unique_source.json",
        SOURCE_DIRECTORY/"island_host_network_pairs.csv",
        SOURCE_DIRECTORY/"island_host_network_islands.csv")
    if len(pairs)!=1275:
        raise ValueError("Original pair table changed")
    return {key:spherical_mean([coords[port] for port in groups[key]]) for key in ids}


def render_map(raster,centers,out,field):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    positive=raster[np.isfinite(raster)&(raster>0)]
    if not len(positive):
        raise ValueError("No positive water currents")
    lower=max(float(np.quantile(positive,.05)),1e-12)
    upper=float(np.quantile(positive,.995))
    if upper<=lower:
        upper=float(positive.max())
    pixels=np.ma.masked_invalid(np.where(raster>0,raster,np.nan))
    fig,ax=plt.subplots(figsize=(15,6.8),constrained_layout=True)
    ax.set_facecolor("#e9e8e5")
    cmap=plt.get_cmap("magma").copy()
    cmap.set_bad("#e9e8e5")
    image=ax.imshow(pixels,extent=(-180,180,-85,85),origin="upper",
                    interpolation="nearest",aspect="auto",
                    cmap=cmap,norm=LogNorm(vmin=lower,vmax=upper))
    ax.scatter([c["longitude"] for c in centers.values()],
               [c["latitude"] for c in centers.values()],
               s=18,c="none",edgecolors="white",linewidths=.8,
               label="51 original WPI port-group centroids",zorder=3)
    ax.set(xlim=(-180,180),ylim=(-85,85),xlabel="Longitude",ylabel="Latitude")
    ax.set_xticks(np.arange(-180,181,60))
    ax.set_yticks(np.arange(-60,61,30))
    ax.grid(alpha=.12,linewidth=.5)
    ax.set_title("Potential over-water corridors among 51 original island groups",fontsize=14)
    ax.text(.5,-.11,field.upper()+" conductance | 1,275 equally weighted island pairs | no observed voyages or stopovers",
            transform=ax.transAxes,ha="center",va="top",fontsize=9)
    cb=fig.colorbar(image,ax=ax,fraction=.026,pad=.015,extend="max")
    cb.set_label("Mean absolute model current per unit-current island pair")
    ax.legend(loc="lower left",fontsize=8)
    out=Path(out)
    for ext in ("png","svg","pdf"):
        fig.savefig(out.with_suffix("."+ext),dpi=220)
    plt.close(fig)
    return {"color_scale_p05":lower,"color_scale_p995":upper,
            "positive_display_cells":int(len(positive))}


def run(input_root,out_dir,field="uniform"):
    if field not in ("uniform","ais","wind"):
        raise ValueError("Field must be a source-frozen ocean conductivity")
    sys.path.insert(0,str(REPO/"scripts"))
    from run_island51_observed_ocean_circuits import prepare_inputs
    root,out=Path(input_root),Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    _,ids,_,topology,edge_fields,weights,receipt=prepare_inputs(
        root,root/"coast_0p0625",out/"_source_receipt",out/"_network_exports",
        marine_ports_only=True)
    marine,u,v,g,w=marine_network(
        topology,edge_fields[field]["conductance"],weights,receipt["marine_component"])
    p,maxres=solve_island_potentials(len(marine),u,v,g,w)
    checks=verify_against_resistance(
        p,u,v,g,ids,
        REPO/"results/island51_observed_ocean_circuit_20261106"/
        (field+"_distributed_port_resistance.npz"))
    mean_edge=pair_mean_absolute_edge_current(p,u,v,g)
    current=node_through_current(mean_edge,u,v,len(marine))
    rr=np.asarray(topology["node_block_row"],int)[marine]
    cc=np.asarray(topology["node_block_col"],int)[marine]
    sums=np.zeros((340,720))
    hits=np.zeros((340,720),dtype=int)
    np.add.at(sums,(rr,cc),current)
    np.add.at(hits,(rr,cc),1)
    raster=np.where(hits>0,sums,np.nan)
    centers=original_island_centers(root)
    if set(centers)!=set(ids):
        raise ValueError("Original island coordinates differ from network endpoints")
    with (out/"original51_centroids.csv").open("w",newline="") as fh:
        writer=csv.writer(fh)
        writer.writerow(["landmass_id","WPI2019_longitude","WPI2019_latitude"])
        for key in ids:
            writer.writerow([key,centers[key]["longitude"],centers[key]["latitude"]])
    np.savez_compressed(out/(field+"_mean_pair_node_current.npz"),
                        raster=raster,marine_nodes=marine,node_current=current,
                        island_ids=np.array(ids))
    figure=render_map(raster,centers,out/("figure4_"+field+"_marine_current"),field)
    report={
      "schema":"chocho_island51_overwater_pair_current_map_v1",
      "conductivity":field,"islands":len(ids),
      "island_pairs":len(ids)*(len(ids)-1)//2,
      "marine_nodes":len(marine),"marine_edges":len(g),
      "original_marine_ports":int(receipt["marine_ports_included"]),
      "frozen_source_manifest_sha":receipt["frozen_source_bundle"].get("manifest_sha256"),
      "max_linear_residual":maxres,"archived_resistance_energy_checks":checks,
      "figure":figure,
      "interpretation":[
        "This is electrical corridor concentration averaged over modelled unit-current island pairs.",
        "No ship voyages, butterfly movement, successful feeding or colonization are measured.",
        "Ports on the same island are not electrically shorted and island passage is not modelled.",
        "High current near islands cannot identify stepping-stone hubs.",
        "No biological outcome was used to choose coastlines, conductance or colour scales."
      ]
    }
    (out/"execution_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
    (out/"README.md").write_text(
      "# Fig.4: source-locked over-water current corridors\n\n"
      "One unit of electrical current is injected and withdrawn among every "
      "unordered pair of 51 original island port groups. The figure plots the "
      "mean absolute edge current across all 1,275 pairs, with half of each "
      "edge attributed to each endpoint water node. The 0.5-degree figure "
      "cell adds currents of disconnected water patches in that display cell.\n\n"
      "The exact pair-average is calculated from 50 independent island-minus-reference "
      "potential solves and a sorted absolute-difference identity. Three pair "
      "energy checks are compared with the archived source-locked resistance matrix.\n\n"
      "**Not evidence for a stepping-stone island, ships or colonization:** islands "
      "are current-injection sites, not through-land transit nodes. Equal pair "
      "weights do not represent measured voyages or propagule pressure. "
      "No source outcomes, results or manuscripts are altered.\n\n"
      "## Execution receipt\n\n"+json.dumps(report,indent=2)+"\n")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root",type=Path,
        default=REPO/"data/exploratory/island51_observed_ocean_20261106")
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--field",choices=("uniform","ais","wind"),default="uniform")
    args=parser.parse_args()
    print("ISLAND51_MARINE_CURRENT",json.dumps(run(args.input_root,args.out_dir,args.field)),
          flush=True)


if __name__=="__main__":
    main()

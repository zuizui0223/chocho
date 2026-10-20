#!/usr/bin/env python3
"""Response-free contact of MODERN insect sampling anchors with 0BP/21ka
potential native larval-host climate components on ICE-6G_C ice-free terrain.

NOT actual historical insect sites, plant occupancy, dispersal routes, or
genetic barriers. Historical ice-covered anchors remain UNAVAILABLE and
unmatched endpoints are never silently declared disconnected.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from decimal import Decimal
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from pmip4_potential_habitat_components_v02 import components

SOURCE={
 "localities":"037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9",
 "edges":"ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9",
 "host_manifest":"a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044",
 "analog":"169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af",
 "extrapolated":"67c0d99852b8b8896559dbc0da07504bcb07a6ad7c017afa0eb86e31d40acc2a",
 "terrain":"d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db",
 "original513":"838b73df7f9ba7cd23cc33b981d3a434931ab1a6ed61f4bd27621f009d129020"}
SHAPE=(360,720)
GRID_N=SHAPE[0]*SHAPE[1]
R=6371.0088
RADII=(50,100,150)
EPOCHS=("0BP","21ka")
SCENARIOS=("analog","extrapolated")

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fh:
        for b in iter(lambda:fh.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open("r",newline="",encoding="utf8") as fh:
        return list(csv.DictReader(fh))

def original_id(value):
    d=Decimal(str(value).strip())
    if not d.is_finite() or d<=0 or d!=d.to_integral_value():
        raise ValueError("original WCVP accepted plant ID invalid")
    return str(int(d))

def latlon_to_cell(latlon):
    coords=np.asarray(latlon,float)
    if (coords.ndim!=2 or coords.shape[1]!=2 or
        not np.isfinite(coords).all() or np.any(abs(coords[:,0])>90)
        or np.any(abs(coords[:,1])>180)):
        raise ValueError("invalid original modern sampled locality coordinate")
    row=np.clip(np.floor((90-coords[:,0])/.5).astype(int),0,359)
    col=np.mod(np.floor((coords[:,1]+180)/.5).astype(int),720)
    return (row*720+col).astype(int)

def unit_vectors(latlon):
    xy=np.asarray(latlon,float)
    lat=np.deg2rad(xy[:,0]);lon=np.deg2rad(xy[:,1])
    c=np.cos(lat)
    return np.column_stack([c*np.cos(lon),c*np.sin(lon),np.sin(lat)])

def grid_cell_centers(indices):
    ind=np.asarray(indices,dtype=int)
    return np.column_stack([89.75-.5*(ind//720),-179.75+.5*(ind%720)])

def load_packed(path):
    with np.load(path,allow_pickle=False) as z:
        ids=[original_id(v) for v in z["plant_wcvp_ids"]]
        masks={epoch:np.asarray(z["mask_"+epoch],dtype="u1") for epoch in EPOCHS}
        if (len(ids)!=len(set(ids))!=416 or
            masks["0BP"].shape!=masks["21ka"].shape!=(416,32400)
            or tuple(z["grid_shape"].tolist())!=SHAPE
            or int(z["land_cut"])!=5
            or abs(float(z["threshold_quantile"])-.1)>1e-10):
            raise ValueError("source 416 plant mask identity/shape or threshold changed")
        excluded=bool(z["nonanalog_excluded"])
    return {ident:i for i,ident in enumerate(ids)},masks,excluded

def host_union(mask,host_indices):
    ind=np.asarray(host_indices,int)
    if len(ind)==0 or len(set(ind))!=len(ind):
        raise ValueError("known hosts missing or duplicate indices")
    arr=np.unpackbits(np.asarray(mask[ind]),axis=1)[:,:GRID_N]
    return np.any(arr.astype(bool),axis=0).reshape(SHAPE)

def build_terrain_contact_index(host_mask,phys_labels):
    """Per physically connected ice-free landmass: suitable-site tree and host patch ID."""
    host=np.asarray(host_mask,bool)
    physical=np.asarray(phys_labels,int)
    if not (host.shape==physical.shape==SHAPE):
        raise ValueError("inconsistent physical and host climate grid shape")
    if np.any(host & (physical==0)):
        raise ValueError("suitable climate cell found on source ocean/ice")
    host_labels,_=components(host,4)
    land_ids=np.unique(physical[host]); land_ids=land_ids[land_ids!=0]
    out={}
    flatphys=physical.reshape(-1)
    cells=np.flatnonzero(host)
    for land in land_ids:
        chosen=cells[flatphys[cells]==land]
        xyz=unit_vectors(grid_cell_centers(chosen))
        out[int(land)]=(cKDTree(xyz),chosen,host_labels.reshape(-1)[chosen])
    return out

def contact_modern_anchors(coords,phys_labels,source):
    """Nearest host cell on SAME physical landmass, regardless of geographic edge response."""
    positions=np.asarray(coords,float)
    indices=latlon_to_cell(positions)
    original_physical=np.asarray(phys_labels).reshape(-1)[indices]
    nearest_km=np.full(len(indices),np.inf,float)
    habitat_comp=np.zeros(len(indices),dtype=np.int32)
    state=np.full(len(indices),"PHYSICAL_ANCHOR_BLOCKED",dtype=object)
    for land in np.unique(original_physical):
        if land==0:continue
        ix=np.flatnonzero(original_physical==land)
        if land not in source:
            state[ix]="NO_HOST_POTENTIAL_ON_SAME_LANDMASS"
            continue
        tree,cells,hostlabels=source[int(land)]
        d,index=tree.query(unit_vectors(positions[ix]),k=1,workers=1)
        nearest_km[ix]=2*R*np.arcsin(np.clip(d/2,0,1))
        habitat_comp[ix]=hostlabels[np.asarray(index,int)]
        state[ix]="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE"
    return nearest_km,habitat_comp,state

def contact_summary(node_ids,coords,edges,nearest_km,hostcomp,state,radius_km):
    """No comparable edge is fabricated from a physically blocked/unmatched endpoint."""
    eligible=(state=="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE")&(nearest_km<=radius_km)
    ids={i:j for j,i in enumerate(node_ids)}
    if len(ids)!=len(node_ids):raise ValueError("duplicate sample locality index")
    pair=np.array([(ids[a],ids[b]) for a,b in edges],int)
    u,v=pair.T
    both=eligible[u]&eligible[v]
    one=eligible[u]^eligible[v]
    same=both&(hostcomp[u]==hostcomp[v])
    diff=both&~same
    n=len(edges)
    both_count=int(both.sum())
    return {
      "locality_n":int(len(node_ids)),
      "edges_n":int(n),
      "physically_blocked_locality_n":int(np.count_nonzero(state=="PHYSICAL_ANCHOR_BLOCKED")),
      "same_landmass_no_climate_locality_n":int(np.count_nonzero(state=="NO_HOST_POTENTIAL_ON_SAME_LANDMASS")),
      "eligible_locality_n":int(eligible.sum()),
      "land_anchor_unmatched_radius_n":int(np.count_nonzero(
           (state=="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE")&~eligible)),
      "edge_both_anchor_contacts_n":both_count,
      "edge_one_anchor_contact_n":int(one.sum()),
      "edge_neither_anchor_contact_n":int(n-both_count-one.sum()),
      "edge_same_host_climate_component_n":int(same.sum()),
      "edge_different_host_climate_components_n":int(diff.sum()),
      "edge_same_conditional_fraction":float(same.sum()/both_count) if both_count else None,
      "edge_same_all_original_fraction":float(same.sum()/n) if n else None,
      "anchor_contact_median_distance_km":float(np.median(nearest_km[eligible])) if eligible.any() else None
    }

def create_analysis(args):
    for k,path in (("localities",args.localities),("edges",args.edges),
                   ("host_manifest",args.manifest),("analog",args.analog),
                   ("extrapolated",args.extrapolated),("terrain",args.terrain),
                   ("original513",args.original513)):
        if digest(path)!=SOURCE[k]:raise ValueError("original response-free source SHA drift: "+k)
    links=defaultdict(set)
    for record in read_csv(args.manifest):
        links[record["insect_species"]].add(original_id(record["wcvp_id_as_frozen"]))
    original=read_csv(args.original513)
    complete={r["insect_species"] for r in original
               if r["status"]=="ALL_RECORDED_HOST_MODELS_AVAILABLE"}
    if len(links)!=641 or len(complete)!=513:
        raise ValueError("641 original species or 513 host-ready identities changed")
    loc=defaultdict(dict);edges=defaultdict(list)
    for record in read_csv(args.localities):
        loc[record["species"]][int(record["locality_index"])]=(
           float(record["latitude"]),float(record["longitude"]))
    for record in read_csv(args.edges):
        edges[record["species"]].append((
          int(record["node_left"]),int(record["node_right"])))
    if not (len(loc)==len(edges)==642) or not complete<=set(loc):
        raise ValueError("original 642 geometry taxa and 513 source insects don't align")
    groups=defaultdict(list)
    for species in complete:
        groups[tuple(sorted(links[species],key=int))].append(species)
    if len(groups)!=349:raise ValueError("349 original source-known host combinations changed")
    analog_idx,analog_masks,analog_censor=load_packed(args.analog)
    alt_idx,alt_masks,alt_censor=load_packed(args.extrapolated)
    if analog_idx!=alt_idx or not analog_censor or alt_censor:
        raise ValueError("main and climate-nonanalog sensitivity plant IDs/censor status changed")
    for epoch in EPOCHS:
        a=analog_masks[epoch];b=alt_masks[epoch]
        if np.any(np.bitwise_and(a,np.bitwise_not(b))):
            raise ValueError("extrapolated mask removed an original supported pixel")
    with np.load(args.terrain,allow_pickle=False) as grid:
        terrain={epoch:(np.asarray(grid["pmip"+("0" if epoch=="0BP" else "21")+"_icefree_subcells"])>=5)
                 for epoch in EPOCHS}
        if any(a.shape!=SHAPE for a in terrain.values()):
            raise ValueError("PMIP4 frozen majority terrestrial grid misalignment")
    physical_labels={ep:components(terrain[ep],4)[0] for ep in EPOCHS}
    rows=[]
    for gi,(hosts,insects) in enumerate(sorted(groups.items())):
        indices=[analog_idx[h] for h in hosts]
        for scenario in SCENARIOS:
            allm=analog_masks if scenario=="analog" else alt_masks
            for epoch in EPOCHS:
                merged=host_union(allm[epoch],indices)
                if np.any(merged & ~terrain[epoch]):
                    raise ValueError("source host climate graph includes physical ocean/ice")
                ci=build_terrain_contact_index(merged,physical_labels[epoch])
                for sp in sorted(insects):
                    order=sorted(loc[sp])
                    coords=[loc[sp][key] for key in order]
                    dist,comp,state=contact_modern_anchors(coords,physical_labels[epoch],ci)
                    for radius in RADII:
                        result=contact_summary(order,coords,edges[sp],dist,comp,state,radius)
                        rows.append({"species":sp,"hostset":"/".join(hosts),
                            "shared_hostset_insects":len(insects),
                            "host_count":len(hosts),"epoch":epoch,
                            "scenario":scenario,"radius_km":radius,**result})
        if (gi+1)%50==0:
            print(json.dumps({"source_hostsets_processed":gi+1,"out_of":349}),flush=True)
    if len(rows)!=513*len(EPOCHS)*len(SCENARIOS)*len(RADII):
        raise ValueError("all species and 12 combinations were not scored")
    def historical_differences(scenario,radius):
        matching={}
        for r in rows:
            if r["scenario"]==scenario and r["radius_km"]==radius:
                matching[(r["species"],r["epoch"])]=r
        qualified=[]
        for sp in complete:
            a,b=matching[(sp,"0BP")],matching[(sp,"21ka")]
            if a["edge_both_anchor_contacts_n"] and b["edge_both_anchor_contacts_n"]:
                qualified.append(b["edge_same_conditional_fraction"]-a["edge_same_conditional_fraction"])
        return {"eligible_insects_with_comparable_both_era_pairs":len(qualified),
                "share_same_patch_decrease":int(sum(x<0 for x in qualified)),
                "share_same_patch_increase":int(sum(x>0 for x in qualified)),
                "share_same_patch_equal":int(sum(x==0 for x in qualified)),
                "median_within_species_conditional_difference":float(np.median(qualified)) if qualified else None}
    summary={"schema":"chocho_LGM_513_modern_geographic_anchor_potential_host_contact_v05",
        "status":"EXPLORATORY_MODERN_ANCHOR_CONTACT_NOT_ANCIENT_OCCUPANCY_OR_GENETIC_BARRIER",
        "original_animal_species":641,"host_complete_insects":513,"distinct_recorded_host_combinations":349,
        "source_rows":len(rows),"source_original_geometry_taxa":642,
        "comparisons":{f"{mode}_radius{radius}km":historical_differences(mode,radius)
            for mode in SCENARIOS for radius in RADII},
        "source_sha256":dict(SOURCE),
        "limits":["Modern insect sampled localities are artificially fixed geographic anchors in a 21ka climate, NOT inferred ancient observations",
            "Assignment requires contact with a potential host climatic patch center on the same PMIP4 ice-free coarse landmass within threshold radius",
            "Contact on same physical landmass does NOT ensure true access from an anchor along a continuously plant-host-occupied route",
            "Unmatched modern localities at 21ka are unavailable/unknown, not automatically genetic isolation",
            "Same-hostset 349 landscapes are shared among 513 insects, within-species edges are correlated and original genetic distances were never opened",
            "Negative TTF host-identity climate proxy remains frozen; no host-mediated genetic mechanism supported"]}
    return rows,summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("localities","edges","manifest","analog","extrapolated","terrain","original513","out"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    rows,q=create_analysis(a)
    out=a.out/"all_513_modern_anchor_host_climate_contact_v05.csv"
    with out.open("w",newline="",encoding="utf8") as file:
        w=csv.DictWriter(file,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader();w.writerows(sorted(rows,key=lambda r:(r["species"],r["scenario"],r["radius_km"],r["epoch"])))
    (a.out/"summary.json").write_text(json.dumps(q,indent=2,sort_keys=True)+"\n")
    print(json.dumps(q,sort_keys=True))
if __name__=="__main__":main()

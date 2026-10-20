#!/usr/bin/env python3
"""Compare EXACT SAME original modern insect sampling edges in two glacial climates.

A modern sampling edge is evaluable only if BOTH endpoints can be assigned to a
source-pinned potential host climatic patch in BOTH 0BP and 21ka. Source PMIP4
ice, landmass and 50/100/150 km snapping must match prior v05 results exactly.

This is modern geographic sampling-anchor counterfactual, NOT historic insect
or plant occupation, actual larval dispersal, or an observed genetic barrier.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from pmip4_potential_habitat_components_v02 import components
from audit_modern_anchor_paleo_host_graph_v05 import (
    SOURCE,SHAPE,EPOCHS,SCENARIOS,RADII,read_csv,original_id,load_packed,
    host_union,build_terrain_contact_index,contact_modern_anchors,contact_summary
)
PRIOR_V05_SHA="6272f2ac132ef3713356a4566b8b21007343931e0e05cd967150343d4b6194b7"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as inp:
        for chunk in iter(lambda:inp.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def source_qualified(original513,links,loc,edges):
    ready={r["insect_species"] for r in read_csv(original513)
           if r["status"]=="ALL_RECORDED_HOST_MODELS_AVAILABLE"}
    if (len(ready)!=513 or len(links)!=641 or
        not ready<=set(loc) or not (len(loc)==len(edges)==642)):
        raise ValueError("original 641 species/642 geo source 513 fully-host-modelled cohort drift")
    return ready

def paired_transition_counts(ordered_nodes,pairs,a,b,radius):
    """Same edge denominator in both ages; zero shared edges becomes None."""
    ids={k:i for i,k in enumerate(ordered_nodes)}
    if len(ids)!=len(ordered_nodes):
        raise ValueError("duplicate original modern sample locality")
    edges=np.array([(ids[x],ids[y]) for x,y in pairs],int)
    u,v=edges.T
    a_dist,a_comp,a_status=a
    b_dist,b_comp,b_status=b
    ok_a=(a_status=="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE")&(a_dist<=radius)
    ok_b=(b_status=="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE")&(b_dist<=radius)
    full_a=ok_a[u]&ok_a[v]
    full_b=ok_b[u]&ok_b[v]
    both=full_a&full_b
    same_a=(a_comp[u]==a_comp[v])
    same_b=(b_comp[u]==b_comp[v])
    n=int(both.sum())
    s2s=int(np.count_nonzero(both&same_a&same_b))
    s2d=int(np.count_nonzero(both&same_a&~same_b))
    d2s=int(np.count_nonzero(both&~same_a&same_b))
    d2d=int(np.count_nonzero(both&~same_a&~same_b))
    if s2s+s2d+d2s+d2d!=n:raise ValueError("incomplete paired exact-edge comparison")
    return {
      "all_original_edges_n":len(pairs),
      "both_endpoints_available_at_0BP_n":int(full_a.sum()),
      "both_endpoints_available_at_21ka_n":int(full_b.sum()),
      "edge_contactable_in_both_eras_n":n,
      "edge_contactable_only_0BP_n":int(np.count_nonzero(full_a&~full_b)),
      "edge_contactable_only_21ka_n":int(np.count_nonzero(~full_a&full_b)),
      "edge_contactable_neither_n":int(np.count_nonzero(~full_a&~full_b)),
      "same_to_same":s2s,"same_to_different":s2d,
      "different_to_same":d2s,"different_to_different":d2d,
      "same_host_component_fraction_0BP_on_exact_common_pairs":float((s2s+s2d)/n) if n else None,
      "same_host_component_fraction_21ka_on_exact_common_pairs":float((s2s+d2s)/n) if n else None,
      "same_component_paired_difference":float((d2s-s2d)/n) if n else None
    }

def load_input(args):
    for k,path in (("localities",args.localities),("edges",args.edges),
                   ("host_manifest",args.manifest),("analog",args.analog),
                   ("extrapolated",args.extrapolated),("terrain",args.terrain),
                   ("original513",args.original513)):
        expected=SOURCE["host_manifest" if k=="host_manifest" else k]
        if sha(path)!=expected:raise ValueError("original response-free source SHA mismatch: "+k)
    if sha(args.prior_v05)!=PRIOR_V05_SHA:
        raise ValueError("source v05 513 modern-anchor classification table modified")
    previous={(r["species"],r["scenario"],r["epoch"],int(r["radius_km"])):r
              for r in read_csv(args.prior_v05)}
    if len(previous)!=6156:
        raise ValueError("prior 513x2x2x3 contact classification count drift")
    links=defaultdict(set)
    for r in read_csv(args.manifest):
        links[r["insect_species"]].add(original_id(r["wcvp_id_as_frozen"]))
    loc=defaultdict(dict);edges=defaultdict(list)
    for r in read_csv(args.localities):
        loc[r["species"]][int(r["locality_index"])]=(float(r["latitude"]),float(r["longitude"]))
    for r in read_csv(args.edges):
        edges[r["species"]].append((int(r["node_left"]),int(r["node_right"])))
    ready=source_qualified(args.original513,links,loc,edges)
    grouped=defaultdict(list)
    for species in ready:
        grouped[tuple(sorted(links[species],key=int))].append(species)
    if len(grouped)!=349:
        raise ValueError("original 349 source host combinations changed")
    plant_ids,analog,main_flag=load_packed(args.analog)
    alt_ids,extrap,alt_flag=load_packed(args.extrapolated)
    if plant_ids!=alt_ids or not main_flag or alt_flag:
        raise ValueError("source main/extrapolated original 416 plant order changed")
    for ep in EPOCHS:
        if np.any(analog[ep]&~extrap[ep]):
            raise ValueError("source extrapolated habitat no longer contains previous analog-only habitat")
    with np.load(args.terrain,allow_pickle=False) as z:
        physical={
          "0BP":np.asarray(z["pmip0_icefree_subcells"])>=5,
          "21ka":np.asarray(z["pmip21_icefree_subcells"])>=5}
    if any(a.shape!=SHAPE for a in physical.values()):
        raise ValueError("original PMIP4 0/21ka ice-free physical terrain mask invalid")
    landcomps={e:components(physical[e],4)[0] for e in EPOCHS}
    return previous,links,loc,edges,ready,grouped,plant_ids,analog,extrap,landcomps

def audit(args):
    previous,links,loc,edges,ready,groups,plants,analog,extrap,landcomps=load_input(args)
    rows=[];prior_checked=0
    for j,(hosts,insects) in enumerate(sorted(groups.items())):
        pi=[plants[h] for h in hosts]
        for scen in SCENARIOS:
            source=analog if scen=="analog" else extrap
            # Patch graph and physical landmass must be computed independently per epoch.
            contact_index={}
            for epoch in EPOCHS:
                landscape=host_union(source[epoch],pi)
                contact_index[epoch]=build_terrain_contact_index(landscape,landcomps[epoch])
            for species in sorted(insects):
                order=sorted(loc[species]);positions=[loc[species][i] for i in order]
                arr={}
                for epoch in EPOCHS:
                    arr[epoch]=contact_modern_anchors(positions,landcomps[epoch],contact_index[epoch])
                for radius in RADII:
                    for epoch in EPOCHS:
                        check=contact_summary(order,positions,edges[species],*arr[epoch],radius)
                        original=previous[(species,scen,epoch,radius)]
                        fields=["eligible_locality_n","physically_blocked_locality_n",
                                "edge_both_anchor_contacts_n","edge_one_anchor_contact_n",
                                "edge_neither_anchor_contact_n","edge_same_host_climate_component_n",
                                "edge_different_host_climate_components_n"]
                        if any(int(check[k])!=int(original[k]) for k in fields):
                            raise ValueError("source v05 per-epoch original site/edge counts NOT reproduced: "+species)
                        prior_checked+=1
                    pair=paired_transition_counts(order,edges[species],
                        arr["0BP"],arr["21ka"],radius)
                    rows.append({"insect_species":species,"hostset":"/".join(hosts),
                                 "recorded_host_count":len(hosts),
                                 "insects_with_same_hostset":len(insects),
                                 "scenario":scen,"radius_km":radius,**pair})
        if (j+1)%75==0:
            print(json.dumps({"complete_source_hostsets":j+1,"total":349}),flush=True)
    if len(rows)!=513*2*3 or prior_checked!=6156:
        raise ValueError("source 513 species and exact 6156 v05 comparison rows not fully replayed")
    summaries={}
    for scen in SCENARIOS:
        for radius in RADII:
            case=[r for r in rows if r["scenario"]==scen and r["radius_km"]==radius]
            good=[r for r in case if r["edge_contactable_in_both_eras_n"]>0]
            robust=[r for r in good if r["edge_contactable_in_both_eras_n"]>=5]
            d=[r["same_component_paired_difference"] for r in good]
            counts=Counter()
            for r in good:
                if r["same_to_different"]>r["different_to_same"]:counts["more_same_to_different"]+=1
                elif r["same_to_different"]<r["different_to_same"]:counts["more_different_to_same"]+=1
                else:counts["equal_or_no_change"]+=1
            summaries[f"{scen}_radius{radius}km"]={
                "species":len(case),
                "species_with_common_contactable_original_edges":len(good),
                "species_with_at_least_five_common_edges":len(robust),
                "unique_hostsets_with_common_edges":len({r["hostset"] for r in good}),
                "original_edges_total":int(sum(r["all_original_edges_n"] for r in case)),
                "common_edges_total":int(sum(r["edge_contactable_in_both_eras_n"] for r in case)),
                "same_to_same_edges":int(sum(r["same_to_same"] for r in good)),
                "same_to_different_edges":int(sum(r["same_to_different"] for r in good)),
                "different_to_same_edges":int(sum(r["different_to_same"] for r in good)),
                "different_to_different_edges":int(sum(r["different_to_different"] for r in good)),
                "species_sign_counts":dict(counts),
                "median_species_within_same_edge_difference":float(np.median(d)) if d else None,
                "median_difference_species_with_at_least_five_edges":float(np.median(
                    [r["same_component_paired_difference"] for r in robust])) if robust else None
            }
    return rows,{"schema":"chocho_513_original_MODERN_anchor_exact_common_edges_0BP_21ka_v06",
        "status":"EXPLORATORY_PAIRED_MODERN_ANCHOR_POTENTIAL_PATCHES_NOT_REAL_GLACIAL_GENETIC_CONNECTIVITY",
        "source_original_species":641,"source_complete_host_insects":513,
        "source_shared_hostsets":349,
        "source_v05_per_epoch_rows_exactly_replayed":6156,
        "source_513x2x3_pair_records":len(rows),
        "matched_edge_comparisons":summaries,
        "source_sha256":dict(SOURCE,prior_v05=PRIOR_V05_SHA),
        "limitations":[
          "Anchors are MODERN insect sample coordinates artificially projected onto glacial landscapes, not actual ancient populations",
          "Both-epoch contactability conditions on a selected subset of original pairs and excludes glacier-blocked and uncontactable anchors",
          "A patch split is a difference in predicted climatic host habitat components, not measured larval resource continuity or gene flow",
          "Same 349 known hostsets are shared across 513 insects; edges within species are correlated, and genetic response was never opened",
          "Sources include coarse 0.5degree land-ice topology and model uncertainty; source TTF v0.5 true-host specificity result remains negative"]}
        
def main():
    p=argparse.ArgumentParser()
    for name in ("localities","edges","manifest","analog","extrapolated","terrain",
                 "original513","prior-v05","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    rows,out=audit(a)
    with (a.out/"all_513_original_edges_paired_contact_v06.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader()
        w.writerows(sorted(rows,key=lambda r:(r["insect_species"],r["scenario"],r["radius_km"])))
    (a.out/"summary.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,sort_keys=True))
if __name__=="__main__":main()

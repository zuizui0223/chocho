#!/usr/bin/env python3
"""Exploratory source-bound OTHER REAL PLANT null for glacial host-climate patch contact.

For each insect, compare its recorded native-host climate landscape with 4-8
OTHER REAL PLANT groups frozen by an earlier 0BP-only matching program.
Every actual-vs-alternative contrast uses EXACT SAME modern sampling pairs,
whose endpoints contact hypothetical host/other-plant climate patches in BOTH
0BP and 21ka, on original PMIP4 physically ice-free land, never across sea.

This is not ancient plant occupancy, insect dispersal, or genetic response.
The original TTF v0.5 identity control remains negative/frozen.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from pmip4_potential_habitat_components_v02 import components
from audit_modern_anchor_paleo_host_graph_v05 import (
  SHAPE,EPOCHS,SCENARIOS,RADII,read_csv,load_packed,host_union,
  build_terrain_contact_index,contact_modern_anchors
)
from audit_common_modern_sampling_edges_v06 import paired_transition_counts

SHA={
"model_plan":"bb791eff085ba46ba6012efff0310a3a01d88144f4d0d18ef5b1160ff2c0aa10",
"old_v06":"7b5a30433d0d20107041a24b22bcd6d19f4f708e8552cfc553f8afb105559d09",
"original_geo_nodes":"037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9",
"original_geo_edges":"ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9",
"old_primary_416":"169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af",
"old_extrap_416":"67c0d99852b8b8896559dbc0da07504bcb07a6ad7c017afa0eb86e31d40acc2a",
"PMIP4_05deg":"d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db"
}
SOURCE_CONTACT="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE"

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for buf in iter(lambda:stream.read(1024*1024),b""):h.update(buf)
    return h.hexdigest()

def full_quad_comparison(nodes,edges,actual0,actual21,alternate0,alternate21,radius):
    """One common denominator: BOTH endpoints contact all FOUR landscapes."""
    positions={k:i for i,k in enumerate(nodes)}
    if len(positions)!=len(nodes) or not edges:
        raise ValueError("original source locality IDs not valid")
    uv=np.array([(positions[a],positions[b]) for a,b in edges],dtype=int)
    u,v=uv.T
    def usable(pair):
        distances,labels,states=pair
        if not(len(distances)==len(labels)==len(states)==len(nodes)):
            raise ValueError("missing source per-anchor contact classification")
        return (np.asarray(states)==SOURCE_CONTACT)&(np.asarray(distances)<=radius)
    choices=[usable(s) for s in (actual0,actual21,alternate0,alternate21)]
    edge_eligible=np.ones(len(edges),bool)
    for ok in choices:
        edge_eligible &= ok[u]&ok[v]
    n=int(edge_eligible.sum())
    result={"all_original_edges_n":len(edges),
            "fourway_exact_common_original_edges_n":n}
    if n==0:
        result.update({"actual_0BP_same_fraction":None,"actual_21ka_same_fraction":None,
            "alternate_0BP_same_fraction":None,"alternate_21ka_same_fraction":None,
            "actual_paired_delta":None,"alternate_paired_delta":None,
            "actual_minus_alternate_paired_delta":None,
            "modern_0BP_baseline_same_patch_gap":None})
        return result
    actual_same0=np.asarray(actual0[1])[u]==np.asarray(actual0[1])[v]
    actual_same21=np.asarray(actual21[1])[u]==np.asarray(actual21[1])[v]
    other_same0=np.asarray(alternate0[1])[u]==np.asarray(alternate0[1])[v]
    other_same21=np.asarray(alternate21[1])[u]==np.asarray(alternate21[1])[v]
    m=edge_eligible
    actual0_s=int(np.count_nonzero(m&actual_same0))
    actual21_s=int(np.count_nonzero(m&actual_same21))
    other0_s=int(np.count_nonzero(m&other_same0))
    other21_s=int(np.count_nonzero(m&other_same21))
    da=(actual21_s-actual0_s)/n
    dc=(other21_s-other0_s)/n
    result.update({"actual_0BP_same_fraction":actual0_s/n,
        "actual_21ka_same_fraction":actual21_s/n,
        "alternate_0BP_same_fraction":other0_s/n,
        "alternate_21ka_same_fraction":other21_s/n,
        "actual_paired_delta":da,"alternate_paired_delta":dc,
        "actual_minus_alternate_paired_delta":da-dc,
        "modern_0BP_baseline_same_patch_gap":(actual0_s-other0_s)/n,
        "actual_same_to_different_n":int(np.count_nonzero(m&actual_same0&~actual_same21)),
        "actual_different_to_same_n":int(np.count_nonzero(m&~actual_same0&actual_same21)),
        "alternate_same_to_different_n":int(np.count_nonzero(m&other_same0&~other_same21)),
        "alternate_different_to_same_n":int(np.count_nonzero(m&~other_same0&other_same21))})
    return result

def source_load(args):
    for name,path in (("model_plan",args.model_plan),("old_v06",args.old_v06),
                      ("original_geo_nodes",args.localities),("original_geo_edges",args.edges),
                      ("old_primary_416",args.primary_mask),("old_extrap_416",args.extrap_mask),
                      ("PMIP4_05deg",args.terrain)):
        if digest(path)!=SHA[name]:raise ValueError("original source SHA changed: "+name)
    plan=json.loads(args.model_plan.read_text())
    selected=plan["matched_species"]
    if (plan["status"]!="ELIGIBILITY_ONLY_NO_NEW_LGM_PATCH_OUTCOME_READ"
        or not plan["no_new_21ka_climate_or_genetic_patch_response_accessed"]
        or len(selected)!=292 or plan["valid_frozen_alternative_group_count"]!=2108):
        raise ValueError("pre-glacial model-eligibility source qualified only 292 insects")
    old={(r["insect_species"],r["scenario"],int(r["radius_km"])):r for r in read_csv(args.old_v06)}
    if len(old)!=513*len(SCENARIOS)*len(RADII):
        raise ValueError("source v06 exact same-edge 513 original insects/age/radii drift")
    loc=defaultdict(dict);edges=defaultdict(list)
    for r in read_csv(args.localities):
        loc[r["species"]][int(r["locality_index"])]=(float(r["latitude"]),float(r["longitude"]))
    for r in read_csv(args.edges):
        edges[r["species"]].append((int(r["node_left"]),int(r["node_right"])))
    if not(len(loc)==len(edges)==642):
        raise ValueError("original geographic 642 taxa source changed")
    if not set(selected)<=set(loc):
        raise ValueError("not all source-qualified insects have original modern geographic anchors")
    modelids,main,mflag=load_packed(args.primary_mask)
    altids,extended,eflag=load_packed(args.extrap_mask)
    if modelids!=altids or not mflag or eflag or len(modelids)!=416:
        raise ValueError("original main/extrapolated WCVP plant model index changed")
    all_groups=sum(len(j["alternative_groups"]) for j in selected.values())
    if all_groups!=2108:raise ValueError("source alternative group count changed")
    with np.load(args.terrain,allow_pickle=False) as src:
        land={"0BP":np.asarray(src["pmip0_icefree_subcells"])>=5,
              "21ka":np.asarray(src["pmip21_icefree_subcells"])>=5}
    if any(x.shape!=SHAPE for x in land.values()):
        raise ValueError("source PMIP4 physically icefree grid changed")
    physical={e:components(land[e],4)[0] for e in EPOCHS}
    return selected,old,loc,edges,modelids,main,extended,physical

def analysis(args):
    selected,prior,loc,edges,ids,original,extrap,physical=source_load(args)
    groups=[];species_rows=[];replayed=0
    for n,(species,entry) in enumerate(sorted(selected.items())):
        indices=[ids[v] for v in entry["known_host_taxa"]]
        order=sorted(loc[species]);xy=[loc[species][i] for i in order]
        for scenario,source in (("analog",original),("extrapolated",extrap)):
            real={}
            for epoch in EPOCHS:
                habitat=host_union(source[epoch],indices)
                idx=build_terrain_contact_index(habitat,physical[epoch])
                real[epoch]=contact_modern_anchors(xy,physical[epoch],idx)
            # Original real-host results must first match exactly v06 source.
            for radius in RADII:
                old=prior[(species,scenario,radius)]
                got=paired_transition_counts(order,edges[species],real["0BP"],real["21ka"],radius)
                for key in ("edge_contactable_in_both_eras_n","same_to_same","same_to_different",
                            "different_to_same","different_to_different"):
                    if int(got[key])!=int(old[key]):
                        raise ValueError("source original host-control same-edge v06 replay failed "+species)
                replayed+=1
            by_radius=defaultdict(list)
            for g in entry["alternative_groups"]:
                donor_ids=g["donor_taxa"]
                if set(donor_ids)&set(entry["known_host_taxa"]):
                    raise ValueError("original source alternative includes known recorded host")
                control={}
                for epoch in EPOCHS:
                    habitat=host_union(source[epoch],[ids[v] for v in donor_ids])
                    idx=build_terrain_contact_index(habitat,physical[epoch])
                    control[epoch]=contact_modern_anchors(xy,physical[epoch],idx)
                for radius in RADII:
                    measured=full_quad_comparison(order,edges[species],
                        real["0BP"],real["21ka"],control["0BP"],control["21ka"],radius)
                    groups.append({"insect_species":species,"scenario":scenario,
                        "radius_km":radius,"other_plant_group_index":g["group_index"],
                        "other_plant_ids":";".join(donor_ids),
                        "known_actual_host_ids":";".join(entry["known_host_taxa"]),
                        "original_0BP_climate_moment_gap":g["modern_climate_moment_gap"],
                        **measured})
                    by_radius[radius].append(measured)
            for radius in RADII:
                options=by_radius[radius]
                good=[x for x in options if x["fourway_exact_common_original_edges_n"]>0]
                stable=[x for x in good if x["fourway_exact_common_original_edges_n"]>=5]
                summary={"insect_species":species,"scenario":scenario,"radius_km":radius,
                   "source_frozen_other_plant_group_n":len(options),
                   "other_plant_groups_with_common_original_edges_n":len(good),
                   "other_plant_groups_with_at_least_5_common_edges_n":len(stable),
                   "source_original_edge_n":len(edges[species]),
                   "status":"COMPUTABLE_PAIRWISE_OTHER_PLANT_NEGATIVE_CONTROL" if good else
                            "NOT_EVALUABLE_FOURWAY_COMMON_EDGE_DENOMINATOR"}
                if good:
                    diff=np.array([x["actual_minus_alternate_paired_delta"] for x in good])
                    summary.update({"median_actual_minus_otherplant_same_edge_difference":float(np.median(diff)),
                        "median_0BP_actual_vs_otherplant_baseline_gap":float(np.median(
                         [x["modern_0BP_baseline_same_patch_gap"] for x in good])),
                        "median_fourway_same_edge_n_per_group":float(np.median(
                         [x["fourway_exact_common_original_edges_n"] for x in good])),
                        "median_difference_on_groups_with_5plus_edges":float(np.median(
                         [x["actual_minus_alternate_paired_delta"] for x in stable])) if stable else None})
                species_rows.append(summary)
        if (n+1)%35==0:
            print(json.dumps({"source_original_insects_completed":n+1,
                "source_target":292,"source_v06_host_replay_rows":replayed}),flush=True)
    if len(groups)!=2108*2*3 or len(species_rows)!=292*2*3 or replayed!=292*2*3:
        raise ValueError("source 292/2108/3 radii/2 scenarios did not finish")
    result={"schema":"chocho_v08_original_0BP_matched_real_OTHER_plants_21ka_patch_sameedges",
       "status":"EXPLORATORY_HOST_SPECIFIC_PATCH_CLIMATE_NEGATIVE_CONTROL_NOT_GENETIC_OR_OCCUPANCY",
       "original_panel_species":641,"prior_0BP_matched_real_other_plant_insects":335,
       "source_qualified_actual_hosts_plus_four_otherplant_groups":292,
       "original_0BP_only_altplant_control_group_count":2108,
       "distinct_source_416_model_native_plant_taxa":416,
       "source_replayed_original_host_v06_scenario_radius_rows":replayed,
       "source_group_scenario_radius_rows":len(groups),
       "source_species_scenario_radius_rows":len(species_rows),
       "source_sha256":dict(SHA),
       "by_scenario_radius":{},
       "limits":[
        "Exact original modern anchor pairs must be contactable in 0 BP and21ka for BOTH the actual known host and the preselected other real plants, selecting a possibly small/nonrepresentative subset",
        "0BP-matched unlisted real plants are NOT biologically proven unused hosts; modern WCVP native-region model training can differ greatly between species",
        "Source other plants were chosen in frozen earlier 0 BP modern-climate moments; fitted 0 BP climatic patch structure may remain imbalanced and is explicitly reported",
        "Repeated 4-8 plant groups are not independent insect species; 292 insects also share known host plant combinations and climate source data",
        "This is not actual 21ka plant occupancy, observed ancestral insect localities, dispersal or genetic breaks; old TTF v0.5 real-other-plant negative control remains unchanged"]}
    for scenario in SCENARIOS:
        for radius in RADII:
            subset=[x for x in species_rows if x["scenario"]==scenario and x["radius_km"]==radius]
            q=[x for x in subset if x["status"].startswith("COMPUTABLE")]
            stable=[x for x in q if x["other_plant_groups_with_at_least_5_common_edges_n"]>=4]
            values=np.array([x["median_actual_minus_otherplant_same_edge_difference"] for x in q],float)
            result["by_scenario_radius"][f"{scenario}_{radius}km"]={
                "original_source_qualified_insects":len(subset),
                "species_with_at_least_one_comparable_otherplant_control":len(q),
                "species_with_at_least_four_groups_of_five_common_pairs":len(stable),
                "source_groups_contactable":sum(x["other_plant_groups_with_common_original_edges_n"] for x in subset),
                "negative_species_actual_host_loses_more":int((values<0).sum()),
                "positive_species_otherplant_loses_more":int((values>0).sum()),
                "equal_species":int((values==0).sum()),
                "median_species_group_median_actual_minus_otherplant_delta":float(np.median(values)) if len(values) else None,
                "median_modern_baseline_gap_actual_vs_otherplant":float(np.median(
                      [x["median_0BP_actual_vs_otherplant_baseline_gap"] for x in q])) if len(q) else None,
                "median_species_group_median_on_stable_5edge_controls":float(np.median(
                      [x["median_difference_on_groups_with_5plus_edges"] for x in stable])) if len(stable) else None}
    return groups,species_rows,result

def to_csv(path,rows):
    cols=sorted(set().union(*(x.keys() for x in rows)))
    with Path(path).open("w",newline="",encoding="utf8") as file:
        w=csv.DictWriter(file,fieldnames=cols,lineterminator="\n")
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    for name in ("model-plan","old-v06","localities","edges","primary-mask",
                 "extrap-mask","terrain","out"):
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    groups,species,summary=analysis(args)
    to_csv(args.out/"all_source_plant_group_fourway_contact_v08.csv",groups)
    to_csv(args.out/"all_292_insect_summary_v08.csv",species)
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()

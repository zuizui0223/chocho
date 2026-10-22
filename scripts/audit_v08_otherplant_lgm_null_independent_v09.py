#!/usr/bin/env python3
"""Independent, hash-locked audit of the completed v08 other-real-plant LGM control.

NO source model or genetic response is opened. Audits real-versus-other-plant
same-edge difference-in-differences, availability, modern baseline mismatch,
species/grain independence and the predeclared >=4 groups with >=5 edges
stability condition. This is exploratory, NOT a confirmatory equivalence test.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

SOURCES={
    "plan":"bb791eff085ba46ba6012efff0310a3a01d88144f4d0d18ef5b1160ff2c0aa10",
    "group":"2574cd7c6b346720733c232b2044e9c4d88407b0d597264b5ffcf115ef7cb261",
    "species":"d5e5e374e89eeac701c449e05edd775e7b0f0014c714b965a03378af0c7fa4b4",
    "result":"7fcd67aa367d138960aa7d0d94e6482139a5678a1c2a4be08a38dcbe4264bc55",
}
MODES=("analog","extrapolated")
RADII=(50,100,150)
TOL=1e-10

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as src:
        for buf in iter(lambda:src.read(1024*1024),b""):
            h.update(buf)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open("r",encoding="utf8",newline="") as src:
        return list(csv.DictReader(src))

def number(row,key,required=True):
    value=row.get(key)
    if value is None or value=="":
        if required:raise ValueError("missing numeric source column "+key)
        return None
    v=float(value)
    if not math.isfinite(v):raise ValueError("nonfinite numeric source column "+key)
    return v

def sign_counts(values):
    return {"actual_host_more_glacial_patch_loss":int(sum(x < 0 for x in values)),
            "other_plant_more_glacial_patch_loss":int(sum(x > 0 for x in values)),
            "equal_delta":int(sum(x == 0 for x in values))}

def same_edge_recompute(group):
    n=int(group["fourway_exact_common_original_edges_n"])
    if not 0<=n<=int(group["all_original_edges_n"]):
        raise ValueError("contact comparison cannot exceed original source edges")
    metrics=("actual_paired_delta","alternate_paired_delta","actual_minus_alternate_paired_delta",
             "modern_0BP_baseline_same_patch_gap")
    if n==0:
        if any(group.get(k) not in (None,"") for k in metrics):
            raise ValueError("ineligible donor group produced fake genetic-climate difference")
        return None
    get=lambda k:number(group,k)
    for prefix in ("actual","alternate"):
        current=get(prefix+"_0BP_same_fraction")
        historical=get(prefix+"_21ka_same_fraction")
        if not 0<=current<=1 or not 0<=historical<=1:
            raise ValueError("invalid source same-patch conditional fraction")
        split=get(prefix+"_same_to_different_n")
        merge=get(prefix+"_different_to_same_n")
        if abs(((merge-split)/n) - get(prefix+"_paired_delta"))>TOL:
            raise ValueError("paired source split/merge mismatch: "+prefix)
        if abs((historical-current)-get(prefix+"_paired_delta"))>TOL:
            raise ValueError("source net patch change is inconsistent with source 0BP/21ka fractions")
    contrast=get("actual_paired_delta")-get("alternate_paired_delta")
    if abs(contrast-get("actual_minus_alternate_paired_delta"))>TOL:
        raise ValueError("actual-other plant source difference-in-differences inconsistent")
    if abs(get("actual_0BP_same_fraction")-get("alternate_0BP_same_fraction")
           -get("modern_0BP_baseline_same_patch_gap"))>TOL:
        raise ValueError("original current climate patch baselines not compared correctly")
    return contrast

def summarize_prequalified_species(grouped,source_species,source_summary):
    expected={f"{mode}_{radius}km" for mode in MODES for radius in RADII}
    if set(source_summary["by_scenario_radius"])!=expected:
        raise ValueError("source scenario/radius identities have changed")
    metrics={}
    per_species={name:[] for name in source_species}
    for mode in MODES:
        for radius in RADII:
            groupkey=(mode,radius)
            eligible=[]
            stable=[]
            for species in sorted(source_species):
                controls=grouped.get((species,mode,radius),[])
                expected_groups=len(source_species[species]["alternative_groups"])
                if len(controls)!=expected_groups:
                    raise ValueError("frozen alternative group count mismatch: "+species)
                contrasts=[same_edge_recompute(g) for g in controls]
                valid=[x for x in contrasts if x is not None]
                fixed_row=source_species[species]["outcomes"][(mode,radius)]
                if int(fixed_row["source_frozen_other_plant_group_n"])!=len(controls):
                    raise ValueError("species denominator excludes source-frozen donor group")
                if int(fixed_row["other_plant_groups_with_common_original_edges_n"])!=len(valid):
                    raise ValueError("species source 4-way donor evaluability count drift")
                if valid:
                    d=median(valid)
                    if abs(d-number(fixed_row,"median_actual_minus_otherplant_same_edge_difference"))>TOL:
                        raise ValueError("v08 published species effect does not replay")
                    gaps=[number(g,"modern_0BP_baseline_same_patch_gap")
                          for g,c in zip(controls,contrasts) if c is not None]
                    if abs(median(gaps)-number(fixed_row,"median_0BP_actual_vs_otherplant_baseline_gap"))>TOL:
                        raise ValueError("v08 source species modern baseline gap does not replay")
                    eligible.append((species,d,median(gaps)))
                    stability=[c for g,c in zip(controls,contrasts)
                               if c is not None and int(g["fourway_exact_common_original_edges_n"])>=5]
                    if int(fixed_row["other_plant_groups_with_at_least_5_common_edges_n"])!=len(stability):
                        raise ValueError("source stable >=5 edges donor count mismatch")
                    if len(stability)>=4:
                        v=median(stability)
                        if abs(v-number(fixed_row,"median_difference_on_groups_with_5plus_edges"))>TOL:
                            raise ValueError("source stable group >=5 median source output did not replay")
                        stable.append((species,v))
                else:
                    if fixed_row["status"]!="NOT_EVALUABLE_FOURWAY_COMMON_EDGE_DENOMINATOR":
                        raise ValueError("missing all comparison donors not declared unevaluable")
                    if number(fixed_row,"median_actual_minus_otherplant_same_edge_difference",False) is not None:
                        raise ValueError("unevaluable species returned a numeric effect")
            control_group_count=sum(len(grouped.get((s,mode,radius),[])) for s in source_species)
            source=source_summary["by_scenario_radius"][f"{mode}_{radius}km"]
            vals=[x[1] for x in eligible]
            signs=sign_counts(vals)
            if (len(eligible)!=source["species_with_at_least_one_comparable_otherplant_control"] or
                signs["actual_host_more_glacial_patch_loss"]!=source["negative_species_actual_host_loses_more"] or
                signs["other_plant_more_glacial_patch_loss"]!=source["positive_species_otherplant_loses_more"] or
                signs["equal_delta"]!=source["equal_species"] or
                len(stable)!=source["species_with_at_least_four_groups_of_five_common_pairs"] or
                control_group_count!=2108):
                raise ValueError("original v08 summary and independent recount differ")
            baseline=[abs(x[2]) for x in eligible]
            metrics[f"{mode}_{radius}km"]={
                "original_prequalified_species":len(source_species),
                "source_frozen_alternative_plant_groups":control_group_count,
                "species_at_least_one_fourway_pair":len(eligible),
                "unevaluable_zero_fourway_species":len(source_species)-len(eligible),
                "species_eligible_genus_count":len({x[0].split()[0] for x in eligible}),
                "species_median_actual_minus_otherplant":median(vals) if vals else None,
                "signs":signs,
                "species_with_four_donors_each_at_least_five_common_edges":len(stable),
                "stable_signs":sign_counts([x[1] for x in stable]),
                "stable_species_median":median([x[1] for x in stable]) if stable else None,
                "median_abs_modern_baseline_gap":median(baseline) if baseline else None,
                "modern_baseline_gap_exact_zero_species":int(sum(x==0 for x in baseline)),
                "modern_baseline_abs_gap_over_0p1_species":int(sum(x>.1 for x in baseline)),
                "median_source_group_matchable_edges":float(median([
                    number(source_species[s]["outcomes"][(mode,radius)],
                        "median_fourway_same_edge_n_per_group") for s,_,_ in eligible])) if eligible else None
            }
            if metrics[f"{mode}_{radius}km"]["species_median_actual_minus_otherplant"]!=source["median_species_group_median_actual_minus_otherplant_delta"]:
                raise ValueError("original v08 source median not reproduced")
    return metrics

def audit(plan_file,group_file,species_file,result_file):
    for key,path in (("plan",plan_file),("group",group_file),("species",species_file),("result",result_file)):
        if sha(path)!=SOURCES[key]:
            raise ValueError("immutable original result/source SHA256 mismatch: "+key)
    plan=json.loads(Path(plan_file).read_text())
    source=json.loads(Path(result_file).read_text())
    if (plan["status"]!="ELIGIBILITY_ONLY_NO_NEW_LGM_PATCH_OUTCOME_READ" or
        not plan["no_new_21ka_climate_or_genetic_patch_response_accessed"] or
        len(plan["matched_species"])!=292 or
        plan["valid_frozen_alternative_group_count"]!=2108 or
        source["status"]!="EXPLORATORY_HOST_SPECIFIC_PATCH_CLIMATE_NEGATIVE_CONTROL_NOT_GENETIC_OR_OCCUPANCY"):
        raise ValueError("original early frozen donor eligibility or scientific scope changed")
    group_rows=read_csv(group_file)
    species_rows=read_csv(species_file)
    if len(group_rows)!=12648 or len(species_rows)!=1752:
        raise ValueError("not all 292 insects and 2108 frozen donors in both climate variants and radii")
    species={sp:{"alternative_groups":detail["alternative_groups"],"outcomes":{}}
             for sp,detail in plan["matched_species"].items()}
    plan_source=set(species)
    for row in species_rows:
        name=row["insect_species"];scenario=row["scenario"];radius=int(row["radius_km"])
        if name not in species or scenario not in MODES or radius not in RADII:
            raise ValueError("new/ineligible modelled insect or scenario")
        key=(scenario,radius)
        if key in species[name]["outcomes"]:
            raise ValueError("duplicated species scenario radius")
        species[name]["outcomes"][key]=row
    if any(len(x["outcomes"])!=6 for x in species.values()):
        raise ValueError("some original source insect scenario/radius rows missing")
    groups=defaultdict(list)
    for row in group_rows:
        insect=row["insect_species"];scenario=row["scenario"];radius=int(row["radius_km"])
        if insect not in plan_source or scenario not in MODES or radius not in RADII:
            raise ValueError("unspecified other plant donor source row")
        g=int(row["other_plant_group_index"])
        frozen={int(x["group_index"]):x for x in plan["matched_species"][insect]["alternative_groups"]}
        if g not in frozen:
            raise ValueError("group added after source 0BP donor freeze")
        if row["other_plant_ids"]!=";".join(frozen[g]["donor_taxa"]):
            raise ValueError("donor plant identity mismatches original 0BP-only selection")
        if row["known_actual_host_ids"]!=";".join(plan["matched_species"][insect]["known_host_taxa"]):
            raise ValueError("actual accepted host identity changed from plan")
        if abs(float(row["original_0BP_climate_moment_gap"])-
               float(frozen[g]["modern_climate_moment_gap"]))>TOL:
            raise ValueError("original source modern climate moment mismatch was retuned")
        groups[(insect,scenario,radius)].append(row)
    if len(groups)!=292*2*3:
        raise ValueError("not every insect/scenario/radius donor group represented")
    for key,rows in groups.items():
        rows.sort(key=lambda v:int(v["other_plant_group_index"]))
        if len({int(r["other_plant_group_index"]) for r in rows})!=len(rows):
            raise ValueError("duplicate source-plant group result")
    stats=summarize_prequalified_species(groups,species,source)
    primary=stats["analog_100km"]
    # Qualified 292 source insects had source-external plant model availability; 
    # 66 could not support even one complete 4-way same-edge contrast.
    assert (primary["species_at_least_one_fourway_pair"],primary["unevaluable_zero_fourway_species"])==(226,66)
    assert primary["signs"]=={"actual_host_more_glacial_patch_loss":62,
                              "other_plant_more_glacial_patch_loss":79,"equal_delta":85}
    assert primary["stable_signs"]=={"actual_host_more_glacial_patch_loss":56,
                                     "other_plant_more_glacial_patch_loss":67,"equal_delta":53}
    return {"schema":"chocho_v09_independent_source_locked_v08_other_real_plant_null_audit",
       "status":"EXPLORATORY_NULL_FAILS_TO_SUPPORT_HOST_IDENTITY_SPECIFIC_LGM_PATCH_SPLITTING",
       "source_all_insects":641,"source_v05_original_0BP_matched_insects":335,
       "source_qualifying_insects":292,
       "source_0BP_preselected_other_plant_groups":2108,
       "source_group_scenario_radius_rows":12648,"source_species_scenario_radius_rows":1752,
       "sources_sha256":dict(SOURCES),"by_scenario_radius":stats,
       "claim":"Compared with fixed modern-climate-matched other REAL plant identities, recorded larval hosts are not generally more vulnerable to hypothetical 21ka climatic patch splitting at modern insect sampling coordinates.",
       "meaning_constraints":[
            "Median zero does not prove equivalence; many model-eligible species are unevaluable and contact is conditioned on surviving four-way original edges",
            "An unlisted other plant is not a biologically confirmed non-host and matching 0BP climate moments does not guarantee equivalent modern 0BP patch baselines",
            "Repeated 2108 groups are not independent insects; separate insect models reuse shared host plants and source geographies",
            "This is potential climate-land spatial adjacency, not actual LGM occupied host-plant populations, insect dispersal or genetic differentiation",
            "TTF original v0.5 other-real-plant negative climate identity test remains unchanged; no genetic response reopened"]}

def main():
    parser=argparse.ArgumentParser()
    for name in ("plan","group","species","result","out"):
        parser.add_argument("--"+name,required=True,type=Path)
    a=parser.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    report=audit(a.plan,a.group,a.species,a.result)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":report["status"],
        "primary_100km":report["by_scenario_radius"]["analog_100km"],
        "extrapolated_100km":report["by_scenario_radius"]["extrapolated_100km"]},sort_keys=True))
if __name__=="__main__":main()

#!/usr/bin/env python3
"""31-species archived butterfly-indexing screen of port×wingspan association.

This is a *conditional observational index pattern*, not a test of true
establishment, ship transport, host utilization, or the full 239-species panel.
Under 13 observed butterfly×island cells, the presence+interaction model is
necessarily sparse; report exactly what the permutation can/cannot identify.
"""
from __future__ import annotations
import argparse,csv,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path
SEED=20261102
REPS=9999

def rows(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def save(path,rr):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def quantile(vals,q):
    a=sorted(vals)
    return a[int((len(a)-1)*q)] if a else None
def main():
    ap=argparse.ArgumentParser()
    for k in ("q2-edges","q2-audit","q3-observed","q3-all-observed",
              "q3-audit","panel","ports","area","out-dir"):
        ap.add_argument("--"+k,type=Path,required=True)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q2=json.loads(a.q2_audit.read_text(encoding="utf8"))
    q3=json.loads(a.q3_audit.read_text(encoding="utf8"))
    if q2["scope"]["added_strict_island_cells"]!=5741 or q3["validation_panel"]["recorded_butterfly_in_new_opportunity_cells"]!=13:
        raise RuntimeError("Not corrected Q2/Q3 source population")
    edges=rows(a.q2_edges)
    if len(edges)!=11573:raise RuntimeError("Corrected Q2 source edges missing")
    val=rows(a.panel)
    use={r["species"]:r for r in val if int(r["occurrence_records"])>0}
    if len(use)!=31:raise RuntimeError("Wrong original 31 GBIF validation species")
    pos=rows(a.q3_observed)
    posset={(r["species"],r["area_code_l3"]) for r in pos}
    if len(posset)!=13 or len(pos)!=13:raise RuntimeError("Q3 13 records drift")
    full=rows(a.q3_all_observed)
    if len(full)!=159:raise RuntimeError("Source 31-butterfly × island observation matrix drift")
    island_record_totals=Counter()
    focal_record_counts={}
    for r in full:
        count=int(r["observed_records"])
        key=(r["species"],r["area_code_l3"])
        focal_record_counts[key]=count
        island_record_totals[r["area_code_l3"]]+=count
    source={(r["species"],r["island_code"]) for r in edges if r["species"] in use}
    if len(source)!=q3["validation_panel"]["sampled_species_added_island_opportunity_cells"] or len(source)!=920:
        raise RuntimeError("Q2/Q3 denominator mismatch for 31 species")
    if not posset<=source:
        raise RuntimeError("Observed positive outside corrected host-resource opportunities")
    portrows=rows(a.ports)
    port={r["area_code_l3"]:r for r in portrows}
    if len(port)!=369 or sum(r["geographic_group"]=="island" for r in portrows)!=126:
        raise RuntimeError("Not corrected 126-island WPI input")
    arearows=rows(a.area)
    areamap={r["code"]:r for r in arearows}
    if len(areamap)<120:raise RuntimeError("Area join missing")
    spans={sp:float(use[sp]["wing_size_proxy"]) for sp in use}
    if any(not math.isfinite(x) for x in spans.values()):raise RuntimeError("Invalid wing span")
    wingmean=statistics.mean(spans.values())
    wingsd=statistics.pstdev(spans.values())
    if wingsd<=0:raise RuntimeError("Wing variance zero")
    wingz={sp:(x-wingmean)/wingsd for sp,x in spans.items()}
    df=[]
    for sp,isle in sorted(source):
        p=port[isle];ar=areamap[isle]
        if p["geographic_group"]!="island":raise RuntimeError("Non island opportunity in source")
        nport=int(p["WPI_port_count"])
        effort_other=max(0,island_record_totals[isle]-focal_record_counts.get((sp,isle),0))
        area=float(ar["area_km2"])
        if area<=0:raise RuntimeError("Invalid area")
        df.append({
          "species":sp,"wgsrpd3_code":isle,"level1":p["level1"],
          "archival_butterfly_record_at_added_resource":int((sp,isle) in posset),
          "wing_size_proxy_original":spans[sp],
          "wing_z":wingz[sp],"indexed_WPI_port_count_2023":nport,
          "log1p_indexed_ports":math.log1p(nport),
          "WGSRPD3_area_km2":area,"log_area":math.log(area),
          "other_30_butterfly_record_count":effort_other,
          "log1p_other_butterfly_records":math.log1p(effort_other),
          "had_other_butterfly_observation":int(effort_other>0),
          "sampled_occurrence_count_for_species":int(use[sp]["occurrence_records"]),
          "not_recorded_in_this_DOI_not_true_absence":1
        })
    if len(df)!=920 or sum(r["archival_butterfly_record_at_added_resource"] for r in df)!=13:
        raise RuntimeError("Indexing source population inconsistent")
    # A standardized across-species wing trait cannot itself be exchanged within
    # species. The interaction score tests whether port association differs
    # across species after holding all per-species index totals fixed.
    # Observed means over all recorded cells, not probability of colonization.
    features={
      "port_index_only":[r["log1p_indexed_ports"] for r in df],
      "wing_x_port":[r["wing_z"]*r["log1p_indexed_ports"] for r in df],
      "log_region_area":[r["log_area"] for r in df],
      "other_species_GBIF_index_effort":[r["log1p_other_butterfly_records"] for r in df],
    }
    positives=[i for i,r in enumerate(df) if r["archival_butterfly_record_at_added_resource"]]
    n_event=len(positives)
    def test(groups):
        blocks=defaultdict(list)
        for i,label in enumerate(groups):blocks[label].append(i)
        counts={k:sum(df[i]["archival_butterfly_record_at_added_resource"] for i in ix)
                for k,ix in blocks.items()}
        movable=[indices for k,indices in blocks.items() if 0<counts[k]<len(indices)]
        movable_positive=sum(counts[k] for k,ix in blocks.items() if 0<counts[k]<len(ix))
        obs={name:sum(values[i] for i in positives)/n_event for name,values in features.items()}
        rng=random.Random(SEED)
        null={name:[] for name in features}
        frozen=[i for k,indices in blocks.items() if counts[k]==len(indices) for i in indices]
        for _ in range(REPS):
            sel=frozen.copy()
            for k,indices in blocks.items():
                if 0<counts[k]<len(indices):
                    sel.extend(rng.sample(indices,counts[k]))
            if len(sel)!=n_event:raise RuntimeError("Index event marginal drift")
            for name,values in features.items():
                null[name].append(sum(values[i] for i in sel)/n_event)
        stats={}
        for name,obsvalue in obs.items():
            d=null[name]
            mean=statistics.mean(d)
            stats[name]={
                "observed_mean_at_indexed_resource_cells":obsvalue,
                "conditional_randomization_mean":mean,
                "observed_minus_null_mean":obsvalue-mean,
                "null_q025":quantile(d,.025),"null_q975":quantile(d,.975),
                "p_high":(1+sum(x>=obsvalue-1e-12 for x in d))/(1+len(d)),
                "p_low":(1+sum(x<=obsvalue+1e-12 for x in d))/(1+len(d)),
                "p_two_centered":(1+sum(abs(x-mean)>=abs(obsvalue-mean)-1e-12 for x in d))/(1+len(d)),
                "note":"Difference in event-weighted predictor score, NOT a regression coefficient or colonization effect"
            }
        return {
          "n_butterfly_resource_index_events":n_event,
          "n_blocks":len(blocks),
          "exchangeable_blocks":len(movable),
          "exchangeable_index_events":movable_positive,
          "permutations":REPS,
          "scores":stats,
        }
    area_cut=statistics.median(r["log_area"] for r in df)
    blocks={
      "within_species":[r["species"] for r in df],
      "within_species_level1":[r["species"]+"|"+r["level1"] for r in df],
      "within_species_level1_and_area_half":[r["species"]+"|"+r["level1"]+"|"+
               ("small" if r["log_area"]<=area_cut else "large") for r in df],
      "within_species_level1_and_other_butterfly_records":[r["species"]+"|"+
               r["level1"]+"|"+("recorded" if r["had_other_butterfly_observation"] else "unrecorded") for r in df],
    }
    stats={name:test(block) for name,block in blocks.items()}
    outcomes=[r for r in df if r["archival_butterfly_record_at_added_resource"]]
    observed_by_species=Counter(x["species"] for x in outcomes)
    obs_wing=[x["wing_size_proxy_original"] for x in outcomes]
    report={
       "schema":"chocho_Q3_31_species_wing_WPI_indexing_screen_v0.1",
       "status":"SPARSE_GBIF_INDEXING_PATTERNS_NOT_ESTABLISHMENT_OR_ARRIVAL",
       "source":{
         "GEB_species":239,"Q2_corrected_added_island_cells":5741,
         "observed_DOI":"10.15468/dl.pp5nc9",
         "source_31_observed_species":31,
         "31_species_modern_host_resource_opportunities":len(df),
         "source_observed_species_x_resource_island_cells":n_event,
         "WPI_modern_published_port_snapshot_year":2023,
         "wing_size_proxy":"GEB frozen LepTraits continuous wing-size proxy",
       },
       "observed":{
         "observed_index_events":n_event,
         "species_with_any_recorded_added_opportunity":len(observed_by_species),
         "index_events_at_0_ports":sum(r["indexed_WPI_port_count_2023"]==0 for r in outcomes),
         "index_events_at_1plus_ports":sum(r["indexed_WPI_port_count_2023"]>=1 for r in outcomes),
         "all_opportunities_at_0_ports":sum(r["indexed_WPI_port_count_2023"]==0 for r in df),
         "all_opportunities_at_1plus_ports":sum(r["indexed_WPI_port_count_2023"]>=1 for r in df),
         "wing_median_in_detected_cells":statistics.median(obs_wing) if obs_wing else None,
         "wing_median_in_all_opportunity_cells":statistics.median(r["wing_size_proxy_original"] for r in df),
       },
       "conditional_index_permutations":stats,
       "stop_gate":{
         "min_positive_cells_to_fit_four_covariates_and_3_interactions":40,
         "sufficient_independent_species_or_cells":False,
         "reason":"13 observed index cells distributed over too few independent species, nonrandom archive and unequal island reporting; effect coefficient would be unstable.",
         "no_multivariable_port_x_wing_establishment_model_fitted":True,
       },
       "limitations":[
          "The 31 sampled species from a frozen historical GBIF validation archive are not the target 239- or 224-species full island-year risk set.",
          "GBIF index non-record is not absence; tests concern where a nonrandom subset was OBSERVED, not where it colonized.",
          "Modern WPI 2023 index does not measure shipping volume, transport modes, or port status at any first butterfly arrival.",
          "Wing size is not direct flight capacity and may itself influence detectability; conditioning on each species index count only partially addresses bias.",
          "Other-butterfly archived records indicate observation activity but also ecological butterfly richness and taxonomic geography, not independent effort.",
          "Within-level1 permutations preserve source species and coarse biogeography but not local climate, land cover or observer preferences.",
          "13 observations cannot identify a wing×port colonization mechanism or a ship-vs-flight distinction.",
       ],
    }
    save(a.out_dir/"fixed_31_species_added_resource_island_index_status.csv",df)
    (a.out_dir/"wing_port_observation_screen_audit.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({key:val for key,val in report.items() if key!="conditional_index_permutations"},indent=2,ensure_ascii=False))
    print("WING_X_PORT_NULLS",json.dumps({k:v["scores"]["wing_x_port"] for k,v in stats.items()},ensure_ascii=False))
if __name__=="__main__":
    main()

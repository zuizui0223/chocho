#!/usr/bin/env python3
"""Offline audit/recovery of first-run paired GBIF island port/nonport site indices.
Never refetches GBIF or fills missing observations with ecological absences.
Original CI run 42601343703 computed all 48 sites but failed on an int64 JSON
serialization at final write; source aggregate CSVs are in its artifact.
"""
from __future__ import annotations
import csv,hashlib,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import rankdata,binomtest

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"results/island_port_nonport_lag_20261109"
RUN_ID=42601343703
ARTIFACT_ID=13666735338

def get(name):
    path=DIR/name
    with path.open(encoding="utf8",newline="") as f:rows=list(csv.DictReader(f))
    return rows,hashlib.sha256(path.read_bytes()).hexdigest()

def num(d,k):
    return float(d[k])
def mean(xs):
    return float(np.mean(xs)) if len(xs) else None
def interval(xs):
    v=np.quantile(xs,[.025,.975])
    return [float(v[0]),float(v[1])]

def sp(x,y,draws=10000):
    if len(x)<6 or len(set(x))<2 or len(set(y))<2:
        return {"n_islands":len(x),"rho":None,"p_island_shuffle":None,"ci_bootstrap":None}
    ax=rankdata(x).astype(float);ay=rankdata(y).astype(float)
    ax-=ax.mean();ay-=ay.mean()
    norm=np.linalg.norm(ax)*np.linalg.norm(ay)
    obs=float(ax@ay/norm)
    rng=np.random.default_rng(20261109)
    more=0
    for _ in range(draws):
        more+=bool(abs(float(ax@rng.permutation(ay)/norm))>=abs(obs)-1e-12)
    boot=[]
    rng=np.random.default_rng(20261110)
    n=len(x);x=np.asarray(x);y=np.asarray(y)
    for _ in range(3000):
        q=rng.integers(0,n,n)
        if len(set(x[q]))<2 or len(set(y[q]))<2:continue
        xx=rankdata(x[q]); yy=rankdata(y[q])
        xx=xx-xx.mean();yy=yy-yy.mean()
        boot.append(float(xx@yy/(np.linalg.norm(xx)*np.linalg.norm(yy))))
    return {"n_islands":n,"rho":obs,"p_island_shuffle":(1+more)/(1+draws),
      "ci_bootstrap":interval(boot) if boot else None,"method":"island-level rank-permutation; island bootstrap"}

def describe(x):
    x=np.asarray(x,dtype=float)
    rng=np.random.default_rng(20261112)
    boot=[float(np.mean(x[rng.integers(0,len(x),len(x))])) for _ in range(5000)]
    return {"n_islands":len(x),"island_equal_mean":float(x.mean()),
       "island_median":float(np.median(x)),
       "island_bootstrap_95pct_mean":interval(boot),
       "positive_islands":int(sum(x>0)),"negative_islands":int(sum(x<0)),
       "tied_islands":int(sum(x==0)),"paired_sign_test_p":float(binomtest(int(sum(x>0)),int(sum(x!=0)),.5).pvalue) if int(sum(x!=0)) else None}

def main():
    sites,sitehash=get("port_nonport_plant_butterfly_site_indices.csv")
    pairs,pairhash=get("matched_pair_host_lag_differences.csv")
    islands,islandhash=get("island_equal_weight_host_lag_differences.csv")
    focal,focalhash=get("focal_named_butterfly_host_time_order.csv")
    fixed,fixedhash=get("fixed_same_island_port_nonport_sites.csv")
    assert len(sites)==48 and len(pairs)==24 and len(fixed)==24
    assert len(islands)==17 and len(set(x["landmass_id"] for x in islands))==17
    assert len(focal)==len(sites)*3==144
    side={(s["landmass_id"],s["port_id"],s["side"]):s for s in sites}
    assert len(side)==len(sites)
    grouped=defaultdict(list)
    descriptive_pairs=[]
    for p in pairs:
        key=(p["landmass_id"],p["port_id"])
        a=side[key+("port",)]
        b=side[key+("nonport",)]
        hpa=num(a,"gbif_host_candidate_pre_specieskeys")
        hco=num(b,"gbif_host_candidate_pre_specieskeys")
        apa=num(a,"gbif_plant_pre_specieskeys")
        aco=num(b,"gbif_plant_pre_specieskeys")
        rpa=num(a,"gbif_plant_pre_total_records")
        rco=num(b,"gbif_plant_pre_total_records")
        assert hpa<=apa and hco<=aco
        # Proportions: candidate Lepidoptera hosts among all indexed
        # plant species, not confirmed butterfly-specific food plants.
        shpa=hpa/apa if apa else 0.
        shco=hco/aco if aco else 0.
        hrel=np.log1p(hpa)-np.log1p(hco)
        arel=np.log1p(apa)-np.log1p(aco)
        early=np.log1p(num(a,"gbif_butterfly_new_2011_2017_225_species"))-np.log1p(num(b,"gbif_butterfly_new_2011_2017_225_species"))
        late=np.log1p(num(a,"gbif_butterfly_new_2018_2022_225_species"))-np.log1p(num(b,"gbif_butterfly_new_2018_2022_225_species"))
        old=np.log1p(num(a,"gbif_butterfly_pre_225_species"))-np.log1p(num(b,"gbif_butterfly_pre_225_species"))
        grouped[key[0]].append({
         "port_id":key[1],
         "host_delta_log1p":hrel,
         "all_plant_delta_log1p":arel,
         "record_effort_delta_log1p":np.log1p(rpa)-np.log1p(rco),
         "relative_host_to_all_plant_log_delta":hrel-arel,
         "host_fraction_port_minus_nonport":shpa-shco,
         "early_new_butterfly_delta_log1p":early,
         "late_new_butterfly_delta_log1p":late,
         "prior_butterfly_delta_log1p":old,
         "botanical_both_sides_10records":int(rpa>=10 and rco>=10),
         "host_port":hpa,"host_nonport":hco,"plant_port":apa,"plant_nonport":aco,
         "butterfly_early_port":num(a,"gbif_butterfly_new_2011_2017_225_species"),
         "butterfly_early_nonport":num(b,"gbif_butterfly_new_2011_2017_225_species"),
         "site_pre2011_plant_records_port":rpa,
         "site_pre2011_plant_records_nonport":rco
        })
        descriptive_pairs.append((a,b))
    fields=list(next(iter(grouped.values()))[0].keys())
    fields=[k for k in fields if k!="port_id" and k!="botanical_both_sides_10records"]
    rows=[]
    for lid,ps in sorted(grouped.items(),key=lambda z:int(z[0])):
        d={"landmass_id":lid,"paired_ports":len(ps),"botanical_pairs_both_sides_10records":sum(v["botanical_both_sides_10records"] for v in ps)}
        d.update({k:mean([float(v[k]) for v in ps]) for k in fields})
        rows.append(d)
    assert len(rows)==17
    h=[r["host_delta_log1p"] for r in rows]
    hrel=[r["relative_host_to_all_plant_log_delta"] for r in rows]
    hshare=[r["host_fraction_port_minus_nonport"] for r in rows]
    early=[r["early_new_butterfly_delta_log1p"] for r in rows]
    late=[r["late_new_butterfly_delta_log1p"] for r in rows]
    plant=[r["all_plant_delta_log1p"] for r in rows]
    effort=[r["record_effort_delta_log1p"] for r in rows]
    prior=[r["prior_butterfly_delta_log1p"] for r in rows]
    robust=[r for r in rows if r["botanical_pairs_both_sides_10records"]>0]
    fk=[r for r in focal if r["compatible_plant_recorded_2001_2010"].lower()=="true"]
    def first_after(f):
        return (f["butterfly_indexed_pre2001_2010"].lower()!="true" and
                (f["butterfly_indexed_2011_2017"].lower()=="true" or f["butterfly_indexed_2018_2022"].lower()=="true"))
    source={
      "schema":"chocho_island_port_nonport_20261109_source_artifact_offline_audit_v1",
      "status":"RECOVERED_ALL_PAIRS_FROM_PREVIOUS_COMPLETE_48SITE_GBIF_RUN",
      "gbif_source_workflow_run":RUN_ID,"original_artifact_id":ARTIFACT_ID,
      "source_csv_sha256":{"sites":sitehash,"matched_pairs":pairhash,
                           "islands":islandhash,"focal":focalhash,"fixed_sites":fixedhash},
      "cohort":{"matched_islands":len(rows),"matched_port_nonport_pairs":len(pairs),
               "individual_GBIF_site_indexes":len(sites),
               "focal_named_pair_site_rows":len(focal),
               "paired_ports_with_both_sides_ge10_plant_records":sum(p["botanical_both_sides_10records"] for ps in grouped.values() for p in ps),
               "islands_with_ge1_botanical_both_sides_ge10":len(robust)},
      "host_index_advantage":{"raw_candidate_species_log1p":describe(h),
            "candidate_fraction_among_all_indexed_plants":describe(hshare),
            "candidate_minus_all_plant_log1p":describe(hrel),
            "all_plant_richness_log1p":describe(plant),
            "plant_occurrence_record_count_log1p":describe(effort)},
      "subsequent_butterfly_record_advantage":{"new_2011_2017_log1p":describe(early),
                                               "new_2018_2022_log1p":describe(late)},
      "associations": {
        "raw_host_delta_vs_new_2011_2017":sp(h,early),
        "relative_host_enrichment_vs_new_2011_2017":sp(hrel,early),
        "candidate_fraction_delta_vs_new_2011_2017":sp(hshare,early),
        "relative_host_enrichment_vs_new_2018_2022":sp(hrel,late),
        "all_plant_delta_vs_new_2011_2017_negative_control":sp(plant,early),
        "all_plant_effort_delta_vs_new_2011_2017_negative_control":sp(effort,early),
        "prior_butterfly_delta_vs_new_2011_2017_negative_control":sp(prior,early)},
      "both_sides_ge10_botanical_records_sensitivity":{
          "islands":len(robust),
          "relative_host_enrichment_vs_new_2011_2017":sp(
                [r["relative_host_to_all_plant_log_delta"] for r in robust],
                [r["early_new_butterfly_delta_log1p"] for r in robust])},
      "focal_species_identity_check":{
           "fixed_source_links":["Colias croceus / Medicago sativa","Colias hyale / Medicago sativa","Nymphalis antiopa / Alnus incana"],
           "compatible_plant_pre2011_indexed_site_rows":len(fk),
           "earlier_host_and_new_butterfly_coindexed_site_rows":sum(first_after(r) for r in fk),
           "earlier_host_and_new_butterfly_coindexed_port_rows":sum(first_after(r) and r["side"]=="port" for r in fk),
           "earlier_host_and_new_butterfly_coindexed_nonport_rows":sum(first_after(r) and r["side"]=="nonport" for r in fk)
      },
      "limitations":["Selection conditioned on >=30 earlier GBIF plant records at PORT, so raw port host richness is biased upward by construction",
        "No independent native-versus-alien status of focal hosts on each island in the historic baseline",
        "687 potential host species catalogue is global Lepidoptera, not a 225-butterfly-by-plant interaction network",
        "This is a new 15-125km remote coastline matching, NOT old ESA 2021 remote coast locations",
        "Port and remote coast climate, built-up area, survey effort and habitat were not jointly matched",
        "Index is neither true absence nor first arrival, butterfly establishment, reproduction or actual host consumption",
        "Only 17 islands and 24 pairs are included; overlapping 5km discs, different early GBIF sampling",
        "Earlier flora and butterfly periods were queried on 2026 index; GBIF revisions can change historical counts",
        "Island-equal means and permutation treat island as sampling unit, not 24 independent port replicates"]
    }
    with (DIR/"recovered_island_equal_weight_host_lag.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (DIR/"audit.json").write_text(json.dumps(source,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("RECOVERED_ISLAND_PORT_NONPORT_FINAL",json.dumps({
       "cohort":source["cohort"],"host_advantage":source["host_index_advantage"],
       "associations":source["associations"],
       "focal":source["focal_species_identity_check"]},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()

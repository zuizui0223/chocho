#!/usr/bin/env python3
"""Q4 observational-age sampling-bias diagnostic for introduced host/island GBIF indices.

The earliest GBIF plant year is a *dated documentation index*. Older first
records can mechanically be more frequent in well surveyed regions. We
quantify this association before treating record-derived plant residence
length as a butterfly-colonization lag predictor.

Frozen 500 response-blind host-island source; 219 positive dated plant pairs;
one WGSRPD3 regional outcome equal to the MEDIAN earliest plant indexed year.
No plant arrival date, no causal inference, no original butterfly index data.
"""
from __future__ import annotations
import argparse,csv,json,math,random,statistics
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np

SEED=20261102; B=4999

def load(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def num(x):
    try:
        y=float(x);return y if math.isfinite(y) else None
    except (ValueError,TypeError):return None
def spearman(x,y):
    def ranks(a):
        pairs=sorted(enumerate(a),key=lambda k:k[1])
        rr=[0.0]*len(a)
        k=0
        while k<len(pairs):
            h=k+1
            while h<len(pairs) and pairs[h][1]==pairs[k][1]:h+=1
            avg=(k+h-1)/2.+1.
            for t in range(k,h):rr[pairs[t][0]]=avg
            k=h
        return np.array(rr)
    xx=ranks(x);yy=ranks(y)
    return float(np.corrcoef(xx,yy)[0,1])

def fwl(y,x,z):
    q,_=np.linalg.qr(z,mode="reduced")
    yr=y-q@(q.T@y);xr=x-q@(q.T@x)
    den=float(xr@xr)
    return (None if den<=1e-9 else float((xr@yr)/den)),q,yr

def fit(rows,correction="total"):
    n=len(rows)
    levels=sorted({r["level1"] for r in rows})
    if n<45 or len(levels)<6:raise RuntimeError(f"Not enough island effort records for independent region analysis: n={n}, groups={levels}")
    y=np.array([r["median_earliest_year"] for r in rows])
    x=np.array([r["ln_effort"] if correction=="total" else r["ln_pre1900_effort"] for r in rows])
    Z=np.column_stack([
       np.ones(n),
       np.array([r["ln_area"] for r in rows]),
       np.array([r["ln_dated_species_count"] for r in rows]),
       np.array([r["ln_introduced_flora"] for r in rows]),
       *[np.array([1 if r["level1"]==lev else 0 for r in rows]) for lev in levels[1:]]
    ])
    if np.linalg.matrix_rank(Z)!=Z.shape[1]:
        raise RuntimeError("Singular detection-effort model")
    beta,q,yr=fwl(y,x,Z)
    if beta is None:raise RuntimeError("Insufficient residual botanical effort")
    strata=defaultdict(list);bylvl=defaultdict(list)
    for i,r in enumerate(rows):bylvl[r["level1"]].append(i)
    for lev,idx in bylvl.items():
        median=float(np.median([rows[i]["ln_area"] for i in idx]))
        for i in idx:
            stratum=lev+"|"+("small" if rows[i]["ln_area"]<=median else "large")
            strata[stratum].append(i)
    movable=[ids for ids in strata.values() if len(ids)>=2 and len({float(x[i]) for i in ids})>1]
    rng=random.Random(SEED)
    null=[]
    for _ in range(B):
        xx=x.copy()
        for ids in movable:
            shuffled=ids.copy()
            rng.shuffle(shuffled)
            xx[ids]=x[shuffled]
        xr=xx-q@(q.T@xx);den=float(xr@xr)
        if den>1e-9:null.append(float(xr@yr/den))
    if len(null)<4900:raise RuntimeError("Too few valid exposure-shuffles")
    vals=np.array(null);center=float(np.mean(vals))
    return {
      "n_independent_WGSRPD3_regions":n,
      "response":"region median of first indexed year among fixed 500 originally sampled host-island plants with dated GBIF records",
      "exposure":correction,
      "controls":["ln_WGSRPD3_area","ln1p_number_dated_host_taxa_in_original_pilot","ln1p_introduced_WCVP_flora","Level1 fixed effects"],
      "regression_year_change_per_ln1p_effort_unit":beta,
      "descriptive_raw_spearman":spearman([r["ln_effort" if correction=="total" else "ln_pre1900_effort"] for r in rows],[r["median_earliest_year"] for r in rows]),
      "permutation":{
         "shuffle":"botanical effort labels within Level1 x within-Level1 area half",
         "draws":len(vals),"exchangeable_regions":len({i for g in movable for i in g}),
         "exchangeable_strata":len(movable),
         "null_mean":center,"null_q025":float(np.quantile(vals,.025)),
         "null_q975":float(np.quantile(vals,.975)),
         "observed_minus_null":beta-center,
         "two_sided_centered_p":(1+int(np.sum(np.abs(vals-center)>=abs(beta-center)-1e-12)))/(1+len(vals))
      }
    }
def main():
    p=argparse.ArgumentParser()
    for k in ("plant-dates","original-audit","botanical-effort","all-flora","area","out-dir"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    dates=load(a.plant_dates)
    original=json.loads(a.original_audit.read_text(encoding="utf-8"))
    effort={r["island_code"]:r for r in load(a.botanical_effort)}
    flora={r["code"]:r for r in load(a.all_flora)}
    area={r["code"]:r for r in load(a.area)}
    if len(dates)!=500 or original["pilot_dated_pairs"]!=219 or len(flora)!=126 or len(effort)!=122:
        raise RuntimeError("Frozen index source drift")
    grouped=defaultdict(list)
    for r in dates:
        if r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND":
            y=num(r["earliest_gbif_indexed_year"])
            if y is None:raise RuntimeError("Dated host row lacks indexed first year")
            grouped[r["island_code"]].append(y)
    if sum(len(z) for z in grouped.values())!=219 or len(grouped)!=82:
        raise RuntimeError(f"Unexpected dated plant source: n={sum(map(len,grouped.values()))}, regions={len(grouped)}")
    cohort=[]
    status=Counter()
    for code,years in sorted(grouped.items()):
        if code not in effort or code not in flora or code not in area:raise RuntimeError("Missing pinned island data "+code)
        f=effort[code]
        if f["gbif_all_plants_index_status"]!="DATED_OCCURRENCE_FOUND":
            status["botanical_effort_missing"]+=1
            continue
        total=num(f["gbif_all_plants_total_records"])
        pre=num(f["gbif_all_plants_pre1900_count"])
        flor=num(flora[code]["introduced"])
        if None in (total,pre,flor) or total<=0 or flor<0:
            status["invalid_effort_counts"]+=1;continue
        all_first=num(f["gbif_all_plants_first_record_year"])
        anomalies=sum(y<all_first for y in years) if all_first is not None else 0
        if anomalies:status["host_year_before_all_plant_first_year"]+=anomalies
        cohort.append({
           "code":code,"level1":flora[code]["level1"],
           "ln_area":float(area[code]["ln_area"]),
           "ln_effort":math.log1p(total),
           "ln_pre1900_effort":math.log1p(pre),
           "ln_introduced_flora":math.log1p(flor),
           "ln_dated_species_count":math.log1p(len(years)),
           "n_dated_plants_in_original_fixed_pilot":len(years),
           "first_index_year_min":min(years),
           "median_earliest_year":statistics.median(years),
           "first_index_year_max":max(years),
           "all_plants_first_index_year":all_first if all_first is not None else "",
           "all_plants_indexed_records":int(total),
           "all_plants_pre1900_indexed_records":int(pre)
        })
    # An island with only one dated host has an unstable median year.
    # Keep an independent >=2-host result and geographic leave-one-region tests.
    at_least_two=[r for r in cohort if r["n_dated_plants_in_original_fixed_pilot"]>=2]
    two_hosts=fit(at_least_two,"total") if len(at_least_two)>=45 else {"status":"TOO_FEW_INDEPENDENT_REGIONS","n":len(at_least_two)}
    geographic_leverage=[]
    for lev in sorted({r["level1"] for r in cohort}):
        without=[r for r in cohort if r["level1"]!=lev]
        try:
            m=fit(without,"total")
            geographic_leverage.append({
                "excluded_level1":lev,"n":len(without),
                "beta":m["regression_year_change_per_ln1p_effort_unit"],
                "p_two":m["permutation"]["two_sided_centered_p"]
            })
        except RuntimeError as e:
            geographic_leverage.append({"excluded_level1":lev,"n":len(without),"status":"UNIDENTIFIABLE","reason":str(e)})
    output={
       "schema":"chocho_Q4_host_first_index_year_botanical_effort_geographic_audit_v0.2",
       "status":"OBSERVATION_INDEX_SAMPLING_BIAS_TEST_NOT_INTRODUCTION_DATING",
       "source":{
          "fixed_500_response_blind_host_island_pairs":len(dates),
          "dated_pairs":219,"islands_with_any_dated_host":len(grouped),
          "islands_with_complete_all_botanical_effort":len(cohort),
          "regions_with_at_least_two_dated_hosts":sum(r["n_dated_plants_in_original_fixed_pilot"]>=2 for r in cohort),
          "reconciliation_statuses":dict(status)
       },
       "primary":fit(cohort,"total"),
       "pre1900_effort_sensitivity":fit(cohort,"pre1900"),
       "sensitivity_2plus_dated_hosts_per_region":two_hosts,
       "leave_one_biogeographic_region_out":geographic_leverage,
       "limits":[
          "The first indexed GBIF host year is an upper bound on the true date of plant introduction, NOT plant arrival or naturalization.",
          "Higher observation effort can move recorded first years earlier even if plants arrived at identical dates.",
          "Only 219 of the response-blind 500 host-island pairs have GBIF dated records; all failed/missing queries are excluded, not imputed absent.",
          "All-botanical-effort count is a contemporaneous cumulative total that may include records after the first indexed host and can cause retrospective selection bias.",
          "Pre-1900 botanical total is an incomplete earlier-collection sensitivity, not an exact era-specific offset at host first record.",
          "Median among detected plant first indexed dates is a biased subsample of all introduced host taxa; the number of dated hosts is included as an analytic control.",
          "WGSRPD3 geographic region can contain multiple distinct physical islands; association is not a host-butterfly temporal risk model.",
          "Regional botanical effort is not independent of ports, settlement, plant introductions or true host distributions.",
          "This test is observational and descriptive and cannot estimate consumer lag or which of three establishment gates opened last."
       ],
    }
    (a.out_dir/"host_first_record_botanical_effort_audit.json").write_text(
         json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    with (a.out_dir/"island_first_year_effort_model_inputs.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(cohort[0]))
        writer.writeheader();writer.writerows(cohort)
    print(json.dumps({k:v for k,v in output.items() if k not in ("limits",)},ensure_ascii=False,indent=2),flush=True)
if __name__=="__main__":
    main()

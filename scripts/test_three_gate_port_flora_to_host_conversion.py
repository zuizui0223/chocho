#!/usr/bin/env python3
"""Source-backed island port-resource pathway screen, NOT causal mediation.

Test whether current WPI-indexed port regions contain more larval-host
taxa than expected from their *total introduced WCVP flora*. Then ask whether
port association with butterfly resource gain persists once current larval
host availability has been held constant.

All quantities are current, are at WGSRPD3 regional grain, and arise from
plant datasets and a frozen 239-butterfly HOSTS panel. Even with positive
results, this cannot establish a historical shipping route, butterfly
establishment, or any temporal order of these contemporaneous variables.
"""
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict
from pathlib import Path
import numpy as np

SEED=20261102
B=4999

def rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def num(s):
    try:
        v=float(s)
        if not math.isfinite(v):return None
        return v
    except (ValueError,TypeError):return None

def coeff(y,x,z):
    qr,_=np.linalg.qr(z,mode="reduced")
    yr=y-qr@(qr.T@y)
    xr=x-qr@(qr.T@x)
    den=float(xr@xr)
    return None if den<=1e-10 else float(xr@yr/den)

def regress(items, response, controls, run_perm=True):
    n=len(items)
    if n<60:raise RuntimeError(f"Small analysis cohort: {n}")
    levels=sorted({r["level1"] for r in items})
    if len(levels)<5:raise RuntimeError(f"Too few botanical biogeographic regions: {levels}")
    y=np.array([r[response] for r in items],dtype=float)
    x=np.array([r["log2_ports"] for r in items],dtype=float)
    all_controls=["const","ln_area",*controls]
    cols=[np.ones(n),np.array([r["ln_area"] for r in items])]
    cols.extend(np.array([r[k] for r in items]) for k in controls)
    for lev in levels[1:]:
        cols.append(np.array([1.0 if r["level1"]==lev else 0.0 for r in items]))
        all_controls.append("level1="+lev)
    z=np.column_stack(cols)
    rank=int(np.linalg.matrix_rank(z))
    if rank!=z.shape[1]:raise RuntimeError(f"Rank deficient controls {rank}/{z.shape[1]} {all_controls}")
    actual=coeff(y,x,z)
    if actual is None:raise RuntimeError("Port exposure not identifiable")
    output={"n":n,"level1_regions":len(levels),"response":response,
            "exposure":"log2(1+2023 WPI indexed ports)","controls":all_controls,
            "coefficient":actual,"per_port_doubling_multiplicative_change_in_1plus_outcome":math.expm1(actual)}
    # Pseudo-outcome checks rely on region-level independent geographic units,
    # not individual plant taxa or species×island cells.
    strata=defaultdict(list)
    bylev=defaultdict(list)
    for i,r in enumerate(items):bylev[r["level1"]].append(i)
    for lev,ii in bylev.items():
        median=float(np.median([items[i]["ln_area"] for i in ii]))
        for i in ii:
            area_label="small" if items[i]["ln_area"]<=median else "large"
            strata[lev+"|"+area_label].append(i)
    movable=[arr for arr in strata.values() if len(arr)>=2 and len({float(x[i]) for i in arr})>1]
    rng=random.Random(SEED)
    if run_perm:
        values=[]
        nullidentical=0
        # QR for controls fixed across all permutations.
        qr,_=np.linalg.qr(z,mode="reduced")
        yr=y-qr@(qr.T@y)
        for rep in range(B):
            xp=x.copy()
            for ids in movable:
                shuffled=list(ids)
                rng.shuffle(shuffled)
                xp[ids]=x[shuffled]
            if np.array_equal(xp,x):nullidentical+=1
            xr=xp-qr@(qr.T@xp)
            den=float(xr@xr)
            if den>1e-10:values.append(float(xr@yr/den))
        if len(values)<4900:raise RuntimeError("Too many invalid permuted port indices")
        perm=np.array(values)
        center=float(np.mean(perm))
        output["within_level1_and_area_half_label_permutation"]={
            "permutations":len(perm),"exchangeable_islands":len({i for x in movable for i in x}),
            "exchangeable_strata":len(movable),"same_exposure_draws":nullidentical,
            "null_mean":center,"null_q025":float(np.quantile(perm,0.025)),
            "null_q975":float(np.quantile(perm,0.975)),
            "observed_minus_null":actual-center,
            "p_two_centered":(1+int(np.sum(np.abs(perm-center)>=abs(actual-center)-1e-12)))/(len(perm)+1),
            "p_one_direction_positive":(1+int(np.sum(perm>=actual-1e-12)))/(len(perm)+1)}
    else:output["within_level1_and_area_half_label_permutation"]=None
    return output

def main():
    p=argparse.ArgumentParser()
    for name in ("all-flora","host-taxa","botanical-effort","out-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    flora=rows(a.all_flora)
    host=rows(a.host_taxa)
    effort=rows(a.botanical_effort)
    if len(flora)!=126 or len(host)!=123 or len(effort)!=98:
        raise RuntimeError(f"Frozen island cohorts changed: {len(flora)} {len(host)} {len(effort)}")
    botan={r["code"]:r for r in effort}
    flo={r["code"]:r for r in flora}
    if len(flo)!=126 or len(botan)!=98:raise RuntimeError("duplicate region keys")
    rset=[]
    for rr in host:
        code=rr["code"];f=flo[code]
        if f["level1"]!=rr["level1"] or int(f["wpi_port_count"])!=int(rr["ports"]):
            raise RuntimeError("Region/port mismatch "+code)
        alien=int(f["introduced"])
        hosttax=int(rr["introduced_larval_host_taxa"])
        added=int(rr["added"])
        native_flora=int(f["native"])
        native_resource=int(rr["native"])
        if alien<=0 or hosttax<0 or hosttax>alien or added<=0:
            raise RuntimeError(f"Host taxon not within WCVP introduced flora count: {code}")
        out={
          "code":code,"level1":rr["level1"],
          "ln_area":math.log(float(rr["area_km2"])),
          "log2_ports":math.log2(1+int(rr["ports"])),
          "WPI_ports":int(rr["ports"]),
          "introduced_flora":alien,"introduced_larval_host_taxa":hosttax,
          "added_butterfly_resource_cells":added,
          "log_alien":math.log1p(alien),
          "log_host":math.log1p(hosttax),
          "log_added":math.log1p(added),
          "log_native_flora":math.log1p(native_flora),
          "log_native_resource":math.log1p(native_resource),
          "host_taxa_per_alien_flora_species":hosttax/alien,
          "log_host_per_alien":math.log((hosttax+0.5)/(alien+0.5)),
          "has_valid_botanical_effort":int(code in botan)
        }
        if code in botan:
            x=botan[code]
            if int(x["n_introduced_flora"])!=alien or int(x["host_taxa"])!=hosttax or int(x["WPI_ports"])!=int(rr["ports"]):
                raise RuntimeError("GBIF effort+flora full join mismatch "+code)
            out["ln_botanical_effort"]=math.log1p(float(x["gbif_all_plant_indexed_records"]))
            out["ln_pre1900_botanical_effort"]=math.log1p(float(x["gbif_pre1900_plant_records"]))
        rset.append(out)
    if len(rset)!=123 or sum(x["has_valid_botanical_effort"] for x in rset)!=98:
        raise RuntimeError("Effective cohorts drift")
    baseline=[
      ("M0_port_and_entire_alien_flora",rset,"log_alien",["log_native_flora"]),
      ("M1_port_host_taxa_conditional_on_all_alien_flora",rset,"log_host",["log_alien","log_native_flora"]),
      ("M1b_port_relative_host_taxa_per_alien_flora",rset,"log_host_per_alien",["log_native_flora"]),
      ("M2_port_added_butterfly_resources_conditional_on_larval_hosts",rset,"log_added",["log_host","log_native_resource"]),
    ]
    filtered=[r for r in rset if r["has_valid_botanical_effort"]]
    sensitive=[
      ("M0eff_port_entire_alien_flora",filtered,"log_alien",["log_native_flora","ln_botanical_effort"]),
      ("M1eff_port_host_taxa_given_alien_flora",filtered,"log_host",["log_alien","log_native_flora","ln_botanical_effort"]),
      ("M1pre1900eff_port_host_taxa_given_alien_flora",filtered,"log_host",["log_alien","log_native_flora","ln_pre1900_botanical_effort"]),
      ("M2eff_port_added_butterfly_resources_given_hosts",filtered,"log_added",["log_host","log_native_resource","ln_botanical_effort"]),
    ]
    reports={}
    for name,sub,y,z in baseline+sensitive:
        reports[name]=regress(sub,y,z)
        print(name,json.dumps({"n":reports[name]["n"],"beta":reports[name]["coefficient"],
             "two_sided":reports[name]["within_level1_and_area_half_label_permutation"]["p_two_centered"]}),flush=True)
    primary=reports["M1_port_host_taxa_conditional_on_all_alien_flora"]
    effort_primary=reports["M1eff_port_host_taxa_given_alien_flora"]
    # Island leverage check on primary with *original* controls fitted again.
    leave=[]
    for skip in sorted({x["code"] for x in rset}):
        sub=[x for x in rset if x["code"]!=skip]
        try:
            fit=regress(sub,"log_host",["log_alien","log_native_flora"],run_perm=False)
            leave.append({"excluded":skip,"beta":fit["coefficient"]})
        except RuntimeError:
            leave.append({"excluded":skip,"beta":None})
    values=[r["beta"] for r in leave if r["beta"] is not None]
    summary={
      "schema":"chocho_island_ports_alien_flora_larval_host_conversion_v0.1",
      "status":"CROSS_SECTIONAL_OBSERVATIONAL_RESOURCE_CHAIN_NOT_TRANSPORT_OR_COLONIZATION",
      "source":{"WCVP_all_flora_accepted_species_regions":126,
                "WCVP_HOSTS_selected_resource_active_islands":123,
                "successful_GBIF_total_plant_effort_islands":98,
                "butterflies_in_HOSTS_fixed_panel":239,
                "modern_2023_WPI_ports_from_published_3630_port_inventory":True},
      "hypothesis":"Does indexed port availability predict the fraction/number of larval-host plants beyond the entire introduced flora, and added butterfly resource cells beyond the larval host flora?",
      "naive_means":{
        "indexed_port_islands":sum(x["WPI_ports"]>0 for x in rset),
        "no_indexed_port_islands":sum(x["WPI_ports"]==0 for x in rset),
        "mean_larval_host_taxa_with_port":float(np.mean([x["introduced_larval_host_taxa"] for x in rset if x["WPI_ports"]>0])),
        "mean_larval_host_taxa_no_port":float(np.mean([x["introduced_larval_host_taxa"] for x in rset if x["WPI_ports"]==0])),
        "mean_host_taxa_per_all_introduced_flora_with_port":float(np.mean([x["host_taxa_per_alien_flora_species"] for x in rset if x["WPI_ports"]>0])),
        "mean_host_taxa_per_all_introduced_flora_no_port":float(np.mean([x["host_taxa_per_alien_flora_species"] for x in rset if x["WPI_ports"]==0]))
      },
      "primary_model_predeclared":"M1_port_host_taxa_conditional_on_all_alien_flora",
      "primary_effort_model_predeclared":"M1eff_port_host_taxa_given_alien_flora",
      "model_results":reports,
      "sensitivity_leave_one_island_out_primary":{
          "n":len(values),"positive":sum(x>0 for x in values),
          "minimum":min(values),"maximum":max(values),
          "note":"Range of fitted coefficients, NOT a confidence interval",
          "extreme_cases":sorted([x for x in leave if x["beta"] is not None],
                            key=lambda x:x["beta"])[:3]+
                          sorted([x for x in leave if x["beta"] is not None],
                            key=lambda x:x["beta"])[-3:]
      },
      "sample_bias":{
          "effort_complete_islands":len(filtered),
          "effort_missing_islands":len(rset)-len(filtered),
          "cannot_assume_missing_at_random":True
      },
      "limits":[
         "Current 2023 indexed WPI ports do not measure historical ship connectivity or actual larval/egg transport.",
         "WCVP native/introduced plant status is current and not the age of introductions.",
         "The number of HOSTS-matched larval host taxa is derived from selected butterfly species; HOSTS records and underlying accepted plant taxa differ in taxonomic resolution from all WCVP accepted flora species. 'Host taxa per all flora species' is a descriptive index, not an exact plant-species proportion.",
         "Resource gain in butterfly x island cells is algebraically related to the number of introduced hosts, so M2 is a structural comparison not independent confirmation.",
         "Area, Level1 and botanical GBIF plant record count do not control economy, human population, tourism, shipping volume, land use, network connectivity, or spatial autocorrelation.",
         "GBIF botanical reporting itself rises with ports; only 98 island polygons have a successful count and this subset differs geographically/area-wise.",
         "Within region/area label permutations are descriptive and cannot establish causal ship-mediated plant introduction.",
         "These analyses do not measure butterfly island establishment, lag, flight capacity, climate suitability, or the last gate to open.",
         "Do not interpret conditioning on contemporaneous alien flora as a causal mediation model."
      ]
    }
    outfile=a.out_dir/"island_port_flora_host_conversion_audit.json"
    outfile.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("SOURCE_AUDIT_OK "+str(outfile))
if __name__=="__main__":
    main()

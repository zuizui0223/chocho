#!/usr/bin/env python3
"""Source-locked island-level moderation of WPI port–introduced plant association.

Scientific unit is the WGSRPD3 botanical ISLAND REGION, NOT a vessel, plant
immigration event, butterfly establishment, or species×island repeat.
Two pre-specified moderator axes: ln area and ln native WCVP flora.
Two outcomes: ln1p(all introduced WCVP flora), and ln1p(known introduced
larval host taxa) CONDITIONED on ln1p(all introduced flora).

This analysis is observational and exploratory. Holm adjustment and
within-Level1-area-half Freedman–Lane residual permutations are descriptive,
not claims of causal randomization or phylogenetic independence.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

SOURCE_SHA={
    "flora":"53b618fa51a1a418ccc6e4363f4c866c389c90cc",
    "hosts":"cfe82bd8223ea453c015c47cd9659a935d5ba3ab",
    "effort":"23f5e5c7f7c875af41ccb5c6094cb4f48f1071fe",
}
SEED=20261117
N_PERM=1499
N_BOOT=999

def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as fh:
        return list(csv.DictReader(fh))

def key(rows,column):
    d={}
    for r in rows:
        if r[column] in d:raise RuntimeError(f"Duplicate source key {r[column]}")
        d[r[column]]=r
    return d

def vec_data(flora,host,eff):
    if (len(flora),len(host),len(eff))!=(126,123,98):
        raise RuntimeError(f"Original source island cohorts drifted: {len(flora)}, {len(host)}, {len(eff)}")
    hosts=key(host,"code"); effs=key(eff,"code")
    codes=set()
    values=[]
    for f in flora:
        code=f["code"]
        if code in codes:raise RuntimeError("Duplicate WCVP island code")
        codes.add(code)
        total=int(f["total"]);native=int(f["native"]);alien=int(f["introduced"])
        if native+alien!=total or total<0:raise RuntimeError("Invalid native+introduced whole flora counts")
        area=float(f["area_km2"]);ports=int(f["wpi_port_count"])
        if area<=0 or ports<0 or abs(math.log(area)-float(f["log_area"]))>1e-8:
            raise RuntimeError("Area/WPI source changed")
        d=dict(code=code,level1=str(f["level1"]),area=area,lnarea=math.log(area),
               native_flora=native,alien_flora=alien,lnnative=math.log1p(native),
               lnalien=math.log1p(alien),ports=ports,ln2ports=math.log2(ports+1),
               has_host=code in hosts,has_effort=code in effs)
        if code in hosts:
            h=hosts[code]
            if (str(h["level1"])!=d["level1"] or int(h["ports"])!=ports or
                abs(float(h["area_km2"])-area)>1e-6*max(1,area)):
                raise RuntimeError(f"Original HOSTS resource island doesn't match WCVP flora and port: {code}")
            hosttax=int(h["introduced_larval_host_taxa"])
            if alien<=0 or hosttax<0 or hosttax>alien:
                raise RuntimeError(f"Impossible accepted introduced host taxa {code}")
            d.update(hosttax=hosttax,lnhost=math.log1p(hosttax),
                     log_native_resource=math.log1p(int(h["native"])),
                     added_resources=int(h["added"]))
        if code in effs:
            e=effs[code]
            if (str(e["level1"])!=d["level1"] or int(e["n_native_flora"])!=native or
                int(e["n_introduced_flora"])!=alien or int(e["WPI_ports"])!=ports):
                raise RuntimeError(f"GBIF botanical effort source not aligned to WCVP plant/port island {code}")
            d["lneffort"]=math.log1p(float(e["gbif_all_plant_indexed_records"]))
            d["lnpre1900effort"]=math.log1p(float(e["gbif_pre1900_plant_records"]))
        values.append(d)
    if len(codes)!=126 or set(hosts)-codes or set(effs)-codes:
        raise RuntimeError("Global island/host/effort support not nested")
    if sum(r["has_host"] for r in values)!=123 or sum(r["has_host"] and r["has_effort"] for r in values)!=98:
        raise RuntimeError("Original 123/98 resource-botanical-effort nesting changed")
    return values

def build_design(items,axis,*,specific=False,effort=False):
    n=len(items)
    levs=sorted({r["level1"] for r in items},key=int)
    if n<70 or len(levs)<7:raise RuntimeError("Insufficient units or regional contrast for moderator")
    Y=np.array([r["lnhost" if specific else "lnalien"] for r in items],dtype=float)
    port=np.array([r["ln2ports"] for r in items])
    area=np.array([r["lnarea"] for r in items]);native=np.array([r["lnnative"] for r in items])
    const=np.ones(n)
    marea=(area-REFERENCE["area_median"])/REFERENCE["area_sd"]
    mnative=(native-REFERENCE["native_median"])/REFERENCE["native_sd"]
    # Main port exposure is included in both the null and alternative.
    # This tests MODERATION, not whether ports associate with flora overall.
    X=[const,area,native,port]
    names=["intercept","ln_area_km2","ln1p_native_flora","log2_1plus_WPI_ports"]
    if specific:X.append(np.array([r["lnalien"] for r in items]));names.append("ln1p_all_introduced_flora")
    if effort:X.append(np.array([r["lneffort"] for r in items]));names.append("ln1p_GBIF_all_plant_records")
    for lev in levs[1:]:
        X.append(np.array([float(r["level1"]==lev) for r in items]))
        names.append("TDWG_Level1_"+lev)
    reduced=np.column_stack(X)
    if np.linalg.matrix_rank(reduced)!=reduced.shape[1]:raise RuntimeError("Source main-effect model rank deficient")
    m=marea if axis=="island_area" else mnative
    full=np.column_stack([reduced,port*m])
    if np.linalg.matrix_rank(full)!=full.shape[1]:raise RuntimeError("Interaction not identifiable")
    return Y,reduced,full,names,m

def fit(y,x):
    b=np.linalg.lstsq(x,y,rcond=None)[0]
    r=y-x@b
    return b,r

def quantile(a,q):
    return float(np.quantile(a,q))

def blocks(items,within_area_half=False):
    by=defaultdict(list)
    if within_area_half:
        group=defaultdict(list)
        for i,r in enumerate(items):group[r["level1"]].append(i)
        med={lv:float(np.median([items[i]["lnarea"] for i in indices]))
             for lv,indices in group.items()}
        for i,r in enumerate(items):
            k=r["level1"]+"|"+str(int(r["lnarea"]>med[r["level1"]]))
            by[k].append(i)
    else:
        for i,r in enumerate(items):by[r["level1"]].append(i)
    return [np.array(v,dtype=int) for v in by.values() if len(v)>=2]

def describe(items,axis,*,specific=False,effort=False,permutations=True,bootstrap=True):
    y,reduced,full,names,m=build_design(items,axis,specific=specific,effort=effort)
    b,r=fit(y,full)
    beta=float(b[-1]);main=float(b[names.index("log2_1plus_WPI_ports")])
    # Descriptive partial R2 from the additional interaction.
    b0,r0=fit(y,reduced)
    partial_r2=1-float(r@r)/float(r0@r0) if r0@r0>1e-12 else None
    den=np.sum(r*r)
    rng=np.random.default_rng(SEED+len(items)+(1 if specific else 0)+(15 if effort else 0)+(42 if axis=="native_flora" else 0))
    boot=[]
    if bootstrap:
        # Preserve number of botanical biogeographic Level1 groups and their
        # original sizes; resample ISLANDS, never species records.
        groups=blocks(items)
        for _ in range(N_BOOT):
            idx=np.concatenate([rng.choice(g,len(g),replace=True) for g in groups])
            yb,_,xb,_,_=build_design([items[int(i)] for i in idx],axis,specific=specific,effort=effort)
            if np.linalg.matrix_rank(xb)!=xb.shape[1]:continue
            bb,_=fit(yb,xb)
            boot.append(float(bb[-1]))
    ci=[quantile(boot,.025),quantile(boot,.975)] if len(boot)>=400 else None

    perm=[]
    if permutations:
        # Freedman-Lane conditional on the existing port MAIN effect:
        # fitted null+permuted residuals, exchangeable only within
        # independent, source-defined TDWG Level1 × area half strata.
        fit_null,reduced_resid=fit(y,reduced)
        mu=reduced@fit_null
        groups=blocks(items,within_area_half=True)
        for _ in range(N_PERM):
            yp=mu+reduced_resid
            yp=yp.copy()
            for g in groups:
                yp[g]=mu[g]+reduced_resid[rng.permutation(g)]
            bp,_=fit(yp,full)
            perm.append(float(bp[-1]))
    pval=(1+int(np.sum(np.abs(np.asarray(perm))>=abs(beta)-1e-12)))/(len(perm)+1) if perm else None
    # Strictly do not make causal claim from permuted observational records.
    leave=[]
    for level in sorted({r["level1"] for r in items},key=int):
        sub=[r for r in items if r["level1"]!=level]
        try:
            yy,_,xx,_,_=build_design(sub,axis,specific=specific,effort=effort)
            coeff=float(fit(yy,xx)[0][-1])
            leave.append({"excluded_TDWG_Level1":level,"n":len(sub),"beta":coeff})
        except RuntimeError as e:
            leave.append({"excluded_TDWG_Level1":level,"status":"NOT_IDENTIFIABLE","reason":str(e)})
    eligible=[x["beta"] for x in leave if "beta" in x]
    out={
      "n_islands":len(items),"level1_groups":len({r["level1"] for r in items}),
      "model":("HOSTS_known_larval_host_taxa_given_all_alien_flora" if specific else "ALL_WCVP_introduced_flora"),
      "moderator":axis,"includes_GBIF_botanical_effort":effort,
      "main_port_beta_at_moderator_reference":main,
      "port_times_moderator_standardized_beta":beta,
      "percent_change_per_port_log2_unit_at_moderator_1SD_above_vs_reference":
         float(math.expm1(beta)*100),
      "interaction_partial_R2_relative_to_same_main_effect_model":partial_r2,
      "island_within_Level1_bootstrap_95pct_descriptive":ci,
      "valid_bootstrap_draws":len(boot),
      "Freedman_Lane_permute_residual_within_Level1_and_area_half_null_95pct":
         [quantile(perm,.025),quantile(perm,.975)] if perm else None,
      "Freedman_Lane_two_sided_exploratory_p":pval,
      "leave_one_TDWG_Level1_out":leave,
      "leave_one_Level1_beta_range":[min(eligible),max(eligible)] if eligible else None,
      "limitations":[
       "This is an interaction in contemporary botanical distribution data, not real cargo movement or an observed introduction.",
       "The moderator was selected for scientific plausibility before the model fit; tests are exploratory, not registered independent confirmations.",
       "The Level1+area-half residual permutation is diagnostic and requires exchangeability not guaranteed by spatial, temporal or human-census selection.",
       "Bootstrap draws resample botanical regions conditional on their TDWG Level1, not independent plant species or vessels.",
       "WCVP accepted native/introduced status and WPI 2023 indexed port counts have distinct, unknown historical timing.",
       "Within-region area/natural richness and GBIF plant indexing effort do not exhaust geographic or observational confounding.",
      ]
    }
    return out

def holm(values):
    ps=sorted((v,i) for i,v in enumerate(values))
    out=[None]*len(ps);last=0.
    for j,(p,i) in enumerate(ps):
        last=max(last,min(1.,(len(ps)-j)*p))
        out[i]=last
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flora",required=True,type=Path)
    ap.add_argument("--hosts",required=True,type=Path)
    ap.add_argument("--effort",required=True,type=Path)
    ap.add_argument("--out-dir",required=True,type=Path)
    a=ap.parse_args()
    for name,path in (("flora",a.flora),("hosts",a.hosts),("effort",a.effort)):
        actual=gitblob(path)
        if actual!=SOURCE_SHA[name]:
            raise RuntimeError(f"Fail-closed original {name} SHA mismatch {actual}")
    z=vec_data(rows(a.flora),rows(a.hosts),rows(a.effort))
    global REFERENCE
    REFERENCE={
      "area_median":float(np.median([x["lnarea"] for x in z])),
      "area_sd":float(np.std([x["lnarea"] for x in z],ddof=1)),
      "native_median":float(np.median([x["lnnative"] for x in z])),
      "native_sd":float(np.std([x["lnnative"] for x in z],ddof=1)),
    }
    h=[r for r in z if r["has_host"]]
    e=[r for r in h if r["has_effort"]]
    # Original model M0 on same 123 islands is a mandatory validity gate.
    yy,red,full,_,_=build_design(h,"island_area",specific=False,effort=False)
    coef,_=fit(yy,red)
    orig=float(coef[3])
    if abs(orig-0.23347121981392502)>1e-8:
        raise RuntimeError(f"Original 123-island port-native-allflora source estimate not reproduced {orig}")
    tests=[]
    for scope,sub,eff in (("whole_botanical_panel",z,False),("matched_host_source_panel",h,False),
                          ("GBIF_effort_matched_host_source_panel",e,True)):
        for response in (False,True):
            if response and scope=="whole_botanical_panel":continue
            for axis in ("island_area","native_flora"):
                d=describe(sub,axis,specific=response,effort=eff)
                d["source_cohort"]=scope
                tests.append(d)
                print("PORT_FLORA_MODERATOR",json.dumps({
                  "scope":scope,"outcome":d["model"],"axis":axis,"n":d["n_islands"],
                  "coef":d["port_times_moderator_standardized_beta"],
                  "region_perm_p":d["Freedman_Lane_two_sided_exploratory_p"],
                  "bootstrap95":d["island_within_Level1_bootstrap_95pct_descriptive"],
                  "loo_range":d["leave_one_Level1_beta_range"],
                }),flush=True)
    # Holm exploratory adjustment only for the four primary size/native and
    # two botanical outcomes in the exact SAME 123 host-active islands.
    primary=[t for t in tests if t["source_cohort"]=="matched_host_source_panel"]
    if len(primary)!=4:raise RuntimeError("Original four primary flora-vs-host moderator contrasts changed")
    ps=holm([x["Freedman_Lane_two_sided_exploratory_p"] for x in primary])
    for i,t in enumerate(primary):
        t["Holm_adjusted_exploratory_p_across_four_matched_cohort_tests"]=ps[i]
    a.out_dir.mkdir(parents=True,exist_ok=True)
    report={
      "schema":"chocho_original_WCVP_all_flora_host_specific_port_island_area_native_moderation_v1",
      "status":"SOURCE_VERIFIED_OBSERVATIONAL_MODERATOR_SCREEN_NOT_SHIPPING_OR_BUTTERFLY_COLONISATION",
      "original_Git_blob_SHA":SOURCE_SHA,
      "n_original_all_WCVP_islands":126,
      "n_original_HOSTS_active_islands":123,
      "n_original_GBIF_botanical_effort_matched_islands":98,
      "n_original_TDWG_Level1_regions":len(set(r["level1"] for r in z)),
      "original_source_123_island_M0_port_alien_flora_beta_reproduced":orig,
      "moderator_reference":REFERENCE,
      "permutations_per_test":N_PERM,"island_bootstrap_attempts_per_test":N_BOOT,
      "one_paper_note":"Existing globally positive +54.9% resource; islands +118.9%; known-butterfly GBIF 4.12% forecast remain unchanged and are different estimands.",
      "tests":tests,
      "causal_warning":"Neither the regional port coefficient nor its interactions reconstruct historical vessels, individual seed movements, new plants, host suitability, butterfly feeding, populations or establishment.",
    }
    out=a.out_dir/"source_verified_port_alien_flora_island_moderators.json"
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("PORT_ALIEN_FLORA_MODERATION_SOURCE_GATE",json.dumps({
      "status":report["status"],"source_123_model_beta":orig,
      "matched_cohort": [{
        "outcome":x["model"],"axis":x["moderator"],
        "beta":x["port_times_moderator_standardized_beta"],
        "Holm_p":x["Holm_adjusted_exploratory_p_across_four_matched_cohort_tests"]}
        for x in primary],
      "source_test_count":len(tests)
    }),flush=True)

if __name__=="__main__":main()

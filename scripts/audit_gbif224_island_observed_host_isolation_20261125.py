#!/usr/bin/env python3
"""Independent recorded-butterfly bridge to the frozen 239 x 126 host map.
Nonrecords are NOT absences; current host maps NOT dated plant arrivals."""
import argparse, csv, gzip, hashlib, json, sys
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t
from audit_island_port_isolation_239_20261120 import make_frame, design

FOLDER=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
SHA={"cells":"f1bd777e2faffdeb0d689b93dca193e8612385fa264dd3bbdb6e104c075f6525",
     "effort":"310f8f58ff1b4bedc9bf3f9b39b85c4fb11248916f2b13793795af7ec8925e82",
     "matrix":"ec6a5530d65ac0b04f230e0c51fd68bf543555641a17b94795996bf11ce257d0",
     "receipt":"ccda0bc669fb5ffe7326175fcc0501dbb1447790925c06add92501d1d30f255b"}
PERIODS={"all":(2001,2024),"prior":(2001,2010),"late":(2018,2024)}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def median(x):return float(x.median()) if len(x) else None
def fraction(a,b):return float(a/b) if b else None
def z(s):return (s-s.mean())/s.std(ddof=1)
def describe(f,tag):
    n,i,zz=[int(f[tag+"_"+k].sum()) for k in ("N","I","Z")]
    pn,pi=[int(f["potential_"+k].sum()) for k in ("N","I")]
    active=f[(f[tag+"_N"]+f[tag+"_I"])>0]
    return {"island_regions":len(f),"resource_record_active_regions":len(active),
            "introduced_only_recorded_regions":int((f[tag+"_I"]>0).sum()),
            "recorded_N":n,"recorded_I":i,"recorded_Z":zz,
            "recorded_I_share_known":fraction(i,n+i),"recorded_I_share_all":fraction(i,n+i+zz),
            "median_equal_island_I_share_known":median(active[tag+"_I"]/(active[tag+"_I"]+active[tag+"_N"])),
            "potential_N":pn,"potential_I":pi,"potential_I_share":fraction(pi,pn+pi),
            "indexed_fraction_I":fraction(i,pi),"indexed_fraction_N":fraction(n,pn)}
def models(f,tag,period):
    f=f[(f[tag+"_N"]+f[tag+"_I"]>0)&(f[period+"_days"]>0)].copy().reset_index(drop=True)
    if len(f)<40 or f.level1.nunique()<5:return {"status":"INSUFFICIENT_INDEPENDENT_REGION_SUPPORT","n":len(f),"level1":int(f.level1.nunique())}
    X=design(f); X["effort"]=z(np.log1p(f[period+"_days"]))
    obs=np.log((f[tag+"_I"]+.5)/(f[tag+"_N"]+.5))
    potential=np.log((f.potential_I+.5)/(f.potential_N+.5))
    def fit(y,x):
        if np.linalg.matrix_rank(x)!=x.shape[1] or len(y)-x.shape[1]<10:return {"status":"INSUFFICIENT_FULL_RANK"}
        m=sm.OLS(y,x).fit(cov_type="HC3")
        cl=sm.OLS(y,x).fit(cov_type="cluster",cov_kwds={"groups":f.level1})
        b=float(m.params["distance"]);se=float(cl.bse["distance"]);df=int(f.level1.nunique())-1
        lo,hi=m.conf_int().loc["distance"]
        return {"n":len(y),"level1_clusters":df+1,"beta":b,"HC3_p":float(m.pvalues["distance"]),
                "HC3_95":[float(lo),float(hi)],"cluster_t_p":float(2*t.sf(abs(b/se),df)) if se>0 else None,
                "cluster_t_95":[b-t.ppf(.975,df)*se,b+t.ppf(.975,df)*se]}
    return {"observed_share_logodds":fit(obs,X),"potential_share_logodds":fit(potential,X),
            "observed_minus_potential_logodds":fit(obs-potential,X),
            "observed_share_conditioning_on_potential":fit(obs,X.assign(potential=z(potential)))}
def main():
    p=argparse.ArgumentParser()
    for k in ("annual-root","resource-root","current-root","tdwg-geojson","out-dir"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    pp={ "cells":a.annual_root/FOLDER/"gbif_global224_annual_v1.species_region_year.csv.gz",
         "effort":a.annual_root/FOLDER/"gbif_global224_annual_v1.region_year_effort.csv",
         "matrix":a.resource_root/FOLDER/"beta_baseline_v3.matrix.csv.gz",
         "receipt":a.resource_root/FOLDER/"beta_baseline_v3.json"}
    for k,v in pp.items():need(digest(v)==SHA[k],"Source SHA drift: "+k)
    annual=pd.read_csv(pp["cells"],dtype={"wgsrpd3_code":str})
    eff=pd.read_csv(pp["effort"],dtype={"wgsrpd3_code":str})
    need(len(annual)==91882 and len(eff)==5390 and annual.indexed_gbif_records.sum()==13587452,"Source annual totals")
    need(not annual.duplicated(["butterfly_species","wgsrpd3_code","year"]).any() and not eff.duplicated(["wgsrpd3_code","year"]).any(),"Duplicate source cells")
    by=annual.groupby(["wgsrpd3_code","year"]).indexed_gbif_records.agg(["sum","size"]).reset_index()
    chk=by.merge(eff,on=["wgsrpd3_code","year"],validate="one_to_one")
    need(len(chk)==5390 and (chk["sum"]==chk.indexed_candidate_butterfly_records).all() and (chk["size"]==chk.indexed_candidate_butterfly_species).all(),"Annual effort reconstruction")
    taxa=json.loads(pp["receipt"].read_text())["resource_matrix"]["species_order"]
    mat=pd.read_csv(pp["matrix"],dtype=str)
    need(len(taxa)==len(set(taxa))==239 and len(mat)==369 and mat.wgsrpd3_code.is_unique,"Matrix source")
    masks={r.wgsrpd3_code:(int(r.native_bits_hex,16),int(r.added_bits_hex,16),int(r.current_bits_hex,16)) for r in mat.itertuples()}
    need(all(n&i==0 and n|i==c for n,i,c in masks.values()),"Nested resources")
    need([sum(v[j].bit_count() for v in masks.values()) for j in range(3)]==[26530,14553,41083],"Source matrix totals")
    species=sorted(annual.butterfly_species.unique()); need(len(species)==224 and set(species)<=set(taxa),"Common cohort not exact 224")
    ix={s:1<<taxa.index(s) for s in species}
    original_regions=pd.read_csv(a.current_root/FOLDER/"gbif_global224_observed_beta_region_holdout_v1.csv")
    prim=set(original_regions.loc[original_regions.model=="M0_all_flora_effort_prior","wgsrpd3_code"])
    need(len(prim)==169,"Original 169-region roster")
    f=make_frame(a.current_root/"results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv",
                 a.current_root/"results/island-mainland-resource-q1/region_metrics.csv",a.tdwg_geojson)
    need(len(f)==126 and set(f.code)<=set(masks),"Island source roster")
    for row in f.itertuples():
        m=masks[row.code]
        need(m[0].bit_count()==row.native_host_239 and m[2].bit_count()==row.current_host_239,"Original island count mismatch")
    q=f.continent_distance_km.quantile([.25,.75])
    f["distance_quartile"]=np.where(f.continent_distance_km<=q.iloc[0],"near",np.where(f.continent_distance_km>=q.iloc[1],"far","middle"))
    need(f.distance_quartile.value_counts().to_dict()=={"middle":62,"near":32,"far":32},"Quartiles changed")
    # Dates measure record recurrence, not breeding. Keep all positive years,
    # and never use unpublished individual occurrence identities.
    rec=annual.groupby(["wgsrpd3_code","butterfly_species"]).year.agg(set).to_dict()
    ef=eff.groupby("wgsrpd3_code")
    results={}
    outrows={}
    monarch={}
    for label,panel in (("all224",species),("without_monarch",[s for s in species if s!="Danaus plexippus"])):
        rows=[];mask=sum(ix[s] for s in panel)
        for r in f.itertuples():
            n,i,c=masks[r.code];e=ef.get_group(r.code) if r.code in ef.groups else eff.iloc[:0]
            d={"code":r.code,"primary169":r.code in prim,"potential_N":(n&mask).bit_count(),"potential_I":(i&mask).bit_count(),"potential_Z":len(panel)-(c&mask).bit_count()}
            for period,(lo,hi) in PERIODS.items():
                d[period+"_days"]=int(e.loc[e.year.between(lo,hi),"distinct_gbif_datasetkey_calendar_days"].sum())
                for k in (1,2,3):
                    counts=Counter()
                    for s in panel:
                        years=rec.get((r.code,s),set()); yrs=sum(lo<=y<=hi for y in years)
                        if yrs>=k:counts["N" if n&ix[s] else "I" if i&ix[s] else "Z"]+=1
                    for st in ("N","I","Z"):d[f"{period}{k}_{st}"]=counts[st]
            rows.append(d)
        g=f.merge(pd.DataFrame(rows),on="code",validate="one_to_one")
        cohort={"all126":g,"original169_islands":g[g.primary169],
                "all_period_min20_days":g[g.all_days>=20],"late_period_min20_days":g[g.late_days>=20]}
        sums={name:{tag:describe(d,tag) for tag in ("all1","all2","all3","late1","late2","late3","prior1")} for name,d in cohort.items()}
        groups={name:{q:{tag:describe(d[d.distance_quartile==q],tag) for tag in ("all1","all2","late1","late2")} for q in ("near","far")} for name,d in cohort.items()}
        fits={}
        for name,d in cohort.items():
            fits[name]={tag:models(d,tag,period) for tag,period in (("all1","all"),("all2","all"),("late1","late"),("late2","late"))}
        results[label]={"summary":sums,"quartiles":groups,"models":fits}
        outrows[label]=g
    n_codes=[k for k in f.code if masks[k][0]&ix["Danaus plexippus"]]
    i_codes=[k for k in f.code if masks[k][1]&ix["Danaus plexippus"]]
    monarch={"native_regions":len(n_codes),"introduced_only_regions":len(i_codes),
             "recorded_introduced_only_regions":[k for k in i_codes if len(rec.get((k,"Danaus plexippus"),set()))>=1],
             "recurrently_recorded_introduced_only_regions":[k for k in i_codes if len(rec.get((k,"Danaus plexippus"),set()))>=2]}
    # Independent source-derived calculation in the conversation: lock the count-level replay.
    base=results["all224"]["summary"]["all126"]["all1"]
    need((base["recorded_N"],base["recorded_I"],base["recorded_Z"],base["introduced_only_recorded_regions"])==(916,193,62,69),"Independent direct-count replay mismatch")
    need(monarch["native_regions"]==22 and monarch["introduced_only_regions"]==42 and len(monarch["recorded_introduced_only_regions"])==20 and len(monarch["recurrently_recorded_introduced_only_regions"])==18,"Monarch direct-count replay mismatch")
    report={"status":"SOURCE_GATED_RECORDED_OPPORTUNITY_ASSOCIATIONS_NOT_OCCUPANCY",
        "source_sha256":SHA,"GBIF_DOI":"10.15468/dl.umqmsf","original_species":224,"frozen_resource_species":239,"island_regions":126,
        "periods":PERIODS,"monarch":monarch,"panels":results,
        "limits":["Current plant maps are not historical arrival maps.","Recorded species-region combinations are not larval-use, dependence, or establishment evidence.",
                  "Z means no documented mapped host in the original source, not biological absence of suitable hosts.",
                  "The original 169-region eligibility cohort contains only a small set of islands; far quartile may be poorly supported.",
                  "Equal-island log-odds models use positive recorded known-host units; effort proxies are observational, not detection probabilities.",
                  "Exploratory thresholds and multiple endpoints; no significance-based cohort selection."]}
    (a.out_dir/"recorded_island_host_isolation.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    outrows["all224"].to_csv(a.out_dir/"island_aggregates.csv",index=False)
    print(json.dumps(report,allow_nan=False))
if __name__=="__main__":main()

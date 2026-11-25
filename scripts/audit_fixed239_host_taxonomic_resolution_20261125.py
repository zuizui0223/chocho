#!/usr/bin/env python3
"""Fixed-239, fixed-126 taxonomic-resolution audit; no inferred novel hosts."""
import argparse, hashlib, json
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t
from audit_island_port_isolation_239_20261120 import design, make_frame, y_standardized
from audit_island_port_spatial_239_20261120 import make_weights, preselected_mem

EXPECTED = {
 "edges": "0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9",
 "native": "b66350a5d965e707cbfc5cb1de1d6fb72a62b6a01e5632eed2267b80cb02d091",
 "current": "bdf3ad670c435f570a715b23cb31ca7ab6add6e3241bee2105dc101fc675ebb7",
}
def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    return pd.read_csv(p, dtype=str, keep_default_na=False)
def require(ok, why):
    if not ok:
        raise RuntimeError(why)
def ranks(ids, tax):
    return {k:{tax[p][k] for p in ids} for k in ("genus","family")}
def aggregate(hosts, ids, tax):
    local = ranks(ids,tax)
    result = {}
    for rank in ("species","genus","family"):
        target = ids if rank=="species" else local[rank]
        known = [h[rank] for h in hosts.values()]
        result[rank] = sum(bool(target & h) for h in known)
        result[rank+"_fraction"] = float(np.mean([len(target & h)/len(h) for h in known]))
    require(result["species"]<=result["genus"]<=result["family"]<=len(hosts),
            "Taxonomic nesting failed")
    return result
def term(y, X, groups):
    require(np.isfinite(y).all(), "Nonfinite response")
    require(np.linalg.matrix_rank(X)==X.shape[1], "Rank deficient design")
    m=sm.OLS(y,X).fit(cov_type="HC3")
    cl=sm.OLS(y,X).fit(cov_type="cluster",cov_kwds={"groups":groups})
    beta=float(m.params["distance"])
    se=float(cl.bse["distance"])
    df=len(set(groups))-1
    lo,hi=m.conf_int().loc["distance"]
    return {"beta":beta,"HC3_p":float(m.pvalues["distance"]),
            "HC3_95":[float(lo),float(hi)],
            "cluster_t_df":df,"cluster_t_p":float(2*t.sf(abs(beta/se),df)) if se>0 else None,
            "cluster_t_95":[beta-t.ppf(.975,df)*se,beta+t.ppf(.975,df)*se]}
def suite(f, prefix, X, Xmem):
    ys={}
    for origin in ("native","current"):
        for rank in ("species","genus","family"):
            col=prefix+origin+"_"+rank
            ys[origin+"_"+rank]=y_standardized(f,col)
            # Equal-weight mean fraction of the documented taxon list present per butterfly.
            ys[origin+"_"+rank+"_fraction"]=f[col+"_fraction"].astype(float)
        ys[origin+"_species_minus_family_z"]=ys[origin+"_species"]-ys[origin+"_family"]
        ys[origin+"_species_minus_genus_z"]=ys[origin+"_species"]-ys[origin+"_genus"]
        # Explicit nesting; logit of species-covered among family/genus-covered.
        for rank in ("family","genus"):
            n=f[prefix+origin+"_species"]; d=f[prefix+origin+"_"+rank]
            ys[origin+"_species_given_"+rank+"_logodds"]=np.log((n+.5)/(d-n+.5))
    for rank in ("species","genus","family"):
        ys[rank+"_current_minus_native_z"]=ys["current_"+rank]-ys["native_"+rank]
    out={}
    for name,y in ys.items():
        out[name]={"base":term(y,X,f.level1),"MEM10":term(y,Xmem,f.level1)}
    return out,ys

def main():
    p=argparse.ArgumentParser()
    for arg in ("edges","native","current","taxonomy","metrics","port-flora","region-hosts","tdwg-geojson","out-dir"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args(); a.out_dir.mkdir(parents=True,exist_ok=True)
    for k,sha in EXPECTED.items():
        require(digest(getattr(a,k))==sha,"Original source SHA mismatch: "+k)
    metrics=read(a.metrics)
    panel=set(metrics.species)
    require(len(metrics)==len(panel)==239,"Wrong frozen panel")
    taxdf=read(a.taxonomy)
    require(taxdf.accepted_plant_name_id.is_unique,"Ambiguous accepted taxon IDs")
    tax=taxdf.set_index("accepted_plant_name_id")[["genus","family"]].to_dict("index")
    n,c=read(a.native),read(a.current)
    require(len(n)==1202784 and len(c)==1338575,"Original allflora source row counts changed")
    for d in (n,c):
        require(not d.duplicated(["accepted_plant_name_id","area_code_l3"]).any(),"Duplicate plant-region")
        require(set(d.accepted_plant_name_id)<=set(tax),"Missing accepted taxon mapping")
    edges=read(a.edges)
    edges=edges[edges.insect_species.isin(panel) & edges.family.ne("")]
    e=edges[["insect_species","accepted_plant_name_id"]].drop_duplicates()
    require(len(e)==2596 and e.accepted_plant_name_id.nunique()==1706 and set(e.insect_species)==panel,
            "Frozen accepted interaction edges changed")
    require(set(e.accepted_plant_name_id)<=set(tax),"Host is outside accepted species taxonomy")
    hosts={}
    for b,g in e.groupby("insect_species"):
        ids=set(g.accepted_plant_name_id)
        hosts[b]={"species":ids,**ranks(ids,tax)}
    no_big={b:h for b,h in hosts.items() if not h["family"] & {"Fabaceae","Poaceae"}}
    require(len(no_big)>0,"No non-Fabaceae/Poaceae butterflies")
    by_n=n.groupby("area_code_l3").accepted_plant_name_id.agg(set).to_dict()
    by_c=c.groupby("area_code_l3").accepted_plant_name_id.agg(set).to_dict()
    single=set(n.groupby("accepted_plant_name_id").area_code_l3.nunique().loc[lambda s:s==1].index)
    f=make_frame(a.port_flora,a.region_hosts,a.tdwg_geojson).reset_index(drop=True)
    require(len(f)==126 and f.code.is_unique and f.level1.nunique()==9,"Wrong PR188 cohort")
    require(set(f.code)<=set(by_n) and set(f.code)<=set(by_c),"Missing botanical region")
    used=set().union(*[by_c[k] for k in f.code],*list(h["species"] for h in hosts.values()))
    require(all(tax[k]["family"] and tax[k]["genus"] for k in used),"Missing family/genus in used taxa")
    records=[]
    for row in f.itertuples(index=False):
        ns,cs=by_n[row.code],by_c[row.code]
        require(ns<=cs,"Native not nested in contemporary")
        rec={"code":row.code,"single_native_region_share":len(ns&single)/len(ns),
             "allflora_native_n":len(ns),"allflora_current_n":len(cs)}
        for prefix,cohort in (("",hosts),("no_big_",no_big)):
            for origin,ids in (("native",ns),("current",cs)):
                vals=aggregate(cohort,ids,tax)
                rec.update({prefix+origin+"_"+k:v for k,v in vals.items()})
        records.append(rec)
    f=f.merge(pd.DataFrame(records),on="code",validate="one_to_one",sort=False)
    f["native_species_delta"]=f.native_species-f.native_host_239
    f["current_species_delta"]=f.current_species-f.current_host_239
    f.to_csv(a.out_dir/"island_taxonomic_resolution_counts.csv",index=False)
    mismatch=f.loc[(f.native_species_delta!=0)|(f.current_species_delta!=0),
                   ["code","native_species","native_host_239","current_species","current_host_239","native_species_delta","current_species_delta"]]
    receipt={"status":"SOURCE_CHECK_ONLY","n_islands":len(f),"n_butterflies":len(hosts),
             "n_excluding_any_Fabaceae_or_Poaceae":len(no_big),"n_species_count_mismatches":len(mismatch),
             "native_correlation":float(f.native_species.corr(f.native_host_239)),
             "current_correlation":float(f.current_species.corr(f.current_host_239)),
             "native_species_total":[int(f.native_species.sum()),int(f.native_host_239.sum())],
             "current_species_total":[int(f.current_species.sum()),int(f.current_host_239.sum())],
             "input_sha256":{k:digest(getattr(a,k)) for k in ("edges","native","current","taxonomy","metrics","port_flora","region_hosts","tdwg_geojson")}}
    (a.out_dir/"source_gate.json").write_text(json.dumps(receipt,indent=2)+"\n")
    mismatch.to_csv(a.out_dir/"species_count_mismatch.csv",index=False)
    require(mismatch.empty,"Species opportunities must exactly reproduce PR188 for all 126 regions; see source_gate and mismatch CSV")
    w,graph=make_weights(f,a.tdwg_geojson)
    mem,_=preselected_mem(w)
    X=design(f); Xm=pd.concat([X,pd.DataFrame(mem,index=f.index,columns=[f"MEM{i+1}" for i in range(mem.shape[1])])],axis=1)
    full,ys=suite(f,"",X,Xm); reduced,_=suite(f,"no_big_",X,Xm)
    # A geographic range-size proxy, not endemic status or mediation.
    endog=f.single_native_region_share
    Xend=X.assign(single_native_region_share=endog)
    end={"single_region_share_distance":term(endog,X,f.level1)}
    for rank in ("family","genus"):
        key="native_species_given_"+rank+"_logodds"
        end[key+"_adjusted_for_single_region_share"]=term(ys[key],Xend,f.level1)
    nearfar=f.sort_values(["continent_distance_km","code"])
    k=len(f)//4
    med={}
    for name,g in (("nearest_31",nearfar.head(k)),("farthest_31",nearfar.tail(k))):
        med[name]={k:float(g[k].median()) for k in ("continent_distance_km","native_species","native_genus","native_family","single_native_region_share")}
    report={**receipt,"status":"EXACT_SPECIES_COUNT_GATE_PASSED_EXPLORATORY_TAXONOMIC_AUDIT",
            "graph":graph,"full":full,"no_Fabaceae_or_Poaceae_butterflies":reduced,
            "single_native_region_proxy":end,"distance_quartile_medians":med,
            "limits":["Taxonomic co-occurrence is not actual host suitability.",
                      "Single recorded native WGSRPD3 region is a range-size proxy, not verified island endemism.",
                      "Nested family/genus coverage can saturate; fractional coverage and paired contrasts included.",
                      "All 239 butterflies fixed globally, not an observed island assemblage.",
                      "Unknown host links, ecological incompatibility and documentation coverage cannot be separated.",
                      "HC3, nine-cluster t8 and MEM10 are separate descriptive sensitivity checks; exploratory multiple endpoints."]}
    (a.out_dir/"taxonomic_resolution_report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":report["status"],"quartile_medians":med,
                      "full_native_counts":{k:full["native_"+k] for k in ("species","genus","family")}},indent=2))
if __name__=="__main__":
    main()

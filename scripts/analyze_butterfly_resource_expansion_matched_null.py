#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

def ranks(v):
    x=np.asarray(v,float); o=np.argsort(x,kind="mergesort"); r=np.empty(len(x)); i=0
    while i<len(o):
        j=i+1
        while j<len(o) and x[o[j]]==x[o[i]]: j+=1
        r[o[i:j]]=0.5*((i+1)+j); i=j
    return r

def spear(a,b):
    if len(a)<3 or len(a)!=len(b): return None
    x,y=ranks(a),ranks(b)
    if np.std(x)==0 or np.std(y)==0:return None
    return float(np.corrcoef(x,y)[0,1])

def load_desc(p):
    out={}
    with p.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[r["species"]]={"host_family_count":float(r["host_family_count"]),"resolved_host_species":int(r["resolved_host_species"]),"native":int(r["host_wgsrpd3_unit_count"])}
    return out

def load_pairs(p):
    by=defaultdict(dict); pools=defaultdict(set)
    with p.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f); req={"insect_species","accepted_plant_name_id","family"}
        if not req<=set(rd.fieldnames or ()): raise RuntimeError("interaction schema drift")
        for r in rd:
            s,h,fam=r["insect_species"].strip(),r["accepted_plant_name_id"].strip(),r["family"].strip()
            if s and h and fam: by[s][h]=fam; pools[fam].add(h)
    return {s:dict(v) for s,v in by.items()},{k:tuple(sorted(v)) for k,v in pools.items()}

def load_units(p):
    d=defaultdict(set)
    with p.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f); req={"accepted_plant_name_id","area_code_l3"}
        if not req<=set(rd.fieldnames or ()): raise RuntimeError("distribution schema drift")
        for r in rd:
            h,u=r["accepted_plant_name_id"].strip(),r["area_code_l3"].strip()
            if h and u:d[h].add(u)
    return d

def metrics(hosts,native,contemp):
    nu=set(); cu=set(); added={}
    for h in hosts:
        n=set(native.get(h,set())); c=set(contemp.get(h,set()))
        if not n<=c: raise RuntimeError(f"native not subset contemporary: {h}")
        nu|=n; cu|=c; added[h]=c-n
    au=cu-nu
    out={"native_units":len(nu),"contemporary_units":len(cu),"introduced_added_units":len(au),
         "log_expansion":float(math.log1p(len(cu))-math.log1p(len(nu)))}
    if not au:
        out.update(max_share=None,effective=None,contributors=0); return out
    credit=defaultdict(float)
    for u in au:
        hs=[h for h in hosts if u in added[h]]
        q=1/len(hs)
        for h in hs:credit[h]+=q
    p=[v/len(au) for v in credit.values()]
    out.update(max_share=float(max(p)),effective=float(1/sum(x*x for x in p)),contributors=len(p))
    return out

def pick(pool,n,seed):
    if n>len(pool): raise RuntimeError("family pool smaller than focal family count")
    return random.Random(seed).sample(list(pool),n)

def seed_int(tag,species,rep,fam):
    b=f"{tag}|{species}|{rep}|{fam}".encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8],"big")

def qtail(vals,obs,side):
    vals=[x for x in vals if x is not None and math.isfinite(x)]
    if not vals or obs is None:return None
    k=sum(x>=obs for x in vals) if side=="high" else sum(x<=obs for x in vals)
    return (k+1)/(len(vals)+1)

def main():
    ap=argparse.ArgumentParser()
    for x in ["protocol_json","descriptors_csv","insect_host_csv","native_distribution_csv","contemporary_distribution_csv","output_json","output_csv"]:
        ap.add_argument("--"+x.replace("_","-"),type=Path,required=True)
    a=ap.parse_args(); prot=json.loads(a.protocol_json.read_text()); null=prot["null"]; B=int(null["permutations"]); tag=str(null["seed_tag"])
    desc=load_desc(a.descriptors_csv); pairs,pools=load_pairs(a.insect_host_csv); native=load_units(a.native_distribution_csv); contemp=load_units(a.contemporary_distribution_csv)
    pools={fam:tuple(h for h in hs if native.get(h) and native[h] <= contemp.get(h,set())) for fam,hs in pools.items()}
    focal=[]
    for s,d in desc.items():
        if d["host_family_count"]<=0 or d["native"]<=0 or d["resolved_host_species"]<d["host_family_count"]:continue
        hs=pairs.get(s,{})
        if not hs:continue
        obs=metrics(tuple(sorted(hs)),native,contemp)
        if obs["introduced_added_units"]<=0:continue
        comp=Counter(hs.values())
        if any(len(pools.get(fam,()))<n for fam,n in comp.items()):continue
        focal.append((s,d,hs,comp,obs))
    if not focal:raise RuntimeError("no eligible focal species")
    per={s:{"log_expansion":[],"max_share":[],"effective":[]} for s,*_ in focal}
    null_rho_eff=[]; null_rho_max=[]; null_rho_exp=[]
    family=[d["host_family_count"] for _,d,_,_,_ in focal]
    for rep in range(B):
        eff=[]; mx=[]; ex=[]
        for s,d,hs,comp,obs in focal:
            sampled=[]
            for fam,n in sorted(comp.items()): sampled.extend(pick(pools[fam],n,seed_int(tag,s,rep,fam)))
            m=metrics(tuple(sampled),native,contemp)
            per[s]["log_expansion"].append(m["log_expansion"]); per[s]["max_share"].append(m["max_share"]); per[s]["effective"].append(m["effective"])
            ex.append(m["log_expansion"]); eff.append(m["effective"] if m["effective"] is not None else 0.0); mx.append(m["max_share"] if m["max_share"] is not None else 1.0)
        null_rho_exp.append(spear(family,ex)); null_rho_eff.append(spear(family,eff)); null_rho_max.append(spear(family,mx))
    rows=[]
    for s,d,hs,comp,obs in focal:
        z=per[s]
        rows.append({"species":s,"host_family_count":d["host_family_count"],"resolved_host_species":d["resolved_host_species"],
         "observed_log_expansion":obs["log_expansion"],"null_log_expansion_median":float(np.median(z["log_expansion"])),"p_expansion_high":qtail(z["log_expansion"],obs["log_expansion"],"high"),
         "observed_effective":obs["effective"],"null_effective_median":float(np.nanmedian([np.nan if x is None else x for x in z["effective"]])),"p_effective_high":qtail(z["effective"],obs["effective"],"high"),
         "observed_max_share":obs["max_share"],"null_max_share_median":float(np.nanmedian([np.nan if x is None else x for x in z["max_share"]])),"p_max_share_high":qtail(z["max_share"],obs["max_share"],"high"),
         "family_composition":";".join(f"{k}:{v}" for k,v in sorted(comp.items()))})
    obs_eff=spear(family,[x[4]["effective"] for x in focal]); obs_max=spear(family,[x[4]["max_share"] for x in focal]); obs_exp=spear(family,[x[4]["log_expansion"] for x in focal])
    payload={"schema":"chocho_butterfly_family_composition_matched_null_v0.1","status":"POSTHOC_REVIEWER_DEFENSE_SENSITIVITY",
      "permutations":B,"species":len(focal),"null_definition":"For each butterfly, preserve exact WCVP host-family composition and resolved host-species count; sample alternative HOSTS-WCVP hosts without replacement within each family; reconstruct native and contemporary WGSRPD3 unions.",
      "observed":{"rho_family_vs_log_expansion":obs_exp,"rho_family_vs_effective_contributors":obs_eff,"rho_family_vs_max_share":obs_max},
      "null":{"rho_family_vs_log_expansion":{"median":float(np.nanmedian(null_rho_exp)),"p_two_sided":(1+sum(abs(x)>=abs(obs_exp) for x in null_rho_exp))/(B+1)},
              "rho_family_vs_effective_contributors":{"median":float(np.nanmedian(null_rho_eff)),"p_high":qtail(null_rho_eff,obs_eff,"high")},
              "rho_family_vs_max_share":{"median":float(np.nanmedian(null_rho_max)),"p_low":qtail(null_rho_max,obs_max,"low")}},
      "species_with_expansion_above_family_matched_null_p05":sum(r["p_expansion_high"]<=0.05 for r in rows),
      "claim_boundary":"This post-hoc sensitivity controls exact host count and host-family composition; it does not randomize evolutionary host accessibility or prove causal host choice."}
    a.output_json.parent.mkdir(parents=True,exist_ok=True); a.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    with a.output_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(payload,indent=2,sort_keys=True))
if __name__=="__main__":main()

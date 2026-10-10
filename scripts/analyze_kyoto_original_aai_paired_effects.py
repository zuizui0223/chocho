#!/usr/bin/env python3
"""Source-verified exploratory reanalysis of originally PUBLISHED split-leaf AAI.

Effect sizes are diagnostic of precision only, NOT an original biological finding.
Pair resampling preserves common leaf identity and all post-assignment loss codes.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

FILES={
 "larvae":("https://ndownloader.figshare.com/files/41146985","982d4931da306a7ff8e8d550cffe5246"),
 "leafarea":("https://ndownloader.figshare.com/files/41146988","ebe8c2404d9a13b73f08c16ecb9ddc22"),
}
N_BOOT=9999
SEED=20261010
SPECIES={"a":"Atrophaneura alcinous","s":"Sericinus montela"}
METRICS=["growth_ratio","consumption_ratio"]

def fetch(key):
    url,checksum=FILES[key]
    with urlopen(Request(url,headers={"User-Agent":"chocho-public-source-aai-precision-audit/1.0"}),timeout=40) as r:
        raw=r.read(500001)
    if len(raw)>500000 or hashlib.md5(raw).hexdigest()!=checksum:
        raise ValueError("original Figshare AAI source checksum mismatch")
    return raw

def float_exact(value,context):
    try:v=float(value)
    except (TypeError,ValueError):raise ValueError(f"bad numeric in {context}")
    if not math.isfinite(v):raise ValueError("nonfinite recorded original source")
    return v

def parse(larva,leaf):
    a=list(csv.DictReader(io.StringIO(larva.decode("utf-8-sig"))))
    b=list(csv.DictReader(io.StringIO(leaf.decode("utf-8-sig"))))
    if len(a)!=120 or len(b)!=120:raise ValueError("source row count changed")
    A={}
    B={}
    for r in a:
        ident=r["h.id"]
        if ident in A:raise ValueError("larva ID repeated")
        A[ident]=r
    for r in b:
        ident=r["h.id"]
        if ident in B:raise ValueError("leafarea ID repeated")
        B[ident]=r
    if set(A)!=set(B):raise ValueError("larval and leafarea ID sets inconsistent")
    leaves=defaultdict(dict)
    observed=[]
    half_label_mismatch=[]
    recorded_leaf_area_notes=defaultdict(Counter)
    zero_leaf_areas=defaultdict(Counter)
    for ident,row in A.items():
        food=B[ident]
        if row["l.id"]!=food["l.id"] or row["treatment"]!=food["treatment"]:
            raise ValueError("larval and leaf treatment/source identity differ")
        if row["l.or.r"]!=food["l.or.r"]:
            half_label_mismatch.append({
                "h.id":ident, "l.id":row["l.id"],
                "larval_half_label":row["l.or.r"],
                "leaf_area_half_label":food["l.or.r"]
            })
        sp=row["species"]
        t=row["treatment"]
        if sp not in SPECIES or t not in ("a","c"):
            raise ValueError("unexpected butterfly or treatment")
        recorded_leaf_area_notes[sp][food.get("note","").strip()]+=1
        if food.get("consumed.leaf.area","").strip() in ("0","0.0","0.000"):
            zero_leaf_areas[sp][t]+=1
        loss=row["loss"]
        if loss not in ("0","1"):
            raise ValueError("unexpected author loss code")
        if t in leaves[(sp,row["l.id"])]:
            raise ValueError("duplicate original leaf condition")
        result={"id":ident,"species":sp,"treatment":t,"loss":int(loss),
                "leaf":row["l.id"]}
        if loss=="0":
            initial=float_exact(row["initial.mass"],"initial.mass")
            day1=float_exact(row["day1.mass"],"day1.mass")
            consumed=float_exact(food["consumed.leaf.area"],"consumed.leaf.area")
            if initial<=0 or consumed<0:
                raise ValueError("invalid initial mass or feeding area")
            result["growth_ratio"]=(day1-initial)/initial
            result["consumption_ratio"]=consumed/initial
        else:
            result.update({key:None for key in METRICS})
        leaves[(sp,row["l.id"])][t]=result
        observed.append(result)
    for (species,leaf_id),pair in leaves.items():
        if set(pair)!={"a","c"}:
            raise ValueError(f"not original split-leaf AAI-control pair {species} {leaf_id}")
    if len(leaves)!=60:
        raise ValueError("not the expected 60 paired leaves")
    quality={"recorded_leaf_area_notes":{SPECIES[k]:dict(v) for k,v in recorded_leaf_area_notes.items()},
             "exact_zero_consumption_by_treatment":{SPECIES[k]:dict(v) for k,v in zero_leaf_areas.items()}}
    return leaves,observed,half_label_mismatch,quality

def estimate(vals,metric,complete_pairs_only):
    diffs=[]
    aa=[]
    cc=[]
    for pair in vals:
        a=pair["a"][metric];c=pair["c"][metric]
        if complete_pairs_only:
            if a is not None and c is not None:diffs.append(a-c)
        else:
            if a is not None:aa.append(a)
            if c is not None:cc.append(c)
    if complete_pairs_only:
        return sum(diffs)/len(diffs) if diffs else None
    return sum(aa)/len(aa)-sum(cc)/len(cc) if aa and cc else None

def bootstrap(pairs,metric,complete_only,seed,draws=N_BOOT):
    if estimate(pairs,metric,complete_only) is None:return None
    rng=random.Random(seed)
    effects=[]
    for _ in range(draws):
        sample=[pairs[rng.randrange(len(pairs))] for _ in pairs]
        x=estimate(sample,metric,complete_only)
        if x is not None:effects.append(x)
    if len(effects)<int(.95*draws):
        raise ValueError("too few valid original paired resamples")
    effects.sort()
    n=len(effects)
    return [effects[int(.025*n)],effects[int(.975*n)]]

def analyse(original_larva,original_leaf):
    leaves,rows,half_label_mismatch,quality=parse(original_larva,original_leaf)
    bysp=defaultdict(list)
    for (sp,lid),pair in sorted(leaves.items()):
        bysp[sp].append(pair)
    result={"schema":"chocho_kyoto_aai_original_paired_precision_v01",
            "source_doi":"10.6084/m9.figshare.23170898",
            "source_checksums":{key:FILES[key][1] for key in FILES},
            "methods":"Published 24h AAI single-compound split-leaf assay, NOT independent replication",
            "n_original_individuals":len(rows),
            "n_original_leaf_pair_blocks":len(leaves),
            "original_leaf_half_label_discrepancies":half_label_mismatch,
            "original_recorded_leaf_area_quality":quality,
            "species":{},
            "new_biological_effect_established":False,
            "can_rule_out_all_compound_blends":False,
            "can_rule_out_induced_plant_quality":False,
            "GEB_PR38_untouched":True}
    for sp,pp in sorted(bysp.items()):
        loss=Counter((t,int(pair[t]["loss"])) for pair in pp for t in ("a","c"))
        both=sum(pair["a"]["loss"]==0 and pair["c"]["loss"]==0 for pair in pp)
        clean=sum(pair[t]["loss"]==0 for pair in pp for t in ("a","c"))
        metric_summary={}
        for k,metric in enumerate(METRICS):
            metric_summary[metric]={
                "complete_leaf_pair_treated_minus_control":estimate(pp,metric,True),
                "paired_leaf_95ci":bootstrap(pp,metric,True,SEED+17*k+(0 if sp=="a" else 100)),
                "all_observed_treated_minus_control":estimate(pp,metric,False),
                "all_observed_leaf_bootstrap_95ci":bootstrap(pp,metric,False,SEED+17*k+1000+(0 if sp=="a" else 100))}
        result["species"][SPECIES[sp]]={
            "original_leaves":len(pp),"assigned_original_larvae":len(pp)*2,
            "both_treatments_valid_leaves":both,
            "valid_individuals":clean,
            "treatment_counts":{"AAI":len(pp),"control":len(pp)},
            "lost_after_assignment":{"AAI":loss[("a",1)],"control":loss[("c",1)]},
            "outcomes":metric_summary,
            "postrandomization_complete_pair_selection":True,
            "interpretation":"95% interval containing zero does not show equivalence; 24h in 3rd instar does not test induced responses or adult recruitment."
        }
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--receipt",type=Path,required=True);a=p.parse_args()
    r=analyse(fetch("larvae"),fetch("leafarea"))
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps(r,indent=2))
if __name__=="__main__":main()

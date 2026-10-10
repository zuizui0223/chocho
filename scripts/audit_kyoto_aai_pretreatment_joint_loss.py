#!/usr/bin/env python3
"""Constrained exploratory audit of pre-exposure factors in shared AAI assay loss.

30 Atrophaneura source-leaf pairs, source-strain fixed permutation. Does not
discover chemical causality, identify death reason, or retest treatment effect.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

ORIGINAL={
 "larvae":("https://ndownloader.figshare.com/files/41146985","982d4931da306a7ff8e8d550cffe5246"),
 "leafarea":("https://ndownloader.figshare.com/files/41146988","ebe8c2404d9a13b73f08c16ecb9ddc22")
}
SEED=20261010
DRAWS=9999
STATS=("initial_half_leaf_area_mean","initial_larva_logmass_mean","composite_multi_leaf_id")

def download(name):
    url,md5=ORIGINAL[name]
    with urlopen(Request(url,headers={"User-Agent":"chocho-plant-source-identity-audit/1.0"}),timeout=35) as conn:
        raw=conn.read(500001)
    if len(raw)>500000 or hashlib.md5(raw).hexdigest()!=md5:
        raise ValueError("source file byte integrity mismatch")
    return raw

def parse(original_larvae,original_leafarea):
    lar=list(csv.DictReader(io.StringIO(original_larvae.decode("utf-8-sig"))))
    leaf=list(csv.DictReader(io.StringIO(original_leafarea.decode("utf-8-sig"))))
    if len(lar)!=120 or len(leaf)!=120:raise ValueError("original records no longer 120/120")
    leafid={}
    for r in leaf:
        h=r["h.id"]
        if h in leafid:raise ValueError("duplicate leafarea caterpillar ID")
        leafid[h]=r
    groups=defaultdict(dict)
    half_mismatches=[]
    for r in lar:
        if r["species"] not in ("a","s"):raise ValueError("unexpected source species")
        h=r["h.id"]
        v=leafid.get(h)
        if not v or (r["l.id"],r["treatment"])!=(v["l.id"],v["treatment"]):
            raise ValueError("unmatched original leaf and caterpillar record")
        if r["l.or.r"]!=v["l.or.r"]:
            half_mismatches.append({"h.id":h,"source_leaf":r["l.id"],
                                    "larval_sides":r["l.or.r"],"leafarea_sides":v["l.or.r"]})
        key=(r["species"],r["l.id"])
        arm=r["treatment"]
        if arm not in ("a","c") or arm in groups[key]:
            raise ValueError("not a split-leaf original paired treatment")
        if r["loss"] not in ("0","1"):raise ValueError("unknown source loss")
        ia=float(v["initial.leaf.area"])
        mass=float(r["initial.mass"])
        if not math.isfinite(ia) or ia<=0 or not math.isfinite(mass) or mass<=0:
            raise ValueError("invalid PRETREATMENT plant/larval measurements")
        groups[key][arm]={"strain":r["strain"],"loss":int(r["loss"]),
                          "leafarea":ia,"mass":mass}
    native=[]
    source_leaf_counts=Counter()
    source_leaf_composite_count=Counter()
    for (species,lid),pp in groups.items():
        if set(pp)!={"a","c"}:raise ValueError("not two treatment records per source leaf")
        if species!="a":continue
        if pp["a"]["strain"]!=pp["c"]["strain"]:
            raise ValueError("native source leaf halves do not share original strain")
        source_leaf_counts[lid]+=1
        composite=int("," in lid)
        source_leaf_composite_count["composite" if composite else "single"]+=1
        native.append({
            "source_leaf_ID":lid,"strain":pp["a"]["strain"],
            "both_lost":int(pp["a"]["loss"]==1 and pp["c"]["loss"]==1),
            "any_lost":int(pp["a"]["loss"]==1 or pp["c"]["loss"]==1),
            "loss_AAI":pp["a"]["loss"],"loss_control":pp["c"]["loss"],
            "initial_half_leaf_area_mean":(pp["a"]["leafarea"]+pp["c"]["leafarea"])/2,
            "initial_larva_logmass_mean":(math.log(pp["a"]["mass"])+math.log(pp["c"]["mass"]))/2,
            "composite_multi_leaf_id":composite
        })
    if len(native)!=30 or len(source_leaf_counts)!=30 or sum(x["both_lost"] for x in native)!=7:
        raise ValueError("original A.alcinous 30 source leaves with 7 both lost not preserved")
    if len(half_mismatches)!=1 or half_mismatches[0]["h.id"]!="a4":
        raise ValueError("source leaf-side discrepancy changed")
    return native,dict(source_leaf_composite_count),half_mismatches

def contrast(rows,name,flags=None):
    labels=flags if flags is not None else [r["both_lost"] for r in rows]
    yes=[r[name] for r,v in zip(rows,labels) if v]
    no=[r[name] for r,v in zip(rows,labels) if not v]
    if not yes or not no:raise ValueError("no baseline comparator")
    return sum(yes)/len(yes)-sum(no)/len(no)

def randomization_test(rows,seed=SEED,draws=DRAWS):
    strata=defaultdict(list)
    for ix,row in enumerate(rows):strata[row["strain"]].append(ix)
    fixed={s:sum(rows[i]["both_lost"] for i in idx) for s,idx in strata.items()}
    stats={name:contrast(rows,name) for name in STATS}
    grouped_descriptive={
        name:{
            "n_joint_loss":sum(bool(r["both_lost"]) for r in rows),
            "n_other":sum(not bool(r["both_lost"]) for r in rows),
            "joint_loss_mean":sum(r[name] for r in rows if r["both_lost"])/sum(r["both_lost"] for r in rows),
            "other_mean":sum(r[name] for r in rows if not r["both_lost"])/sum(not r["both_lost"] for r in rows)
        }
        for name in STATS
    }
    extreme={name:0 for name in STATS}
    rng=random.Random(seed)
    for _ in range(draws):
        labels=[0]*len(rows)
        for s,idx in sorted(strata.items()):
            chosen=rng.sample(idx,fixed[s])
            for ix in chosen:labels[ix]=1
        for name in STATS:
            delta=contrast(rows,name,labels)
            if abs(delta)>=abs(stats[name])-1e-12:
                extreme[name]+=1
    pvals={name:(extreme[name]+1)/(draws+1) for name in STATS}
    sorted_p=sorted(STATS,key=lambda name:pvals[name])
    adj={};prev=0
    for k,name in enumerate(sorted_p):
        prev=max(prev,min(1.0,(len(STATS)-k)*pvals[name]))
        adj[name]=prev
    return {"original_composite_id_by_joint_loss":{
                "composite_joint_loss":sum(bool(r["composite_multi_leaf_id"]) and bool(r["both_lost"]) for r in rows),
                "composite_other":sum(bool(r["composite_multi_leaf_id"]) and not bool(r["both_lost"]) for r in rows),
                "single_joint_loss":sum(not bool(r["composite_multi_leaf_id"]) and bool(r["both_lost"]) for r in rows),
                "single_other":sum(not bool(r["composite_multi_leaf_id"]) and not bool(r["both_lost"]) for r in rows)},
            "statistics":[{"pretreatment_feature":n,
                            "source_group_descriptive":grouped_descriptive[n],
                            "both_lost_minus_other_source_leaf_mean_difference":stats[n],
                            "conditional_strain_fixed_permutation_two_sided_p":pvals[n],
                            "holms_multiplicity_adjusted_p":adj[n]}
                           for n in STATS],
            "source_strain_counts":{s:{"n_pairs":len(idx),"both_lost":fixed[s]}
                                    for s,idx in sorted(strata.items())},
            "draws":draws,"seed":seed,
            "exploratory_not_independent_pre_registration":True,
            "plant_chemical_or_handling_cause_identified":False}

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--receipt",type=Path,required=True)
    a=arg.parse_args()
    larvae,leafarea=download("larvae"),download("leafarea")
    rows,leafids,mismatch=parse(larvae,leafarea)
    obj={"schema":"chocho_AAI_pretreatment_joint_loss_sensitivity_v01",
         "source_url":{k:v[0] for k,v in ORIGINAL.items()},
         "source_md5":{k:v[1] for k,v in ORIGINAL.items()},
         "source_sha256":{"larvae":hashlib.sha256(larvae).hexdigest(),
                          "leafarea":hashlib.sha256(leafarea).hexdigest()},
         "original_native_leaf_pairs":len(rows),
         "observed_joint_loss_source_leaves":sum(v["both_lost"] for v in rows),
         "original_composite_leaf_id_counts":leafids,
         "original_leaf_side_discrepancy":mismatch,
         "source_pretreatment_factors":randomization_test(rows),
         "GEB_PR38_untouched":True}
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(obj,indent=2,ensure_ascii=False))

if __name__=="__main__":main()

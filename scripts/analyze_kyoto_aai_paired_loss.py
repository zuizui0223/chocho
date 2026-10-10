#!/usr/bin/env python3
"""Source-verified pair-matched loss analysis of an ALREADY PUBLISHED AAI assay.

Original haphazard larvae deaths are not independent evidence of AAI toxicity.
This audits ITT source-leaf pairing; it does not fit other outcome models.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import random
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

SOURCE="https://ndownloader.figshare.com/files/41146985"
SOURCE_MD5="982d4931da306a7ff8e8d550cffe5246"
SEED=20261010
BOOT=9999
SPECIES={"a":"Atrophaneura alcinous","s":"Sericinus montela"}

def authentic_data():
    req=Request(SOURCE,headers={"User-Agent":"chocho-academic-data-aai-loss-audit/1.0"})
    with urlopen(req,timeout=35) as r:
        raw=r.read(500001)
    if len(raw)>500000 or hashlib.md5(raw).hexdigest()!=SOURCE_MD5:
        raise ValueError("author original source checksum failed")
    return raw

def load_pairs(raw):
    records=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(records)!=120:
        raise ValueError("unexpected author bioassay cohort length")
    pairs=defaultdict(dict)
    seen_larva=set()
    for r in records:
        h=r["h.id"]
        if h in seen_larva or not h:
            raise ValueError("duplicate/empty experimental caterpillar")
        seen_larva.add(h)
        species=r["species"]
        leaf=r["l.id"]
        trt=r["treatment"]
        loss=r["loss"]
        if species not in SPECIES or not leaf or trt not in {"a","c"} or loss not in {"0","1"}:
            raise ValueError("author treatment/leaf/fate unknown")
        group=pairs[(species,leaf)]
        if trt in group:
            raise ValueError("duplicate treatment within original source leaf")
        group[trt]=int(loss)
    by_species=defaultdict(list)
    for (species,leaf),pair in pairs.items():
        if set(pair)!={"a","c"}:raise ValueError("original source leaf not paired")
        by_species[species].append((pair["a"],pair["c"]))
    if len(by_species["a"])!=30 or len(by_species["s"])!=30:
        raise ValueError("original 30 source leaf pairs for each butterfly species missing")
    return dict(by_species)

def exact_mcnemar(n_a_only:int,n_c_only:int):
    n=n_a_only+n_c_only
    if n==0:return 1.0
    k=min(n_a_only,n_c_only)
    numerator=sum(math.comb(n,j) for j in range(k+1))
    return min(1.0,2.0*numerator/(2**n))

def bootstrap_paired_risk_difference(pairs,seed,draws=BOOT):
    rng=random.Random(seed)
    n=len(pairs)
    d=[a-c for a,c in pairs]
    res=[]
    for _ in range(draws):
        res.append(sum(d[rng.randrange(n)] for _ in range(n))/n)
    res.sort()
    return [res[int(.025*draws)],res[int(.975*draws)]]

def summarize_pair_data(raw):
    cohort=load_pairs(raw)
    output={"schema":"chocho_kyoto_paired_AAI_author_loss_v01",
            "original_source":"Hashimoto & Ohgushi 2023 Figshare file 41146985",
            "source_md5":SOURCE_MD5,
            "source_sha256":hashlib.sha256(raw).hexdigest(),
            "original_butterflies":sum(2*len(v) for v in cohort.values()),
            "species":{},
            "source_loss_is_haphazard_death_not_proven_AAI_toxicity":True,
            "new_biological_mechanism_estimated":False,
            "no_biological_equivalence_margin_predeclared":True,
            "GEB_PR38_untouched":True}
    for index,species in enumerate(("a","s")):
        pairs=cohort[species]
        b=Counter((a,c) for a,c in pairs)
        both=b[(1,1)]
        aa=b[(1,0)]
        cc=b[(0,1)]
        neither=b[(0,0)]
        n=len(pairs)
        loss_a=both+aa
        loss_c=both+cc
        output["species"][SPECIES[species]]={
            "original_leaf_pairs":n,
            "original_larvae":2*n,
            "pair_status":{
                "neither_lost":neither,"only_AAI_lost":aa,
                "only_control_lost":cc,"both_lost":both},
            "loss_AAI":loss_a,
            "loss_control":loss_c,
            "risk_loss_AAI":loss_a/n,
            "risk_loss_control":loss_c/n,
            "paired_loss_risk_difference_AAI_minus_control":(loss_a-loss_c)/n,
            "paired_leaf_bootstrap_95ci":bootstrap_paired_risk_difference(pairs,SEED+index*1701),
            "discordant_original_leaves":aa+cc,
            "paired_exact_mcnemar_two_sided_p":exact_mcnemar(aa,cc),
            "biological_cause_of_loss_identified":False,
            "interpretation":"Loss imbalance and uncertainty require reporting; no significant McNemar signal is proof of no toxicity or no other host-quality response."}
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    result=summarize_pair_data(authentic_data())
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Authenticate original 2023 split-leaf aristolochic acid I assay.

This inventories pairing/missingness BEFORE any outcome comparison. Published
null effects are prior art, and absence of significance is not equivalence.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import urlopen,Request
FILES={
    "larvae":("https://ndownloader.figshare.com/files/41146985","982d4931da306a7ff8e8d550cffe5246"),
    "leafarea":("https://ndownloader.figshare.com/files/41146988","ebe8c2404d9a13b73f08c16ecb9ddc22"),
}
def fetch(name):
    url,md5=FILES[name]
    with urlopen(Request(url,headers={"User-Agent":"chocho-butterfly-paired-leaf-scientific-source-audit/1.0"}),timeout=40) as r:
        data=r.read(500001)
    if len(data)>500000 or hashlib.md5(data).hexdigest()!=md5:
        raise ValueError("original source checksum mismatch "+name)
    return data
def getrows(raw):
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
def summarize(larva,leaf):
    a=getrows(larva);b=getrows(leaf)
    if len(a)!=120 or len(b)!=120:raise ValueError("original 120+120 assay row counts changed")
    for rows,req in [(a,{"h.id","strain","species","initial.mass","day1.mass","l.id","l.or.r","treatment","loss"}),
                     (b,{"l.id","l.or.r","treatment","initial.leaf.area","h.id","day1.leaf.area","consumed.leaf.area","note"})]:
        if not req.issubset(rows[0]):raise ValueError("source header mismatch")
    def table(rows):
        return {"header":list(rows[0]),"n":len(rows),
                "unique_hid":len({r["h.id"] for r in rows}),
                "unique_lid":len({r["l.id"] for r in rows}),
                "species_counts":dict(Counter(r.get("species","") for r in rows)) if "species" in rows[0] else {},
                "treatment":dict(Counter(r["treatment"] for r in rows)),
                "leaf_side":dict(Counter(r["l.or.r"] for r in rows)),
                "h_id_replicates":dict(Counter(Counter(r["h.id"] for r in rows).values())),
                "l_id_replicates":dict(Counter(Counter(r["l.id"] for r in rows).values())),
                "missing":{k:n for k in rows[0] if (n:=sum(not str(r.get(k) or "").strip() for r in rows))},
                "head_three":rows[:3]}
    t1=table(a);t2=table(b)
    larva_keys=[(r["h.id"],r["l.id"],r["l.or.r"],r["treatment"]) for r in a]
    leaf_keys=[(r["h.id"],r["l.id"],r["l.or.r"],r["treatment"]) for r in b]
    if len(set(larva_keys))!=len(a) or len(set(leaf_keys))!=len(b):
        raise ValueError("row compound IDs not unique")
    unmatched_larva=[r for r in a if (r["h.id"],r["l.id"],r["l.or.r"],r["treatment"]) not in set(leaf_keys)]
    unmatched_leaf=[r for r in b if (r["h.id"],r["l.id"],r["l.or.r"],r["treatment"]) not in set(larva_keys)]
    by_hid_a={r["h.id"]:r for r in a}
    by_hid_b={r["h.id"]:r for r in b}
    if set(by_hid_a)!=set(by_hid_b):
        raise ValueError("original larval and leaf identifiers do not align one-to-one")
    species_by_leaf=defaultdict(set)
    arms_by_leaf=defaultdict(set)
    for r in a:
        species_by_leaf[r["l.id"]].add(r["species"])
        arms_by_leaf[r["l.id"]].add(r["treatment"])
    losses=Counter((r["species"],r["treatment"],r["loss"]) for r in a)
    species_treatment=Counter((r["species"],r["treatment"]) for r in a)
    leaf_design={"leaf_ids":len(species_by_leaf),
                 "two_treatments_per_leaf":sum(arms=={"a","c"} for arms in arms_by_leaf.values()),
                 "one_species_per_leaf":sum(len(sps)==1 for sps in species_by_leaf.values()),
                 "multiple_species_per_leaf":sum(len(sps)>1 for sps in species_by_leaf.values())}
    return {"schema":"chocho_kyoto_aai_paired_leaf_source_inventory_v01",
        "source":"Hashimoto and Ohgushi 2023, Figshare 23170898, not independent experiment",
        "source_sha256":{"larvae":hashlib.sha256(larva).hexdigest(),"leafarea":hashlib.sha256(leaf).hexdigest()},
        "source_md5":{name:v[1] for name,v in FILES.items()},
        "larvae":t1,"leafarea":t2,
        "matching_row_keys":len(set(larva_keys)&set(leaf_keys)),
        "unmatched_larva_keys":unmatched_larva,
        "unmatched_leaf_keys":unmatched_leaf,
        "leaf_paired_design":leaf_design,
        "species_treatment_counts":{"|".join(k):v for k,v in species_treatment.items()},
        "species_treatment_loss_codes":{"|".join(k):v for k,v in losses.items()},
        "source_effect_estimated":False,"can_reject_all_plant_quality_mechanisms":False,
        "GEB_PR38_unchanged":True}
def main():
    p=argparse.ArgumentParser();p.add_argument("--receipt",type=Path,required=True);a=p.parse_args()
    result=summarize(fetch("larvae"),fetch("leafarea"))
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":main()

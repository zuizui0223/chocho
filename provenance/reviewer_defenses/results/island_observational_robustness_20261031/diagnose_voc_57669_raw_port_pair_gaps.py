#!/usr/bin/env python3
"""Diagnose genuine 57,669 VOC intra-Asian origin/destination coverage.

Uses independent source authenticity check, and PRE-OUTCOME historical place
alias lock. No source ship identity/record key, raw dates or individual voyages
are published. A zero observed city-to-city island-pair does not prove
historical absence. Counts distinguish same-place, one-end, unmapped and
two-distinct-WPI51-island source rows for both strict and broad aliases.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter
from pathlib import Path
from derive_dataversenl_17c_VOC_island_routes import (
    SOURCE, LOCK, EXPECTED_SHA1, EXPECTED_ROWS, PORT_FROM, PORT_TO,
    DATE_FROM, DATE_TO, SOURCE_KEY, digest, source_year, normalize_port,
    load_aliases)
BASE=Path(__file__).resolve().parent
OUT=BASE/"voc_57669_17c_portpair_diagnostic_v1.json"


def scan(source=SOURCE, lock=LOCK):
    if digest(source,"sha1")!=EXPECTED_SHA1:
        raise RuntimeError("historical 57669 source SHA1 does not match Dataverse")
    strict, strict_only, broad=load_aliases(lock,False)
    extra,_,_=load_aliases(lock,True)
    counts=Counter()
    raw_years=Counter()
    route_all=Counter()
    route_offdiagonal=Counter()
    unmapped_near_known_from=Counter()
    unmapped_near_known_to=Counter()
    strict_ends=Counter()
    broad_ends=Counter()
    strict_links=Counter()
    broad_links=Counter()
    both_strict_same_island=Counter()
    no_year_links=Counter()
    departure_not_valid=Counter()
    arrival_not_valid=Counter()
    for r in csv.DictReader(source.open(encoding="utf-8-sig",newline=""),delimiter="\t"):
        counts["source_rows"]+=1
        a,b=(normalize_port(r.get(k)) for k in (PORT_FROM,PORT_TO))
        if not a or not b:
            counts["blank_endpoint"]+=1
        if a==b:
            counts["same_original_from_and_to_place"]+=1
        else:
            counts["different_original_from_and_to_place"]+=1
        yr1=source_year(r.get(DATE_FROM))
        yr2=source_year(r.get(DATE_TO))
        if yr1 is None:
            departure_not_valid[str(r.get(DATE_FROM) or "")[:20]]+=1
        else:
            raw_years[yr1]+=1
        if yr2 is None:
            arrival_not_valid[str(r.get(DATE_TO) or "")[:20]]+=1
        if a and b:
            route_all[(a,b)]+=1
            if a!=b:
                route_offdiagonal[(a,b)]+=1
        u,v=strict.get(a),strict.get(b)
        A,B=extra.get(a),extra.get(b)
        counts["one_or_more_strict_anchor_endpoint"]+=int(bool(u or v))
        counts["both_strict_endpoints"]+=int(bool(u and v))
        counts["two_distinct_strict_island_endpoints"]+=int(bool(u and v and u!=v))
        counts["one_or_more_broad_anchor_endpoint"]+=int(bool(A or B))
        counts["both_broad_endpoints"]+=int(bool(A and B))
        counts["two_distinct_broad_island_endpoints"]+=int(bool(A and B and A!=B))
        if u: strict_ends[("from",u,a)]+=1
        if v: strict_ends[("to",v,b)]+=1
        if A: broad_ends[("from",A,a)]+=1
        if B: broad_ends[("to",B,b)]+=1
        if u and not v and a!=b and b:
            unmapped_near_known_from[(a,b)]+=1
        if v and not u and a!=b and a:
            unmapped_near_known_to[(a,b)]+=1
        if u and v and u==v:
            both_strict_same_island[(a,b)]+=1
        if u and v and u!=v:
            strict_links[(a,b,u,v)]+=1
        if A and B and A!=B:
            broad_links[(a,b,A,B)]+=1
        if (u and v and u!=v or A and B and A!=B) and yr1 is None:
            no_year_links[(a,b)]+=1
    if counts["source_rows"]!=EXPECTED_ROWS:
        raise RuntimeError("historical 57,669 source row count mismatch")
    def placepairs(c,max_n=40,min_n=2):
        return [{"from":a,"to":b,"source_rows":n} for (a,b),n in c.most_common(max_n)
                if n>=min_n]
    def anchorpairs(c,max_n=50):
        return [{"from":a,"to":b,"source_rows":n} for (a,b),n in c.most_common(max_n) if n>=2]
    def full_links(c):
        return [{"from":a,"to":b,"island_from":u,"island_to":v,"rows":n}
                for (a,b,u,v),n in c.most_common()]
    return {
        "status":"SOURCE_VERIFIED_57669_VOC_PORT_PAIR_DIAGNOSTIC_NO_BIOLOGICAL_ROUTE_EFFECT",
        "source":{"doi":"10.34894/5HGTCD","sha1":EXPECTED_SHA1,
                  "file_bytes":source.stat().st_size,
                  "response_blind_alias_sha256":digest(lock)},
        "source_row_classification":dict(counts),
        "source_years":{"min":min(raw_years,default=None),
                        "max":max(raw_years,default=None),
                        "n_covered":sum(raw_years.values()),
                        "before_1650":sum(v for k,v in raw_years.items() if k<1650),
                        "1650_1699":sum(v for k,v in raw_years.items() if 1650<=k<1700),
                        "1700_or_later":sum(v for k,v in raw_years.items() if k>=1700)},
        "source_departure_bad_date_values_top10":departure_not_valid.most_common(10),
        "source_arrival_bad_date_values_top10":arrival_not_valid.most_common(10),
        "unique_original_port_pairs":len(route_all),
        "unique_original_offdiagonal_port_pairs":len(route_offdiagonal),
        "common_raw_offdiagonal_anchored_or_not":placepairs(route_offdiagonal),
        "common_original_same_port_segments":placepairs(
            Counter({k:n for k,n in route_all.items() if k[0]==k[1]}),20),
        "strict_port_alias_count":len(strict_only),
        "broad_port_alias_count":len(extra),
        "strict_matched_source_endpoint_summary":[
            {"direction":direction,"island":island,"source_historical_place":name,"rows":n}
            for (direction,island,name),n in strict_ends.most_common()],
        "source_routes_from_known_island_to_unmapped_place":anchorpairs(unmapped_near_known_from),
        "source_routes_from_unmapped_place_to_known_island":anchorpairs(unmapped_near_known_to),
        "strict_both_endpoints_same_island_top25":placepairs(both_strict_same_island,25,1),
        "strict_direct_island_links":full_links(strict_links),
        "broad_direct_island_links":full_links(broad_links),
        "date_uninterpretable_link_rows":dict(no_year_links),
        "limitations":[
            "No observed precise-city pair is not a confirmed absence of all historic trade.",
            "Some geographic names are broad island/region labels and may not designate an individual harbor.",
            "Source port labels are not evidence that modern 2019 WPI point harbors had the same identity.",
            "No source cargo species/plant seed dispersal or arrival of any butterfly was examined.",
            "Source file rows are not necessarily distinct ships or independent voyages; no iid pair p-values.",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--out",type=Path,default=OUT)
    args=ap.parse_args()
    r=scan(args.source)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in r.items() if k not in (
        "source_routes_from_known_island_to_unmapped_place",
        "source_routes_from_unmapped_place_to_known_island",
        "common_raw_offdiagonal_anchored_or_not",
        "strict_matched_source_endpoint_summary",
    )},ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

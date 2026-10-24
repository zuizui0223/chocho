#!/usr/bin/env python3
"""Posthoc weighting stability in the 26-country-cell degree-constrained GBIF host reference.

All alternatives enumerate exactly the SAME simple butterfly×host graphs, keeping
each country, butterfly degree, plant degree and GADM1 records fixed. Only the
reporting unit weights change. A shuffled edge may be a NONHOST: reference tails
are NOT biological significance tests, recruitment or host-use evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

SCHEMA = "chocho_all26_country_degree_host_identity_weighting_posthoc_v1"
METRICS = ("record", "case", "species", "country", "case_n2", "case_n5", "case_n10")
MAX_STATES = 200000
MAX_SUPPORT = 700000

def graph_vectors(butterflies, plants):
    b = sorted(butterflies, key=lambda r: r["butterfly"])
    p = sorted(plants, key=lambda r: r["accepted_id"])
    rd = [r["candidate_k"] for r in b]
    remain = [len(z["links"]) for z in p]
    if sum(rd) != sum(remain) or any(n > len(b) for n in remain) or any(n > len(p) for n in rd):
        raise ValueError("Country host graph degree margins changed")
    observed = tuple(r["observed_count"] for r in b)
    if any(r["N"] < 1 or not 0 <= r["observed_count"] <= r["N"] for r in b):
        raise ValueError("Invalid original country butterfly counts")
    vector_counts = Counter()
    chosen = []
    def dfs(i):
        if i == len(rd):
            if any(remain):
                return
            row = []
            for target, ids in zip(b, chosen):
                footprint = set().union(*(p[j]["gadm"] for j in ids))
                row.append(sum(n for area, n in target["gadm"].items() if area in footprint))
            vector_counts[tuple(row)] += 1
            if sum(vector_counts.values()) > MAX_STATES:
                raise RuntimeError("Exact country enumeration exceeds predeclared computational bound")
            return
        possible = [j for j, k in enumerate(remain) if k]
        for ids in combinations(possible, rd[i]):
            for j in ids:
                remain[j] -= 1
            if sum(remain) == sum(rd[i+1:]) and all(k <= len(rd)-i-1 for k in remain):
                chosen.append(ids)
                dfs(i+1)
                chosen.pop()
            for j in ids:
                remain[j] += 1
    dfs(0)
    if not vector_counts or not vector_counts[observed]:
        raise ValueError("Observed original host graph not in enumerated reference")
    return {"butterflies": b, "vectors": vector_counts, "observed": observed,
            "states": sum(vector_counts.values())}

def factor(row, mode, repeat):
    if mode == "record":
        return Fraction(1)
    if mode == "case":
        return Fraction(1, row["N"])
    if mode == "species":
        return Fraction(1, row["N"] * repeat[row["butterfly"]])
    if mode == "country":
        return Fraction(1, row["N"] * row["country_case_count"])
    if mode.startswith("case_n"):
        return Fraction(int(row["N"] >= int(mode[6:])), row["N"])
    raise ValueError("Unknown weighting target")

def per_country(state, metric, repeat):
    b = state["butterflies"]
    b_with_count = [{**r, "country_case_count": len(b)} for r in b]
    f = [factor(r, metric, repeat) for r in b_with_count]
    scores = Counter()
    for vec, count in state["vectors"].items():
        score = sum((n * w for n, w in zip(vec, f)), Fraction(0))
        scores[score] += count
    if sum(scores.values()) != state["states"]:
        raise AssertionError("A country lost probability mass")
    obs = sum((n * w for n, w in zip(state["observed"], f)), Fraction(0))
    E = sum((x * n for x, n in scores.items()), Fraction(0)) / state["states"]
    return {"counts": scores, "states": state["states"], "observed": obs, "expected": E}

def exact_reference(groups, divisor):
    scale = math.lcm(*(value.denominator for d in groups.values() for value in d["counts"]))
    scores = {0:1}
    for _, g in sorted(groups.items(), key=lambda pair:(len(pair[1]["counts"]),pair[0])):
        updated = defaultdict(int)
        for x, n in g["counts"].items():
            scaled = x * scale
            if scaled.denominator != 1:
                raise AssertionError("Nonintegral reference score")
            shift = int(scaled)
            for prev, count in scores.items():
                updated[prev + shift] += count * n
        scores = updated
        if len(scores) > MAX_SUPPORT:
            raise RuntimeError("Exact conditional support too large; refuse approximation")
    mass = sum(scores.values())
    if mass != math.prod(r["states"] for r in groups.values()):
        raise ValueError("Exact reference graph count changed")
    observed = sum((v["observed"] for v in groups.values()), Fraction(0))
    expected = sum((v["expected"] for v in groups.values()), Fraction(0))
    observed_scaled = observed * scale
    if observed_scaled.denominator != 1:
        raise AssertionError("Observed statistic not an integer scaled support point")
    tail = Fraction(sum(n for x,n in scores.items() if x >= int(observed_scaled)), mass)
    return {
        "observed_fraction":round(float(observed / divisor),10),
        "reference_expected_fraction":round(float(expected / divisor),10),
        "observed_minus_reference_pp":round(100*float((observed-expected)/divisor),7),
        "exact_label_upper_tail_NOT_biological_p":round(float(tail),12),
        "source_unit_denominator":divisor,
        "graph_states":str(mass),
        "distinct_exact_reference_scores":len(scores),
        "countries_where_statistic_varies":sum(len(g["counts"])>1 for g in groups.values()),
    }

def analyse(botanical, butterflies):
    # Import source-verified original checker only at execution time: offline
    # tests load this standalone script by file location.
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    from analyze_gbif_all26_host_identity_reference_20261023 import build_countries
    b,p = build_countries(botanical,butterflies)
    states = {ct:graph_vectors(b[ct],p[ct]) for ct in sorted(b)}
    source_cases = [r for v in b.values() for r in v]
    species = Counter(r["butterfly"] for r in source_cases)
    if (len(source_cases),len(species),len(b),sum(r["N"] for r in source_cases),
        sum(r["observed_count"] for r in source_cases))!=(26,16,11,453,231):
        raise ValueError("All26 source and butterfly GADM case margins drifted")
    total = math.prod(z["states"] for z in states.values())
    if total != 13826598912:
        raise ValueError("Degree-constrained graph count disagrees with original exact reference")
    outcomes={}
    for mode in METRICS:
        v={ct:per_country(state,mode,species) for ct,state in states.items()}
        denom={"record":453,"case":26,"species":16,"country":11}.get(mode)
        if denom is None:
            depth=int(mode[6:])
            denom=sum(r["N"]>=depth for r in source_cases)
        outcomes[mode]=exact_reference(v,denom)
    if (abs(outcomes["record"]["observed_minus_reference_pp"]-.2023547)>1e-6
        or abs(outcomes["case"]["observed_minus_reference_pp"]-6.8580651)>1e-6):
        raise ValueError("Prior archived original weighted and balanced expectation not replayed")
    bycountry={ct:per_country(state,"case",species) for ct,state in states.items()}
    omitted=[]
    for ct,group in bycountry.items():
        keep={name:row for name,row in bycountry.items() if name!=ct}
        denom=26-len(states[ct]["butterflies"])
        observed=sum((v["observed"] for v in keep.values()),Fraction(0))
        E=sum((v["expected"] for v in keep.values()),Fraction(0))
        row={"omitted_country":ct,"retained_cases":denom,
             "observed_minus_reference_pp":round(100*float((observed-E)/denom),7)}
        if ct in ("CU","DO"):
            row["exact_label_upper_tail_NOT_biological_p"]=exact_reference(
                keep,denom)["exact_label_upper_tail_NOT_biological_p"]
        omitted.append(row)
    return {
        "schema":"chocho_all26_country_degree_host_identity_weighting_posthoc_v1",
        "status":"POSTHOC_SOURCE_WEIGHTING_DEPENDENCE_NOT_BIOLOGICAL_CONFIRMATION",
        "n_country_butterfly_cases":26,
        "butterfly_species":16,"countries":11,
        "original_source_indexed_butterfly_records":453,
        "original_host_gadm_co_detected_records":231,
        "fixed_degree_graphs":str(total),
        "definitions":{
            "record":"GBIF indexed butterfly records weighted equally",
            "case":"26 species×country cells weighted equally",
            "species":"16 butterfly species weighted equally; each species countries averaged",
            "country":"11 countries weighted equally; each country's butterfly cases averaged",
            "case_n2":"Cases with ≥2 GBIF records only, but graph margins fixed for all 26",
            "case_n5":"Cases with ≥5 GBIF records only, but graph margins fixed for all 26",
            "case_n10":"Cases with ≥10 GBIF records only, but graph margins fixed for all 26"
        },
        "metrics":outcomes,
        "posthoc_leave_one_country_out_case_metric":omitted,
        "biology_warning":"Each exact tail concerns a constructed, often biologically impossible host-label permutation. It is NOT a biological p-value. GBIF co-indexing at administrative-unit scale is not demonstrated plant availability, larval feeding, population growth or recruitment.",
    }

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--botanical",required=True,type=Path)
    ap.add_argument("--butterfly48",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    raw_bot=a.botanical.read_bytes()
    raw_b=a.butterfly48.read_bytes()
    o=analyse(json.loads(raw_bot),json.loads(raw_b))
    o["source_sha256"]={"botanical":hashlib.sha256(raw_bot).hexdigest(),
                        "butterfly48":hashlib.sha256(raw_b).hexdigest()}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ALL26_HOST_WEIGHTING_AUDIT="+json.dumps({
       "metrics":o["metrics"],"posthoc_leave_one_country_out_case_metric":
        o["posthoc_leave_one_country_out_case_metric"]},separators=(",",":")),flush=True)

if __name__=="__main__":
    main()

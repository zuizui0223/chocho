#!/usr/bin/env python3
"""Audit retrospective host recording versus later butterfly record continuity.

This is a diagnostics-only analysis; host_group comes from 2001-2024 port
botanical catalogues and can include host records AFTER butterfly first records.
No estimates in this module identify historic host arrival or butterfly breeding.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent

SOURCES = {
    "butterfly_events": ("provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv",
                         "efbc770469a720e2dea82478db2f6196ba03b353"),
    "traits": ("data/frozen/figure_sources/anthropogenic_species_metrics.csv",
                "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"),
}

ALIEN = "ALIEN_CANDIDATE_RECORDED"
NATIVE = "NATIVE_ONLY_CLASSIFIED"

def read(name):
    path, expected = SOURCES[name]
    raw = (ROOT/path).read_bytes()
    digest = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if digest != expected:
        raise ValueError("Source Git blob changed: " + path + " (" + digest + ")")
    return list(csv.DictReader(raw.decode("utf-8").splitlines()))


def midranks(values):
    values = list(values)
    order = sorted(enumerate(values), key=lambda p: p[1])
    out = [0.0]*len(values)
    i = 0
    while i < len(order):
        j = i+1
        while j < len(order) and order[j][1] == order[i][1]:
            j += 1
        for k in range(i,j):
            out[order[k][0]] = (i+j+1)/2
        i=j
    return out


def correlation(x,y):
    if len(x) < 3:
        return None
    ax,bx = sum(x)/len(x), sum(y)/len(y)
    aa = sum((a-ax)**2 for a in x)
    bb = sum((b-bx)**2 for b in y)
    ab = sum((a-ax)*(b-bx) for a,b in zip(x,y))
    return ab/math.sqrt(aa*bb) if aa>0 and bb>0 else None


def boot(groups, draws=4000, seed=20261031):
    """Fixed-group bootstrap, not island/phylogeny/recording-process uncertainty."""
    s=seed
    vals=[]
    for _ in range(draws):
        total=0.0
        for _ in range(len(groups)):
            s=(s*1664525+1013904223) & 0xffffffff
            total += groups[s % len(groups)]["difference"]
        vals.append(total/len(groups))
    vals.sort()
    return [vals[int(.025*draws)],vals[int(.975*draws)]]


def analyze():
    events=read("butterfly_events")
    traits={x["species"]:x for x in read("traits")}
    assert len(events)==5665
    groups=defaultdict(list)
    views={name: {} for name in ("species","island","island_year")}
    expected={ALIEN:139,NATIVE:80,
              "NOT_EVALUABLE":5333,"UNRESOLVED_ONLY":28,
              "NO_MATCHING_HOST_RECORD":85}
    for event in events:
        label=event["host_group"]
        if label not in expected:
            raise ValueError("Unknown host group")
        groups[label].append(event)
        if label not in (ALIEN,NATIVE):
            continue
        out=event["any_followup_record"]
        if out not in ("True","False"):
            raise ValueError("Non-boolean outcome")
        y=int(out=="True")
        for name, key in (
            ("species", event["species"]),
            ("island", event["landmass_id"]),
            ("island_year", event["landmass_id"]+"|"+event["first_record_year"])
        ):
            dest=views[name].setdefault(key, {"a":[], "n":[]})
            dest["a" if label==ALIEN else "n"].append(y)
    assert {k:len(v) for k,v in groups.items()}==expected

    compared={}
    for name, records in views.items():
        pairs=[]
        for key,v in records.items():
            if not v["a"] or not v["n"]:
                continue
            rateA=sum(v["a"])/len(v["a"])
            rateN=sum(v["n"])/len(v["n"])
            pairs.append({
                "id":key,"n_alien":len(v["a"]),"n_native":len(v["n"]),
                "alien_reindex_rate":rateA,
                "native_reindex_rate":rateN,
                "difference":rateA-rateN
            })
        assert len(pairs)==(17 if name=="species" else 22)
        def mean_diff(rows):
            return sum(p["difference"] for p in rows)/len(rows) if rows else None
        compared[name]={
            "paired_groups":len(pairs),
            "mean_rate_difference":mean_diff(pairs),
            "groups_with_positive_difference":sum(x["difference"]>0 for x in pairs),
            "groups_with_negative_difference":sum(x["difference"]<0 for x in pairs),
            "groups_with_zero_difference":sum(x["difference"]==0 for x in pairs),
            "species_or_group_bootstrap_95":boot(pairs),
            "at_least_2_events_per_group":{
                "paired_groups":len([p for p in pairs if min(p["n_alien"],p["n_native"])>=2]),
                "mean_rate_difference":mean_diff([p for p in pairs if min(p["n_alien"],p["n_native"])>=2])
            }
        }
        if name=="species":
            compared[name]["wing_size_vs_within_species_difference_spearman"]=correlation(
                midranks([float(traits[p["id"]]["wing_size_proxy"]) for p in pairs]),
                midranks([p["difference"] for p in pairs])
            )
    result={
        "schema":"chocho_retrospective_host_5yr_reindex_matched_v1",
        "status":"DESCRIPTIVE_CONFOUNDING_DIAGNOSTIC_NOT_HISTORIC_HOST_OR_MOVEMENT_EVIDENCE",
        "source_git_blobs":{k:{"path":v[0],"sha":v[1]} for k,v in SOURCES.items()},
        "butterfly_first_events_total":len(events),
        "groups":{k:{
            "events":len(v),
            "reindexed_within_5y":sum(x["any_followup_record"]=="True" for x in v),
            "raw_reindex_rate":sum(x["any_followup_record"]=="True" for x in v)/len(v)
        } for k,v in groups.items()},
        "matched_comparisons":compared,
        "limitations":[
            "Host-group assignment pooled 2001-2024 flora; can postdate butterfly index year.",
            "Native-only classified does not establish matched biological resource quality.",
            "All 5665 first-indexed events are a selected presence cohort, not a full at-risk island panel.",
            "Within-species, within-island and within-island-year are overlapping, not independent replications.",
            "Group comparisons retain observation, island habitat, geography and sampling confounding.",
            "Differences change with >=2-per-arm restriction; do not interpret sign causally.",
            "Wing-size interaction is a comparison among just 17 same-species eligible contrasts.",
            "No demonstration of dispersal routes, reproductive establishment, or shipping introductions."
        ]
    }
    return result


def main():
    result=analyze()
    (HERE/"host_reindex_matched_results.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"groups":result["groups"],"matched":result["matched_comparisons"]},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

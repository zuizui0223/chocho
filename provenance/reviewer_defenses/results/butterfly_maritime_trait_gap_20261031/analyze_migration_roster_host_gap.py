#!/usr/bin/env python3
"""Analyze published-author migration roster against global island record gaps.

The paper's original 2021 migration Table S1 could not be fetched from Wiley.
An independently curated author-selected 2025 GBIF taxonKey proxy is already
source-verified in migration_gbif_predicate_audit.json. Those NOT listed in
the proxy are UNKNOWN, not non-migrant controls.

Crucially, both introduced-host regional fraction and port co-indexing are
contemporary/pooled exposures, potentially AFTER butterfly first records.
The endpoint measures distances between GBIF records, NOT actual flight.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
DIR = Path(__file__).resolve().parent
INPUTS = {
    "roster": (
        "provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
        "439dd1034ca16f96ff94222695d98dc39747c2e4"
    ),
    "gaps": (
        "provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/matched_species_screen.csv",
        "eb316611c58d5d08e37345f9f944a0dbf170a41e"
    ),
    "followup": (
        "provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv",
        "efbc770469a720e2dea82478db2f6196ba03b353"
    )
}
ROSTER_MIGRANT = "REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"
ROSTER_UNKNOWN = "NOT_LISTED_OR_TAXON_UNRESOLVED_UNKNOWN"
GROUPS = ("ROSTER_SELECTED_MIGRANT", "NOT_LISTED_UNKNOWN")
ALIEN = "ALIEN_CANDIDATE_RECORDED"
NATIVE = "NATIVE_ONLY_CLASSIFIED"
HISTORICAL_CAUSAL_CLAIM = False


def read_source(kind: str) -> list[dict[str,str]]:
    path, expected_blob = INPUTS[kind]
    raw = (ROOT / path).read_bytes()
    actual = hashlib.sha1(
        b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
    ).hexdigest()
    if actual != expected_blob:
        raise ValueError(f"Git blob changed for {kind}: {actual}")
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))


def label(value: str) -> str:
    if value == ROSTER_MIGRANT:
        return GROUPS[0]
    if value == ROSTER_UNKNOWN:
        return GROUPS[1]
    raise ValueError(f"Unrecognized migration source disposition: {value}")


def join() -> tuple[list[dict],dict[str,str]]:
    roster=read_source("roster")
    assert len(roster)==225
    mapping={}
    for row in roster:
        key=row["speciesKey"]
        if key in mapping:
            raise ValueError("Duplicate butterfly GBIF speciesKey")
        mapping[key]=label(row["author_2025_GBIF_requested_taxon_key_proxy"])
    table=[]
    for r in read_source("gaps"):
        key=r["speciesKey"]
        if key not in mapping:
            raise ValueError("193-species gap panel not represented in all225")
        line={
            "species":r["species"],
            "speciesKey":key,
            "source_migration":mapping[key],
            "n":int(r["n"]),
            "far500":int(r["far500"]),
            "far100":int(r["far100"]),
            "host":float(r["host"]),
            "ports":float(r["port"]),
            "wing":float(r["wing"]),
            "host_species":float(r["host_species"]),
            "native":float(r["native"])
        }
        assert 0 <= line["far500"] <= line["far100"] <= line["n"]
        assert line["n"]>0
        table.append(line)
    assert len(table)==193 and len({r["speciesKey"] for r in table})==193
    assert sum(r["n"] for r in table)==3136
    assert sum(r["far500"] for r in table)==164
    return table,mapping


def ranks(vals):
    order=sorted(enumerate(vals),key=lambda item:item[1])
    result=[0.0]*len(order)
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and order[j][1]==order[i][1]:
            j+=1
        mid=(i+j+1)/2
        for k in range(i,j):
            result[order[k][0]]=mid
        i=j
    return result


def spearman(a,b):
    if len(a)<3:
        return None
    x,y=ranks(a),ranks(b)
    mx=sum(x)/len(x);my=sum(y)/len(y)
    xy=sum((u-mx)*(v-my) for u,v in zip(x,y))
    xx=sum((u-mx)**2 for u in x)
    yy=sum((v-my)**2 for v in y)
    return xy/math.sqrt(xx*yy) if xx>0 and yy>0 else None


def lcg_draw(state,size):
    state=(state*1664525+1013904223)&0xffffffff
    return state,state%size


def quantile(v,p):
    if not v:
        return None
    a=sorted(v)
    pos=(len(a)-1)*p
    low=int(math.floor(pos))
    high=int(math.ceil(pos))
    return a[low]+(a[high]-a[low])*(pos-low)


def stats(rows):
    if not rows:
        raise ValueError("Empty source-migration group")
    events=sum(r["n"] for r in rows)
    gaps=sum(r["far500"] for r in rows)
    fraction=lambda k:[r[k]/r["n"] for r in rows]
    gap500=fraction("far500")
    gap100=fraction("far100")
    measures={
        "rho_host_500":spearman([r["host"] for r in rows],gap500),
        "rho_ports_500":spearman([r["ports"] for r in rows],gap500),
        "rho_wing_500":spearman([r["wing"] for r in rows],gap500),
        "rho_host_100":spearman([r["host"] for r in rows],gap100),
        "rho_ports_100":spearman([r["ports"] for r in rows],gap100),
        "rho_wing_100":spearman([r["wing"] for r in rows],gap100)
    }
    return {
        "species":len(rows),
        "island_first_index_events":events,
        "previous_record_gap_500km_events":gaps,
        "species_with_500km_gap":sum(r["far500"]>0 for r in rows),
        "gap_500_event_weighted":gaps/events,
        "gap_500_species_equal":sum(gap500)/len(gap500),
        "correlations":measures
    }


def bootstrap_rhos(rows,seed,draws=1200):
    names=list(stats(rows)["correlations"])
    samples={k:[] for k in names}
    state=seed
    for _ in range(draws):
        selected=[]
        for _ in range(len(rows)):
            state,i=lcg_draw(state,len(rows))
            selected.append(rows[i])
        est=stats(selected)["correlations"]
        for k,v in est.items():
            if v is not None and math.isfinite(v):
                samples[k].append(v)
    return {
        k:{"ci95":[quantile(samples[k],0.025),quantile(samples[k],0.975)],
           "valid_resamples":len(samples[k])}
        for k in names
    }


def retrospective_followup(mapping):
    events=read_source("followup")
    assert len(events)==5665
    rate=defaultdict(lambda:{"n":0,"reindexed":0})
    species_group=defaultdict(lambda:{"a":[],"n":[]})
    for r in events:
        group=mapping.get(r["speciesKey"])
        if group is None:
            raise ValueError("A butterfly first event cannot be linked to all225")
        status=r["host_group"]
        outcome=r["any_followup_record"]
        assert outcome in ("True","False")
        y=int(outcome=="True")
        rate[(group,status)]["n"]+=1
        rate[(group,status)]["reindexed"]+=y
        if status in (ALIEN,NATIVE):
            species_group[(group,r["species"])][
                "a" if status==ALIEN else "n"
            ].append(y)
    result={}
    for group in GROUPS:
        pairs=[]
        for (g,species),r in species_group.items():
            if g==group and r["a"] and r["n"]:
                pa=sum(r["a"])/len(r["a"])
                pn=sum(r["n"])/len(r["n"])
                pairs.append({"species":species,"n_alien":len(r["a"]),
                              "n_native":len(r["n"]),"difference":pa-pn})
        assert len(pairs)>0
        sample_mean=sum(p["difference"] for p in pairs)/len(pairs)
        ci=[]
        state=20261033 if group==GROUPS[0] else 20261034
        for _ in range(3000):
            score=0.0
            for _ in range(len(pairs)):
                state,i=lcg_draw(state,len(pairs))
                score+=pairs[i]["difference"]
            ci.append(score/len(pairs))
        qualifying=[p for p in pairs if p["n_alien"]>=2 and p["n_native"]>=2]
        result[group]={
            "host_groups":{
                name:{"n":rate[(group,name)]["n"],
                      "reindexed":rate[(group,name)]["reindexed"],
                      "rate":rate[(group,name)]["reindexed"]/rate[(group,name)]["n"]}
                for name in (ALIEN,NATIVE)
            },
            "within_species":{
                "paired_species":len(pairs),
                "mean_alien_minus_native_reindex_rate":sample_mean,
                "ci95_resampled":[quantile(ci,.025),quantile(ci,.975)],
                "positive":sum(p["difference"]>0 for p in pairs),
                "negative":sum(p["difference"]<0 for p in pairs),
                "zero":sum(p["difference"]==0 for p in pairs),
                "at_least_two_each":{"n":len(qualifying),
                    "mean_difference":sum(p["difference"] for p in qualifying)/len(qualifying) if qualifying else None}
            }
        }
    return result


def analyze():
    rows,roster=join()
    migrant=[r for r in rows if r["source_migration"]==GROUPS[0]]
    unknown=[r for r in rows if r["source_migration"]==GROUPS[1]]
    assert (len(migrant),len(unknown))==(106,87)
    assert sum(x["n"] for x in migrant)==2177
    assert sum(x["far500"] for x in migrant)==119
    groups={GROUPS[0]:migrant,GROUPS[1]:unknown}
    report={
        "schema":"chocho_migrant_selected_host_gap_followup_v1",
        "status":"SOURCE_VERIFIED_PROXY_WITH_NO_CAUSAL_HOST_OR_STEPPING_STONE_SUPPORT",
        "source_files":{k:{"path":v[0],"git_blob":v[1]} for k,v in INPUTS.items()},
        "cohort":"both_period_opportunity",
        "full_candidate_butterfly_species":len(roster),
        "gap_comparison_species":len(rows),
        "groups":{},
        "followup_retrospective_2001_2024_plant_pooled":retrospective_followup(roster),
        "not_listed_is_unknown":True,
        "not_actual_butterfly_displacement":True,
        "plant_records_can_postdate_butterfly":True,
        "not_actual_host_feeding_or_establishment":True,
        "notes":[
            "Migration proxy derives indirectly from a later author's GBIF species-key request.",
            "Original 2021 published review Table S1 is not redistributed or claimed independently recovered.",
            "Migration trait is species-level; populations can be resident in particular regions.",
            "Unlisted means undocumented/unresolved here, not biologically non-migratory.",
            "500 km gaps are between past and current GBIF recording sites; not flight distance.",
            "Contemporary host resource fraction is not historically dated island host availability.",
            "Pooled 2001–2024 port host records can follow butterfly first records.",
            "Species bootstrap CIs condition on these data and ignore shared islands/phylogeny.",
            "Event persistence analyses select previously indexed species-island pairs.",
            "Retrospective native-only controls are NOT equivalent biological resource controls.",
            "n thresholds and multiple correlations are exploratory, not independent confirmatory tests."
        ]
    }
    for i,group in enumerate(GROUPS):
        q=groups[group]
        report["groups"][group]={
            **stats(q),
            "bootstrap_95":bootstrap_rhos(q,20261031+i),
            "minimum_record_thresholds":{
                str(n):stats([r for r in q if r["n"]>=n])
                for n in (5,10,20)
            }
        }
    return report,rows


def main():
    report,rows=analyze()
    (DIR/"migrant_roster_host_gap_results.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    with (DIR/"migrant_roster_host_gap_193_species.csv").open("w",encoding="utf8",newline="") as f:
        fields=list(rows[0])
        writer=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    print(json.dumps({k:v for k,v in report.items()
                      if k in ("status","gap_comparison_species","groups")},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

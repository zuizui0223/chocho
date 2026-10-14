#!/usr/bin/env python3
"""Source-verifiable Kyoto 2012 pupation-timing comparison at FIXED native density.

Original Hashimoto & Ohgushi (2023) 30-cage/23-census larval CSV, Figshare
file 41146991. This is a descriptive, post hoc subset of a PUBLISHED density
experiment, not a new causal effect, time-varying resource experiment or adult
emergence finding. The source contains no greenhouse ID, source maternal
family, observed death causes, or longitudinal edible host supply.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

URL = "https://ndownloader.figshare.com/files/41146991"
MD5 = "0e9ff924e4ae38cb78529ec540887c42"
COLUMNS = (
    "plot.no","treat.no","s.density","a.density","date","total.s",
    "total.a","day","s1","s2","s3","s4","s5","spp","sp",
    "cumul.sp","a1","a2","a3","a4","a5","app","ap","cumul.ap",
)
CHECK_DAYS=(7,15,21,29,42)
DOI="10.1002/ece3.10164"


def _integer(cell, description: str) -> int:
    text = str(cell or "").strip()
    try:
        val = float(text)
    except ValueError as exc:
        raise ValueError(f"{description}: missing or invalid source integer") from exc
    if not val.is_integer() or val < 0:
        raise ValueError(f"{description}: invalid nonnegative source count")
    return int(val)


def _load(raw: bytes, *, checksum: str = MD5) -> dict[str, list[dict]]:
    if hashlib.md5(raw).hexdigest() != checksum:
        raise ValueError("Original Figshare MD5 mismatch")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if tuple(reader.fieldnames or []) != COLUMNS:
        raise ValueError("Original Figshare field schema changed")
    rows = list(reader)
    if len(rows) != 690:
        raise ValueError("Expected 690 original repeated cage-by-date rows")
    cages: dict[str, list[dict]] = defaultdict(list)
    combinations=Counter()
    for row in rows:
        cage=row["plot.no"]
        if not cage:
            raise ValueError("Missing original cage identity")
        cages[cage].append(row)
    if len(cages)!=30:
        raise ValueError("Expected precisely 30 independent whole cages")
    all_days=None
    for cid, subset in cages.items():
        ordered=sorted(subset,key=lambda r:_integer(r["day"], "day"))
        dates=[_integer(r["day"], "day") for r in ordered]
        if len(dates)!=23 or len(set(dates))!=23:
            raise ValueError("Missing/duplicate scheduled original census time")
        if all_days is None:
            all_days=dates
        elif dates!=all_days:
            raise ValueError("Different original census grids; date-matching required")
        first=ordered[0]
        s=_integer(first["s.density"], "initial Sericinus")
        a=_integer(first["a.density"], "initial Atrophaneura")
        combinations[(s,a)]+=1
        previous_a=previous_s=-1
        for row in ordered:
            if (_integer(row["s.density"],"S")!=s or
                    _integer(row["a.density"],"A")!=a):
                raise ValueError("Source treatment assignment changed within original cage")
            total_a=_integer(row["total.a"],"visible native A larvae")
            total_s=_integer(row["total.s"],"visible introduced S larvae")
            cumulative_a=_integer(row["cumul.ap"],"cum native pupae")
            cumulative_s=_integer(row["cumul.sp"],"cum introduced pupae")
            if total_a>a or total_s>s:
                raise ValueError("Observed larvae exceed original initial assignment")
            if not previous_a<=cumulative_a<=a:
                raise ValueError(
                    "Native A cumulative pupation dropped or exceeded starting cohort: "
                    f"source_cage={cid},day={row['day']},A_initial={a},S_initial={s},"
                    f"A_previous_pupae={previous_a},A_now_pupae={cumulative_a}"
                )
            # Original introduced S cumulative column is not a required
            # variable for the native A pupation-time estimand. Preserve any
            # anomalous decreases in a separate source-quality receipt; never
            # silently modify, smooth, carry forward, or drop original rows.
            if cumulative_s>s:
                raise ValueError(
                    f"S cumulative pupation exceeds initial assigned source cohort: "
                    f"cage={cid},day={row['day']},observed={cumulative_s},initial={s}"
                )
            previous_a,previous_s=cumulative_a,cumulative_s
        cages[cid]=ordered
    if len(combinations)!=15 or any(n!=2 for n in combinations.values()):
        raise ValueError("Original 15 randomized density combinations x 2 whole cages required")
    return cages


def analyze_source(raw: bytes, *, checksum: str = MD5) -> dict:
    cages=_load(raw,checksum=checksum)
    cage_data={}
    original_grid=sorted({
        _integer(row["day"],"day") for rows in cages.values() for row in rows
    })
    actual_cutoffs={}
    for planned_day in CHECK_DAYS:
        candidates=[day for day in original_grid if day<=planned_day]
        if not candidates:
            raise ValueError("Source starts later than specified observation window")
        actual_cutoffs[str(planned_day)]=max(candidates)
    actual_cutoffs["original_final"]=original_grid[-1]

    for cid, rows in cages.items():
        first=rows[0]
        a=_integer(first["a.density"],"initial A")
        s=_integer(first["s.density"],"initial S")
        visits={_integer(row["day"],"day"):row for row in rows}
        outcomes={}
        for planned, actual in actual_cutoffs.items():
            row=visits[actual]
            pup=_integer(row["cumul.ap"],"cumulative A pupae")
            visible=_integer(row["total.a"],"current A larvae")
            # The unaccounted quantity is descriptive only: visible larvae
            # may be incompletely detected; pupa fate and census counts are
            # not unique individuals followed to an adult.
            residual=a-visible-pup
            outcomes[planned]={
                "observed_day":actual,
                "recorded_cumulative_pupae":pup,
                "pupation_fraction_of_initial_A":pup/a if a>0 else None,
                "living_A_larvae_observed":visible,
                "initial_minus_visible_larvae_minus_cumulative_pupae":residual,
            }
        cage_data[cid]={
            "initial_A":a,"initial_S":s,
            "source_cage_id":cid,
            "outcomes":outcomes,
        }

    strata=defaultdict(lambda:defaultdict(list))
    for cid, row in cage_data.items():
        if row["initial_A"]>0:
            strata[row["initial_A"]][row["initial_S"]].append(row)
    # Choose only truly source-matched native densities; NO outcome selection.
    pairs={}
    for a, s_groups in sorted(strata.items()):
        if len(s_groups)!=2:
            continue
        levels=sorted(s_groups)
        if any(len(s_groups[s])!=2 for s in levels):
            raise ValueError("Fixed native density matched pair must have 2 source cages per S dose")
        low_s, high_s=levels
        comparisons={}
        for planned in actual_cutoffs:
            avgs={}
            for s in levels:
                rr=s_groups[s]
                avgs[str(s)]=sum(
                    r["outcomes"][planned]["pupation_fraction_of_initial_A"]
                    for r in rr
                )/len(rr)
            comparisons[planned]={
                "actual_day":actual_cutoffs[planned],
                "low_S_mean_pupation_fraction":avgs[str(low_s)],
                "high_S_mean_pupation_fraction":avgs[str(high_s)],
                "higher_minus_lower_S":avgs[str(high_s)]-avgs[str(low_s)],
            }
        pairs[str(a)]={
            "matched_original_A_neonates_per_cage":a,
            "lower_S_initial_n":low_s,
            "higher_S_initial_n":high_s,
            "cages_per_density_cell":2,
            "source_cage_ids_lower_S":sorted(r["source_cage_id"] for r in s_groups[low_s]),
            "source_cage_ids_higher_S":sorted(r["source_cage_id"] for r in s_groups[high_s]),
            "by_day":comparisons,
        }
    if not pairs:
        raise ValueError("Original source contains no exact-native-density comparison")
    summaries={}
    for planned in actual_cutoffs:
        diffs=[pair["by_day"][planned]["higher_minus_lower_S"] for pair in pairs.values()]
        summaries[planned]={
            "actual_census_day":actual_cutoffs[planned],
            "number_of_source_matched_A_density_strata":len(diffs),
            "mean_higher_minus_lower_S_pupation_fraction_across_A_density_strata":sum(diffs)/len(diffs),
            "n_strata_high_S_lower_pupation":sum(d<0 for d in diffs),
            "n_strata_high_S_higher_pupation":sum(d>0 for d in diffs),
            "n_strata_tied":sum(d==0 for d in diffs),
            "all_source_matched_differences":diffs,
        }

    negative_residuals=[
        (cid, planned, values["initial_minus_visible_larvae_minus_cumulative_pupae"])
        for cid, row in cage_data.items()
        for planned, values in row["outcomes"].items()
        if values["initial_minus_visible_larvae_minus_cumulative_pupae"]<0
    ]
    s_cumulative_declines=[]
    for cid, rows in sorted(cages.items()):
        previous=-1
        for r in rows:
            value=_integer(r["cumul.sp"],"source cumulative S pupae")
            if value<previous:
                s_cumulative_declines.append({
                    "source_cage_id":cid,
                    "current_day":_integer(r["day"],"day"),
                    "previous_cumulative_S":previous,
                    "current_cumulative_S":value,
                    "original_initial_S":_integer(rows[0]["s.density"],"initial S"),
                })
            previous=value
    original_terminal_rates={
        "A":sum(row["outcomes"]["original_final"]["recorded_cumulative_pupae"]
                for row in cage_data.values())/
            sum(row["initial_A"] for row in cage_data.values()),
        "S":(
            None if s_cumulative_declines else
            sum(_integer(cages[cid][-1]["cumul.sp"],"cum S pupae")
                for cid in cages)/
            sum(_integer(cages[cid][0]["s.density"],"initial S") for cid in cages)
        ),
    }
    return {
        "schema":"chocho_kyoto_2012_matched_A_density_pupation_timeline_v01",
        "status":"AUTHOR_ORIGINAL_SOURCE_DESCRIPTIVE_POSTHOC_NOT_NEW_CAUSAL_MECHANISM",
        "publication_DOI":DOI,
        "original_source_url":URL,
        "original_source_md5":checksum,
        "original_year":2012,
        "whole_cages":len(cages),
        "cage_x_visit_original_rows":sum(len(x) for x in cages.values()),
        "original_census_day_grid":original_grid,
        "planned_days_and_actual_source_days":actual_cutoffs,
        "native_density_matched_strata":pairs,
        "equal_native_stratum_weight_summary_by_day":summaries,
        "source_total_pupation_probability_by_species":original_terminal_rates,
        "introduced_S_cumulative_pupation_decreases_from_original":s_cumulative_declines,
        "introduced_S_terminal_pupation_summary_withheld_if_inconsistent":bool(s_cumulative_declines),
        "all_observed_cage_cutoff_negative_unaccounted_residuals":negative_residuals,
        "source_cage_records":cage_data,
        "warnings":[
            "Pupation success and its calendar progression, NOT adult eclosion, adult recruitment or confirmed individual death times",
            "Introduced S cumulative pupation in source may decrease: retain every anomaly as source quality, withhold species-level S terminal summary",
            "The same 2023 source paper has already reported density effects on final pupation survival and developmental time",
            "Post hoc S high vs S low matched A comparisons use the ORIGINAL density design: total larval headcount also differs",
            "Only 2 whole cages per original density cell; greenhouse IDs and maternal lineage IDs absent from original Figshare CSV",
            "Stage counts are cage-level censuses, not individually tagged event times",
            "No repeated host edible young/mature leaf supply, host age, quality, or timed resource clamp in original file",
            "Differences at early day windows can reflect pupation timing, eventual survival or censoring; do not label mortality onset",
            "An original larva not currently visible or cumulatively pupated is not necessarily confirmed dead",
            "No day-specific p-values, interaction significance or causal food mediator inferred",
            "No changes to GEB main/pr38 or 30-cage future intervention study",
        ],
        "new_butterfly_intervention_or_fitness_observations":False,
        "GEB_PR38_untouched":True,
    }


def download_original() -> bytes:
    with urlopen(Request(URL, headers={"User-Agent":"chocho-exact-original-kyoto-pupation-timeline-audit/1.0"}),timeout=45) as response:
        raw=response.read(200000)
    if len(raw)>=200000:
        raise ValueError("Original larval data unexpectedly large")
    return raw


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receipt",required=True,type=Path)
    args=ap.parse_args()
    result=analyze_source(download_original())
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        k:v for k,v in result.items() if k!="source_cage_records"
    }, indent=2))


if __name__=="__main__":
    main()

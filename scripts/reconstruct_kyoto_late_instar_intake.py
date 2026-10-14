#!/usr/bin/env python3
"""Replicate SOURCE-DESCRIBED 2-day dry-leaf intake chronology.

Hashimoto & Ohgushi (2023), Section 2.2, and its Figshare original:
- leaves are replaced every two days, fresh 'given' then 'taken'
- leaf dry mass: ln(dry mg)=-1.48+0.90 ln(fresh mg)
- single larvae follow the group-reared first/second instars.
The analysis is descriptively conditional on post-third-instar tracked larvae.
It does NOT observe the counterfactual host supply or competitive survival.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import audit_kyoto_calendar_demand_identifiability as source

DRY_CONVERSION_INTERCEPT = -1.48
DRY_CONVERSION_SLOPE = 0.90
PUBLISHED_LATE_DRY_MG_PER_LARVA = {"A": 1003.23, "S": 394.25}
VALID_STAGE = {"1", "2", "3", "4", "5", "pp", "p"}


def date(row):
    return datetime.strptime(row["date"], "%Y/%m/%d").date()


def fresh_g_to_dry_mg(fresh_g: float) -> float:
    if not math.isfinite(fresh_g) or fresh_g < 0:
        raise ValueError("fresh leaf wet mass must be finite, nonnegative grams")
    if fresh_g == 0:
        return 0.0
    # Both fresh and dry masses in the published log-log calibration are
    # expressed in grams; return dry milligrams only AFTER conversion.
    # Converting fresh grams to mg before applying a 0.90 exponent creates
    # an erroneous factor of 1000**0.10 (~1.995) in the inferred mass.
    return 1000.0 * math.exp(
        DRY_CONVERSION_INTERCEPT + DRY_CONVERSION_SLOPE * math.log(fresh_g)
    )


def reconstruct(rows: list[dict[str, str]], *, species: str) -> dict:
    if species not in PUBLISHED_LATE_DRY_MG_PER_LARVA:
        raise ValueError("species must be A or S")
    by_id = defaultdict(list)
    for row in rows:
        if str(row.get("stage") or "").strip() not in VALID_STAGE:
            raise ValueError("unrecognized original stage code")
        by_id[row["id"]].append(row)
    for ident in by_id:
        by_id[ident].sort(key=date)
        raw_by_day = defaultdict(list)
        for record in by_id[ident]:
            raw_by_day[date(record)].append(record)
        collisions = {
            day.isoformat(): [
                {k: rec.get(k, "") for k in ("date", "stage", "taken", "given", "memo", "no.individuals")}
                for rec in recs
            ]
            for day, recs in raw_by_day.items() if len(recs) > 1
        }
        if collisions and not all(source.number(rec["no.individuals"]) == 10
                                  for rec in by_id[ident]):
            raise ValueError(
                "duplicate identity/calendar-day records require original-code reconciliation: "
                + json.dumps({"id": ident, "dates": collisions})
            )
        # The mass analysis excludes the original founder group; duplicate
        # stage-transition snapshots within that group are preserved, never
        # interpreted as two independent feeding intervals.
    group_ids = [
        ident for ident, group in by_id.items()
        if group[0]["stage"] == "1"
        and source.number(group[0]["no.individuals"]) == 10
    ]
    if len(group_ids) != 1:
        raise ValueError("not exactly one verified 10-neonate founder group")
    founder = group_ids[0]
    founder_duplicate_dates = {
        day.isoformat(): [
            {"stage": r["stage"], "taken": r.get("taken", ""), "given": r.get("given", "")}
            for r in by_id[founder] if date(r) == day
        ]
        for day in sorted({date(r) for r in by_id[founder]})
        if sum(date(r) == day for r in by_id[founder]) > 1
    }
    hatch_date = date(by_id[founder][0])
    individual_ids = sorted(ident for ident, rr in by_id.items()
                            if ident != founder and ident.startswith(founder)
                            and source.number(rr[0]["no.individuals"]) == 1)
    if not individual_ids:
        raise ValueError("single-larva follow-up cohort absent")
    skipped = Counter()
    events = []
    by_individual = defaultdict(list)
    for ident in individual_ids:
        sequence = by_id[ident]
        for previous, current in zip(sequence, sequence[1:]):
            if "redundant" in (previous.get("memo") or "").lower() or (
                "redundant" in (current.get("memo") or "").lower()
            ):
                skipped["memo_redundant_adjacent"] += 1
                continue
            if source.number(previous["no.individuals"]) != 1 or source.number(current["no.individuals"]) != 1:
                skipped["not_single_individual"] += 1
                continue
            gap = (date(current) - date(previous)).days
            if gap != 2:
                skipped["not_exactly_2day_pair"] += 1
                continue
            if previous["stage"] not in {"3", "4", "5"}:
                skipped["previous_stage_not_feeding_late_instar"] += 1
                continue
            wet_given = source.number(previous["given"])
            wet_remaining = source.number(current["taken"])
            if wet_given is None or wet_remaining is None:
                skipped["no_matched_given_taken"] += 1
                continue
            dry_cons = fresh_g_to_dry_mg(wet_given) - fresh_g_to_dry_mg(wet_remaining)
            if dry_cons < -1e-8:
                skipped["negative_dry_consumption"] += 1
                continue
            if dry_cons < 0:
                dry_cons = 0.0
            event = {
                "species": species,
                "original_larva_id": ident,
                "feeding_stage_start": previous["stage"],
                "feeding_stage_end": current["stage"],
                "end_day_since_founder_hatch": (date(current) - hatch_date).days,
                "observed_interval_days": gap,
                "consumed_leaf_dry_mg": dry_cons,
            }
            events.append(event)
            by_individual[ident].append(event)
    if not events:
        raise ValueError("no valid original matched food-intake intervals")
    # The earliest cohort date is not the first date of individually tracked third instars.
    n_completed = sum(
        any(row["stage"] == "p" for row in by_id[i]) for i in individual_ids
    )
    summaries = []
    for ident in individual_ids:
        intervals = by_individual.get(ident, [])
        total = sum(x["consumed_leaf_dry_mg"] for x in intervals)
        if not intervals or total <= 0:
            continue
        accumulated = 0.0
        half_age = None
        for event in sorted(intervals, key=lambda e: e["end_day_since_founder_hatch"]):
            accumulated += event["consumed_leaf_dry_mg"]
            if accumulated >= total / 2 and half_age is None:
                half_age = event["end_day_since_founder_hatch"]
        summaries.append({
            "original_larva_id": ident,
            "tracked_2day_intake_intervals": len(intervals),
            "observed_late_dry_intake_total_mg": total,
            "day_since_hatch_reaching_half_own_tracked_late_intake": half_age,
            "pupa_code_recorded": any(r["stage"] == "p" for r in by_id[ident]),
        })
    final_age = defaultdict(lambda: {"sum_dry_mg": 0.0, "intervals": 0})
    for event in events:
        entry = final_age[event["end_day_since_founder_hatch"]]
        entry["sum_dry_mg"] += event["consumed_leaf_dry_mg"]
        entry["intervals"] += 1
    per_day = [{
        "day_since_founder_group_start": age,
        "original_individual_intervals": info["intervals"],
        "sum_dry_mg_across_tracked_individuals": round(info["sum_dry_mg"], 5),
        "mean_dry_mg_per_observed_interval": round(info["sum_dry_mg"] / info["intervals"], 5),
    } for age, info in sorted(final_age.items())]
    individual_means = [row["observed_late_dry_intake_total_mg"] for row in summaries]
    observed_mean = statistics.mean(individual_means) if individual_means else None
    original_published = PUBLISHED_LATE_DRY_MG_PER_LARVA[species]
    return {
        "status": "SOURCE_DESCRIBED_LAB_CONSUMPTION_RECONSTRUCTION_NOT_COMPETITION",
        "species": species,
        "source_trial_founder_id": founder,
        "founder_group_same_day_stage_transitions": founder_duplicate_dates,
        "founder_neonate_group_start": hatch_date.isoformat(),
        "n_single_larva_ids": len(individual_ids),
        "n_pupa_coded_ids": n_completed,
        "n_individuals_with_positive_usable_intake": len(summaries),
        "n_intervals_included": len(events),
        "interval_events": events,  # source-matched events for separate age-alignment audit
        "skipped_adjacent_pairs": dict(skipped),
        "published_table1_late_dry_mg_per_larva": original_published,
        "mean_usable_observed_late_dry_mg": observed_mean,
        "mean_relative_to_published_reference": observed_mean / original_published if observed_mean else None,
        "table1_dry_intake_within_2pct": (
            abs(observed_mean / original_published - 1) <= 0.02
            if observed_mean else False
        ),
        "median_day_since_hatch_half_observed_late_intake": (
            statistics.median(x["day_since_hatch_reaching_half_own_tracked_late_intake"]
                              for x in summaries) if summaries else None
        ),
        "calendar_series": per_day,
        "per_individual": summaries,
        "interpretation": (
            "Within-species, independently reared 2016 larval wet-leaf given/taken "
            "measurement pairs, converted using the published dry-leaf calibration. "
            "Founder group-to-individual splits require special treatment; excluded "
            "intervals and pupation markers are explicitly counted. Intervals missing "
            "or marked redundant are not imputed. Values cannot be assumed identical "
            "to the published Table 1 after filtering without a complete source-code "
            "reconciliation. No larva was experimentally exposed to different "
            "food-pulse dates in these data."
        ),
    }


def audit():
    meta = json.loads(source.fetch_bytes(source.ARTICLE_URL))
    if meta.get("id") != source.ARTICLE_ID:
        raise ValueError("Figshare source article ID mismatch")
    names = {record.get("name"): record for record in meta.get("files", [])}
    provenance, results = {}, {}
    for species, name in source.FOOD_FILES.items():
        if name not in names:
            raise ValueError("original food-demand data missing")
        record = names[name]
        raw = source.fetch_bytes(record["download_url"])
        check_md5 = __import__("hashlib").md5(raw).hexdigest()
        published_md5 = record.get("computed_md5") or record.get("supplied_md5")
        if published_md5 and check_md5 != published_md5:
            raise ValueError("upstream original checksum mismatch")
        columns, rows = source.read_csv(raw)
        if "stage" not in columns:
            raise ValueError("source contains no stage column")
        provenance[species] = {"original_file": name, "original_file_id": record.get("id"), "md5": check_md5}
        results[species] = reconstruct(rows, species=species)
    return {
        "schema": "chocho_kyoto_late_instar_intake_reconstruction_v01",
        "source": "Hashimoto & Ohgushi (2023) Figshare 23170898 and paper section 2.2",
        "mass_conversion": "dry mg = 1000 * exp(-1.48 + 0.90 ln(fresh g))",
        "source_identity": provenance,
        "per_species": results,
        "calendar_start_date_difference_days_A_minus_S": (
            datetime.fromisoformat(results["A"]["founder_neonate_group_start"]).date() -
            datetime.fromisoformat(results["S"]["founder_neonate_group_start"]).date()
        ).days,
        "measured_on_common_contemporaneous_plants": False,
        "original_2012_mixed_cage_supply_measured": False,
        "new_randomized_resource_timing_effect_estimated": False,
        "main_GEB_unchanged": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit()
    except Exception as exc:
        result = {"schema": "chocho_kyoto_late_instar_intake_reconstruction_v01",
                  "status": "SOURCE_BLOCKED_OR_INVALID",
                  "reason": f"{type(exc).__name__}: {str(exc)[:250]}",
                  "new_randomized_resource_timing_effect_estimated": False}
    else:
        result["status"] = "SOURCE_RECONSTRUCTED_DESCRIPTIVE_ONLY"
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

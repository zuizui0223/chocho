#!/usr/bin/env python3
"""Original-source audit of Kyoto butterfly calendar-stage and food-demand records.

The data were published by Hashimoto & Ohgushi (2023), Figshare 23170898.
All outputs are descriptive source support; this does NOT estimate adult
emergence, assigned food-pulse timing, or a causal interaction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import statistics
from datetime import datetime
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

ARTICLE_ID = 23170898
ARTICLE_URL = f"https://api.figshare.com/v2/articles/{ARTICLE_ID}"
LARVAL_FILE = "Hashimoto_and_Ohgushi_2023_larvaldata.csv"
FOOD_FILES = {
    "A": "Hashimoto_and_Ohgushi_2023_AFoodDemand.csv",
    "S": "Hashimoto_and_Ohgushi_2023_SFoodDemand.csv",
}
LARVAL_MD5 = "0e9ff924e4ae38cb78529ec540887c42"
EXPECTED_LARVAL_COLUMNS = (
    "plot.no", "treat.no", "s.density", "a.density", "date", "total.s",
    "total.a", "day", "s1", "s2", "s3", "s4", "s5", "spp",
    "sp", "cumul.sp", "a1", "a2", "a3", "a4", "a5",
    "app", "ap", "cumul.ap",
)
EXPECTED_FOOD_COLUMNS = ("id", "date", "taken", "given", "mass", "no.individuals", "memo")


def fetch_bytes(url: str, limit: int = 2_000_000) -> bytes:
    with urlopen(Request(url, headers={"User-Agent": "chocho-kyoto-source-audit/1.0"}),
                 timeout=40) as response:
        raw = response.read(limit + 1)
    if not raw or len(raw) > limit:
        raise ValueError("empty/oversize upstream original")
    return raw


def read_csv(raw: bytes) -> tuple[list[str], list[dict[str, str]]]:
    decoded = raw.decode("utf-8-sig")
    head = decoded[:2000]
    try:
        dialect = csv.Sniffer().sniff(head, delimiters=",;\t")
        delimiter = dialect.delimiter
    except csv.Error:
        delimiter = ","
    stream = csv.DictReader(io.StringIO(decoded), delimiter=delimiter)
    cols = list(stream.fieldnames or [])
    if not cols or len(set(cols)) != len(cols):
        raise ValueError("no header or duplicated columns")
    rows = list(stream)
    if not rows or any(None in row for row in rows):
        raise ValueError("empty or ragged source CSV")
    return cols, rows


def number(cell: str) -> float | None:
    t = str(cell or "").strip()
    if not t:
        return None
    try:
        value = float(t)
    except ValueError:
        return None
    if not __import__("math").isfinite(value):
        return None
    return value


def _late(row: dict[str, str], prefix: str) -> int:
    counts = [number(row[f"{prefix}{n}"]) for n in (3, 4, 5)]
    if any(x is None or x < 0 or not x.is_integer() for x in counts):
        raise ValueError("missing/invalid 3rd–5th instar count")
    return int(sum(counts))


def _positive_initial(row: dict[str, str], letter: str) -> bool:
    value = number(row[f"{letter.lower()}.density"])
    if value is None or value < 0 or not value.is_integer():
        raise ValueError("missing/invalid assigned initial density")
    return value > 0


def summarize_stages(columns: list[str], rows: list[dict[str, str]]) -> dict:
    if tuple(columns) != EXPECTED_LARVAL_COLUMNS or len(rows) != 690:
        raise ValueError("published larval source structure is not the frozen 690-row schema")
    by_cage: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_cage[str(row["plot.no"])].append(row)
    if len(by_cage) != 30:
        raise ValueError("expected 30 original cages, not 30 larvae/plant replicates")
    paired_cages = []
    cage_counts = Counter()
    paired_day_status = Counter()
    first_day_differences: list[float] = []
    for cage, group in sorted(by_cage.items()):
        ordered = sorted(group, key=lambda x: float(x["day"]))
        days = [number(row["day"]) for row in ordered]
        if len(ordered) != 23 or any(d is None for d in days) or len(set(days)) != 23:
            raise ValueError("cage has missing or duplicate census days")
        a_present = _positive_initial(ordered[0], "a")
        s_present = _positive_initial(ordered[0], "s")
        cage_counts[f"initial_A_{int(a_present)}_S_{int(s_present)}"] += 1
        if not (a_present and s_present):
            continue
        observed_a = [float(row["day"]) for row in ordered if _late(row, "a") > 0]
        observed_s = [float(row["day"]) for row in ordered if _late(row, "s") > 0]
        paired = {
            "cage": cage,
            "first_A_late_observed_day": min(observed_a) if observed_a else None,
            "first_S_late_observed_day": min(observed_s) if observed_s else None,
            "both_late_on_same_census_n": sum(
                _late(row, "a") > 0 and _late(row, "s") > 0 for row in ordered
            ),
            "A_late_census_n": len(observed_a),
            "S_late_census_n": len(observed_s),
        }
        for row in ordered:
            a, s = _late(row, "a") > 0, _late(row, "s") > 0
            if a and s:
                paired_day_status["both_late"] += 1
            elif a:
                paired_day_status["A_only_late"] += 1
            elif s:
                paired_day_status["S_only_late"] += 1
            else:
                paired_day_status["neither_late"] += 1
        if observed_a and observed_s:
            difference = min(observed_s) - min(observed_a)
            paired["first_observed_S_minus_A_late_day"] = difference
            first_day_differences.append(difference)
        else:
            paired["first_observed_S_minus_A_late_day"] = None
        paired_cages.append(paired)
    return {
        "n_source_cages": 30,
        "n_cage_by_initial_species_presence": dict(sorted(cage_counts.items())),
        "n_initially_mixed_cages": len(paired_cages),
        "n_mixed_cages_both_species_observed_late": len(first_day_differences),
        "observed_pairwise_S_minus_A_first_late_day_median": (
            statistics.median(first_day_differences) if first_day_differences else None
        ),
        "paired_late_onset_sign_counts": {
            "S_first_late_earlier": sum(x < 0 for x in first_day_differences),
            "same_first_observed_day": sum(x == 0 for x in first_day_differences),
            "A_first_late_earlier": sum(x > 0 for x in first_day_differences),
        },
        "mixed_cage_census_overlap_counts": dict(sorted(paired_day_status.items())),
        "per_mixed_cage": paired_cages,
        "boundary": (
            "Observed 3rd-5th instar presence is a repeated cage-level proxy for "
            "high-demand developmental stages, not food intake, exact molt-onset "
            "time, emergence success or a randomized supply window. Paired first "
            "onsets condition on both species surviving to late instars."
        ),
    }


def summarize_food(columns: list[str], rows: list[dict[str, str]]) -> dict:
    missing = set(EXPECTED_FOOD_COLUMNS) - set(columns)
    if missing:
        raise ValueError(f"required published feeding-sheet columns absent: {sorted(missing)}")
    ids: dict[str, set[str]] = defaultdict(set)
    dates: Counter[str] = Counter()
    no_individuals: Counter[str] = Counter()
    memo: Counter[str] = Counter()
    stages: Counter[str] = Counter()
    stage_by_group_size: Counter[str] = Counter()
    first_dates: dict[str, tuple[str, str]] = {}
    records_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    incomplete = Counter()
    for row in rows:
        ident = str(row["id"] or "").strip()
        date = str(row["date"] or "").strip()
        if not ident or not date:
            raise ValueError("food-demand observation without larva/group ID or date")
        ids[ident].add(date)
        dates[date] += 1
        records_by_id[ident].append(row)
        stage = str(row.get("stage") or "").strip()
        stages[stage or "blank"] += 1
        if "stage" in columns and not stage:
            incomplete["missing_stage"] += 1
        number_group = number(row["no.individuals"])
        if number_group is None or number_group <= 0 or not number_group.is_integer():
            incomplete["invalid_group_size"] += 1
        else:
            no_individuals[str(int(number_group))] += 1
            stage_by_group_size[f"{stage or 'blank'}|n={int(number_group)}"] += 1
        for field in ("taken", "given", "mass"):
            if number(row[field]) is None:
                incomplete[f"missing_or_non_numeric_{field}"] += 1
        memo[(row.get("memo") or "").strip().lower() or "blank"] += 1
    for ident, seq in sorted(records_by_id.items()):
        try:
            ordered = sorted(seq, key=lambda row: datetime.strptime(row["date"], "%Y/%m/%d"))
        except ValueError:
            ordered = sorted(seq, key=lambda row: row["date"])
        first_dates[ident] = (
            str(ordered[0]["date"]), str(ordered[0].get("stage", ""))
        )
    duplicates = {
        ident: [date for date, n in Counter(row["date"] for row in seq).items() if n > 1]
        for ident, seq in records_by_id.items()
        if any(n > 1 for n in Counter(row["date"] for row in seq).values())
    }
    formatted_dates = []
    for date in dates:
        try:
            formatted_dates.append(datetime.strptime(date, "%Y/%m/%d").date())
        except ValueError:
            try:
                formatted_dates.append(datetime.fromisoformat(date).date())
            except ValueError:
                incomplete["invalid_date"] += 1
    return {
        "raw_rows": len(rows),
        "raw_columns": columns,
        "unique_food_trial_ids": len(ids),
        "distinct_observation_dates": len(dates),
        "dates_with_most_rows": dates.most_common(3),
        "date_range_iso": (
            [min(formatted_dates).isoformat(), max(formatted_dates).isoformat()]
            if len(formatted_dates) == len(dates) else None
        ),
        "distinct_trial_start_dates": dict(Counter(pair[0] for pair in first_dates.values())),
        "stage_code_counts": dict(stages),
        "stage_by_group_size_counts": dict(stage_by_group_size),
        "per_trial_first_date_and_stage": dict(first_dates),
        "duplicate_ID_date_rows": duplicates,
        "group_size_counts": dict(sorted(no_individuals.items())),
        "memo_counts": dict(memo.most_common(12)),
        "data_quality": dict(sorted(incomplete.items())),
        "trial_ids_with_repeated_dates": sum(len(value) > 1 for value in ids.values()),
        "illustrative_original_rows_first_two": [
            {key: row.get(key, "") for key in EXPECTED_FOOD_COLUMNS}
            for row in rows[:2]
        ],
        "boundary": (
            "Source-paper Section 2.2 specifies 2-day leaves given then removed, "
            "converted individually from fresh to dry mass before subtracting. "
            "A valid consumption curve requires previous-given/current-taken "
            "pairing, mixed group-to-individual handling, source unit checks, "
            "unambiguous stage assignment and explicit censoring of redundant "
            "rows. These isolated 2016 food assays are not the 2012 mixed-cage "
            "food availability series and do not identify timing causality."
        ),
    }


def audit() -> dict:
    meta = json.loads(fetch_bytes(ARTICLE_URL))
    if meta.get("id") != ARTICLE_ID:
        raise ValueError("source article identity mismatch")
    all_files = {x.get("name"): x for x in meta.get("files", [])}
    required = [LARVAL_FILE, *FOOD_FILES.values()]
    if any(name not in all_files for name in required):
        raise ValueError("required author-original food or larval CSV missing")
    source = {}
    tables = {}
    for name in required:
        record = all_files[name]
        raw = fetch_bytes(record["download_url"], limit=5_000_000)
        actual = hashlib.md5(raw).hexdigest()
        pinned = record.get("computed_md5") or record.get("supplied_md5")
        if pinned and pinned != actual:
            raise ValueError(f"original Figshare checksum changed: {name}")
        if name == LARVAL_FILE and actual != LARVAL_MD5:
            raise ValueError("existing independent larval MD5 receipt differs")
        source[name] = {
            "figshare_file_id": record.get("id"),
            "source_md5": actual,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "size_bytes": len(raw),
        }
        tables[name] = read_csv(raw)
    stages = summarize_stages(*tables[LARVAL_FILE])
    food = {species: summarize_food(*tables[name]) for species, name in FOOD_FILES.items()}
    return {
        "schema": "chocho_kyoto_calendar_demand_source_audit_v01",
        "source_article_id": ARTICLE_ID,
        "source_doi": "10.6084/m9.figshare.23170898",
        "source_files": source,
        "original_cage_stage_chronology": stages,
        "original_single_species_food_trial_chronology": food,
        "not_independent_empirical_competition_effect": True,
        "source_aligned_feeding_demand_and_cage_calendar_validated": False,
        "randomized_leaf_supply_timing_present": False,
        "flight_capable_adult_emergence_present": False,
        "main_GEB_submission_unchanged": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = audit()
    except Exception as exc:
        result = {
            "schema": "chocho_kyoto_calendar_demand_source_audit_v01",
            "status": "SOURCE_BLOCKED_OR_INVALID",
            "reason": f"{type(exc).__name__}: {str(exc)[:240]}",
            "new_competition_effect_estimated": False,
        }
    else:
        result["status"] = "SOURCE_VERIFIED_DESCRIPTIVE_STAGE_SUPPORT_ONLY"
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

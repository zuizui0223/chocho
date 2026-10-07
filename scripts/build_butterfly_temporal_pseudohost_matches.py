#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


def parse_year(value: str):
    m = re.search(r"(?<!\d)(1[5-9]\d{2}|20[0-2]\d)(?!\d)", str(value or ""))
    return None if not m else int(m.group(1))


def load_sinas(path: Path, keep_names: set[str]):
    stats = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for row in reader:
            name = str(row.get("Taxon") or "").strip()
            if name not in keep_names:
                continue
            s = stats.setdefault(name, {"records": 0, "dated_records": 0, "earliest_year": None})
            s["records"] += 1
            y = parse_year(row.get("eventDate"))
            if y is not None:
                s["dated_records"] += 1
                if s["earliest_year"] is None or y < s["earliest_year"]:
                    s["earliest_year"] = y
    return stats


def sd(vals):
    vals = [float(x) for x in vals if x is not None and math.isfinite(float(x))]
    if len(vals) < 2:
        return 1.0
    m = sum(vals) / len(vals)
    v = sum((x - m) ** 2 for x in vals) / (len(vals) - 1)
    return max(math.sqrt(v), 1e-9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event-host-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--introduced-pool-csv", type=Path, required=True)
    ap.add_argument("--sinas-csv", type=Path, required=True)
    ap.add_argument("--matches-per-host", type=int, default=10)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    events = list(csv.DictReader(args.event_host_csv.open(newline="", encoding="utf-8")))
    if not events:
        raise RuntimeError("empty event-host input")

    focal_hosts = defaultdict(set)
    with args.insect_host_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            sp = str(row["insect_species"]).strip()
            pid = str(row["accepted_plant_name_id"]).strip()
            if sp and pid:
                focal_hosts[sp].add(pid)

    pool_meta = {}
    pool_by_region = defaultdict(list)
    with args.introduced_pool_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pid = str(row["plant_name_id"]).strip()
            name = str(row["accepted_name"]).strip()
            region = str(row["area_code_l3"]).strip()
            if not pid or not name or not region:
                continue
            pool_meta[pid] = {
                "plant_name_id": pid,
                "accepted_name": name,
                "family": str(row.get("family") or "").strip(),
                "introduced_wgsrpd3_count": int(float(row["introduced_wgsrpd3_count"])),
            }
            pool_by_region[region].append(pid)

    pool_names = {m["accepted_name"] for m in pool_meta.values()}
    sinas = load_sinas(args.sinas_csv, pool_names)

    eligible_ids = [
        pid for pid, meta in pool_meta.items()
        if meta["accepted_name"] in sinas
        and sinas[meta["accepted_name"]]["records"] > 0
        and sinas[meta["accepted_name"]]["earliest_year"] is not None
    ]
    if not eligible_ids:
        raise RuntimeError("no WCVP-SInAS matched introduced plants")

    breadth_vals = [math.log1p(pool_meta[pid]["introduced_wgsrpd3_count"]) for pid in eligible_ids]
    record_vals = [math.log1p(sinas[pool_meta[pid]["accepted_name"]]["records"]) for pid in eligible_ids]
    dated_vals = [math.log1p(sinas[pool_meta[pid]["accepted_name"]]["dated_records"]) for pid in eligible_ids]
    year_vals = [sinas[pool_meta[pid]["accepted_name"]]["earliest_year"] for pid in eligible_ids]
    scales = {
        "breadth": sd(breadth_vals),
        "records": sd(record_vals),
        "dated": sd(dated_vals),
        "year": sd(year_vals),
    }

    def features(pid):
        m = pool_meta[pid]
        s = sinas[m["accepted_name"]]
        return (
            math.log1p(m["introduced_wgsrpd3_count"]),
            math.log1p(s["records"]),
            math.log1p(s["dated_records"]),
            float(s["earliest_year"]),
        )

    def distance(a, b):
        fa, fb = features(a), features(b)
        return math.sqrt(
            ((fa[0] - fb[0]) / scales["breadth"]) ** 2
            + ((fa[1] - fb[1]) / scales["records"]) ** 2
            + ((fa[2] - fb[2]) / scales["dated"]) ** 2
            + ((fa[3] - fb[3]) / scales["year"]) ** 2
        )

    out = []
    diagnostics = []
    skipped = []
    for ev in events:
        sp = str(ev["species"]).strip()
        region = str(ev["wgsrpd3_code"]).strip()
        actual_id = str(ev["host_id"]).strip()
        actual_name = str(ev["host_name"]).strip()
        if actual_id not in pool_meta:
            skipped.append({
                "species": sp, "region": region, "actual_host_id": actual_id,
                "actual_host_name": actual_name, "reason": "ACTUAL_HOST_NOT_IN_WCVP_INTRODUCED_POOL"
            })
            continue
        if pool_meta[actual_id]["accepted_name"] not in sinas or sinas[pool_meta[actual_id]["accepted_name"]]["earliest_year"] is None:
            skipped.append({
                "species": sp, "region": region, "actual_host_id": actual_id,
                "actual_host_name": actual_name, "reason": "ACTUAL_HOST_LACKS_SINAS_DATED_CHRONOLOGY"
            })
            continue

        group = f"{sp}|{region}|{actual_id}"
        actual_s = sinas[pool_meta[actual_id]["accepted_name"]]
        out.append({
            **ev,
            "match_group": group,
            "assignment_type": "actual",
            "candidate_plant_id": actual_id,
            "candidate_name": pool_meta[actual_id]["accepted_name"],
            "candidate_family": pool_meta[actual_id]["family"],
            "match_rank": 0,
            "match_distance": 0.0,
            "introduced_wgsrpd3_count": pool_meta[actual_id]["introduced_wgsrpd3_count"],
            "sinas_records": actual_s["records"],
            "sinas_dated_records": actual_s["dated_records"],
            "sinas_earliest_global_year": actual_s["earliest_year"],
        })

        candidates = []
        for pid in sorted(set(pool_by_region.get(region, []))):
            if pid == actual_id or pid in focal_hosts.get(sp, set()):
                continue
            meta = pool_meta.get(pid)
            if meta is None:
                continue
            s = sinas.get(meta["accepted_name"])
            if not s or s["earliest_year"] is None:
                continue
            candidates.append((distance(actual_id, pid), pid))
        candidates.sort(key=lambda x: (x[0], pool_meta[x[1]]["accepted_name"], x[1]))
        chosen = candidates[: args.matches_per_host]
        diagnostics.append({
            "species": sp,
            "region": region,
            "actual_host_id": actual_id,
            "actual_host_name": actual_name,
            "candidate_pool": len(candidates),
            "matches_selected": len(chosen),
            "best_distance": None if not chosen else chosen[0][0],
            "worst_selected_distance": None if not chosen else chosen[-1][0],
        })
        for rank, (dist, pid) in enumerate(chosen, start=1):
            meta = pool_meta[pid]
            s = sinas[meta["accepted_name"]]
            out.append({
                **ev,
                "match_group": group,
                "assignment_type": "pseudo",
                "candidate_plant_id": pid,
                "candidate_name": meta["accepted_name"],
                "candidate_family": meta["family"],
                "match_rank": rank,
                "match_distance": dist,
                "introduced_wgsrpd3_count": meta["introduced_wgsrpd3_count"],
                "sinas_records": s["records"],
                "sinas_dated_records": s["dated_records"],
                "sinas_earliest_global_year": s["earliest_year"],
            })

    if not out:
        raise RuntimeError("no pseudo-host assignments")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out[0]))
        writer.writeheader()
        writer.writerows(out)

    payload = {
        "schema": "chocho_butterfly_temporal_pseudohost_matches_v0.1",
        "status": "POSTHOC_IDENTITY_SPECIFIC_TEMPORAL_NULL_MATCHES",
        "event_host_pairs": len(events),
        "input_actual_host_region_rows": len(events),
        "match_groups": len(diagnostics),
        "skipped_actual_host_region_rows": len(skipped),
        "skipped_reason_counts": {
            reason: sum(x["reason"] == reason for x in skipped)
            for reason in sorted({x["reason"] for x in skipped})
        },
        "matches_per_host_requested": args.matches_per_host,
        "all_plant_introduced_pool_species": len(pool_meta),
        "wcvp_sinas_eligible_pool_species": len(eligible_ids),
        "groups_with_full_requested_matches": sum(d["matches_selected"] == args.matches_per_host for d in diagnostics),
        "minimum_candidate_pool": None if not diagnostics else min(d["candidate_pool"] for d in diagnostics),
        "median_candidate_pool": None if not diagnostics else sorted(d["candidate_pool"] for d in diagnostics)[len(diagnostics)//2],
        "maximum_selected_match_distance": max(
            (d["worst_selected_distance"] for d in diagnostics if d["worst_selected_distance"] is not None),
            default=None
        ),
        "matching_features": [
            "log1p WCVP introduced WGSRPD3 breadth",
            "log1p SInAS alien-region record count",
            "log1p SInAS dated-record count",
            "SInAS earliest global alien-record year",
        ],
        "exclusions": [
            "all documented known hosts of the focal butterfly",
            "plants not introduced in the focal WGSRPD3 region",
            "plants without SInAS chronology",
        ],
        "diagnostics": diagnostics,
        "skipped": skipped,
        "claim_boundary": "These matches create a conditional post-hoc identity/timing null for the nine event cells. They do not make the event-cell sample response-blind.",
    }
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

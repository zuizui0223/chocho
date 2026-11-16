#!/usr/bin/env python3
"""Reverse-time historical record-index benchmark for the host/butterfly bridge.

The focal host record index is from 2001-2010 and contemporary regional WCVP
origins. This historical benchmark estimates its association with 1700-2000
butterfly GBIF indexing on the SAME original butterfly/island population.
It does not prove a causal negative control: native hosts may have existed
before 2001, and historical record retrieval is strongly selective.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from audit_native_introduced_temporal_bridge_20261116 import (
    CATEGORIES, make_panel, model,
)


def indexed_counts(rows: list[dict], outcome: str) -> dict:
    out = {}
    for category in CATEGORIES:
        subset = [r for r in rows if r["category"] == category]
        out[category] = {
            "risk": len(subset),
            "observed_index": sum(r[outcome] for r in subset),
            "fraction": (sum(r[outcome] for r in subset) / len(subset)
                         if subset else None),
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--history-csv", type=Path, required=True)
    parser.add_argument("--history-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path,
                        default=Path("results/native_introduced_bridge_20261116/"
                                     "pre2001_reverse_time_host_control.json"))
    args = parser.parse_args()
    rows, receipts = make_panel(args.repo_root, args.history_csv, args.history_summary)
    enriched = []
    for old in rows:
        row = dict(old)
        row["pre2001_flag"] = int(row["prior_1700_2000_count"] > 0)
        row["butterfly_2001_2010_flag"] = int(row["prior_2001_2010_count"] > 0)
        row["butterfly_2011_2022_flag"] = int(row["post"] > 0)
        enriched.append(row)

    audit = {}
    for minimum in (10, 30):
        subset = [r for r in enriched if r["effort"] >= minimum]
        if len(subset) != (8100 if minimum == 10 else 5850):
            raise ValueError("Original source botanical-coverage denominator changed")
        scope = {}
        for name, response in [
            ("reverse_time_1700_2000", "pre2001_flag"),
            ("same_era_2001_2010", "butterfly_2001_2010_flag"),
            ("later_records_2011_2022", "butterfly_2011_2022_flag"),
        ]:
            scope[name] = {
                "eligible_rows": len(subset),
                "index_count": sum(r[response] for r in subset),
                "host_classes": indexed_counts(subset, response),
                "model": model(subset, response),
            }
        # Within each identical species-island risk cell, compare later
        # source indexing with the earlier index indicators. Because the
        # observation windows and GBIF effort differ, this is only a paired
        # descriptive contrast of source associations, not a time effect.
        scope["paired_within_cell_index_differences"] = {}
        for short_name, old in (
            ("later_minus_1700_2000", "pre2001_flag"),
            ("later_minus_2001_2010", "butterfly_2001_2010_flag"),
        ):
            paired = [
                dict(row, index_difference=(
                    row["butterfly_2011_2022_flag"] - row[old]))
                for row in subset
            ]
            same_model = model(paired, "index_difference")
            same_model["net_index_count"] = same_model.pop("events")
            scope["paired_within_cell_index_differences"][short_name] = {
                "model": same_model,
                "n": len(paired),
                "later_yes_earlier_no": sum(
                    r["index_difference"] == 1 for r in paired
                ),
                "earlier_yes_later_no": sum(
                    r["index_difference"] == -1 for r in paired
                ),
                "same_index_status": sum(
                    r["index_difference"] == 0 for r in paired
                ),
                "interpretation": "difference in source-index associations; time windows not exchangeable",
            }
        audit[f"historical_plant_records_ge_{minimum}"] = scope

    result = {
        "schema": "chocho_reverse_time_host_record_benchmark_20261116_v1",
        "status": "SOURCE_INDEX_BENCHMARK_NOT_CAUSAL_NEGATIVE_CONTROL",
        "same_original_source": receipts,
        "exposure": "2001-2010 indexed compatible host near WPI port; native/introduced origin assigned from contemporary regional WCVP",
        "response": "Observed butterfly source-index presence in 1700-2000, 2001-2010, and 2011-2022, not actual historical first arrivals",
        "windows_identical_species_island_denominator": True,
        "paired_index_difference_modelled_within_original_225x51_species_islands": True,
        "results": audit,
        "limits": [
            "Plant host record is from 2001–2010 and may not be an ecological introduction or a new resource. The native host can exist before 2001.",
            "Reverse-time association detects longstanding geography or sampling selection, but DOES NOT falsify a true stable natural host requirement.",
            "Butterfly record probability differs greatly across the three historical time windows, so effect-size comparisons are descriptive.",
            "Paired index-difference effects are not interpretable as changes in ecological colonization because periods differ in duration, digitization and detection effort.",
            "Prior-date GBIF index zero is unknown detection, NOT butterfly absence.",
            "The classification none means no catalog-linked host indexed in these port buffers, NOT a truly missing host on an island.",
            "Species and island fixed effects do not adjust species-by-island historical sampling effort.",
            "These overlapping record windows are NOT independent test cohorts; CIs are approximate two-way clustered, not trial intervals.",
            "Neither historical index, reindex nor introduced-plant status identifies actual immigration, local larval feeding, establishment, abundance or plant arrival date.",
        ],
    }
    output = args.output if args.output.is_absolute() else args.repo_root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    compact = {threshold: {
        period: {
            "events": record["index_count"],
            "categories": {
                category: [value["observed_index"], value["risk"]]
                for category, value in record["host_classes"].items()
            },
            "adjusted_pp": {
                c: record["model"]["coefficients"][c]["percentage_points"]
                for c in ("introduced_only", "native_only", "both")
            },
            "approx95_pp": {
                c: record["model"]["coefficients"][c]["approximate_95pct_pp"]
                for c in ("introduced_only", "native_only", "both")
            },
        } for period, record in scopes.items()
        if period != "paired_within_cell_index_differences"
    } for threshold, scopes in audit.items()}
    print("REVERSE_TIME_SOURCE_BENCHMARK", json.dumps(compact, ensure_ascii=False),
          flush=True)
    print("REVERSE_TIME_PAIRED_SOURCE_ASSOCIATION",
          json.dumps({
              scope: {
                  name: {
                      "n": val["n"],
                      "later_yes_earlier_no": val["later_yes_earlier_no"],
                      "earlier_yes_later_no": val["earlier_yes_later_no"],
                      "native_only": val["model"]["coefficients"]["native_only"],
                      "introduced_only": val["model"]["coefficients"]["introduced_only"],
                      "both": val["model"]["coefficients"]["both"],
                  } for name, val in group["paired_within_cell_index_differences"].items()
              } for scope, group in audit.items()
          }, ensure_ascii=False), flush=True)
    print("REVERSE_TIME_OUTPUT", str(output), flush=True)


if __name__ == "__main__":
    main()

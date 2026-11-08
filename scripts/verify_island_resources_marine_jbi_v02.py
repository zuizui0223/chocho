#!/usr/bin/env python3
"""Audit integrated JBI v0.2 claims against committed source files.

This is a *reference/numeric consistency* gate, not a fresh ecological model fit.
It never changes source data or manuscripts and never downloads anything.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


NUMBER = re.compile(r"(?<![A-Za-z_])[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[Ee][-+]?\d+)?")
BASE = Path(__file__).resolve().parents[1]
FOLDER = BASE / "manuscript/island_resources_marine_access_20261108"

# These quantities are recorded in the successful workflow artifact rather
# than the committed pre-analysis protocol; they cannot be verified from Git
# source files alone. Do not claim these have been verified by this script.
EXTERNAL_CURRENT_VALUES = {"max_grounded_residual"}


def numbers_in(obj):
    if isinstance(obj, bool) or obj is None:
        return []
    if isinstance(obj, (int, float)):
        return [float(obj)] if math.isfinite(obj) else []
    if isinstance(obj, dict):
        values = []
        for item in obj.values():
            values.extend(numbers_in(item))
        return values
    if isinstance(obj, (list, tuple)):
        values = []
        for item in obj:
            values.extend(numbers_in(item))
        return values
    if isinstance(obj, str):
        out = []
        for match in NUMBER.finditer(obj.replace("−", "-")):
            try:
                num = float(match.group(0).replace(",", ""))
                if math.isfinite(num):
                    out.append(num)
            except ValueError:
                continue
        return out
    return []


def load_numbers(path: Path):
    content = path.read_text(encoding="utf-8")
    if path.suffix == ".json":
        return numbers_in(json.loads(content))
    return numbers_in(content)


def matched(expected, observed):
    # Commits store some values in 4-decimal prose and others at full precision.
    # A 0.00011 tolerance accepts 4-decimal rounding but NOT a 1pp error.
    tolerance = max(0.00011, abs(expected) * 0.000012)
    return any(math.isclose(expected, actual, rel_tol=0.0, abs_tol=tolerance)
               for actual in observed)


def audit():
    ledger = json.loads((FOLDER / "EVIDENCE_LEDGER_V02.json").read_text(encoding="utf-8"))
    manuscript = (FOLDER / "ISLAND_RESOURCES_MARINE_ACCESS_MANUSCRIPT_V02.md").read_text(encoding="utf-8")
    figures = (FOLDER / "FIGURES_CAPTIONS_V02.md").read_text(encoding="utf-8")
    assert ledger["schema"] == "chocho_JBI_island_resource_marine_v02_evidence_ledger"
    assert len(ledger["evidence"]) == 16, "Source ledger was unexpectedly changed"
    assert "239 butterfly" in manuscript and "225" in manuscript
    assert "1,275" in manuscript and "51" in manuscript
    assert "not submitted" in manuscript.lower()
    assert "not a demonstrated ship-mediated" in manuscript.lower()
    assert "10.7" in manuscript and "10.7" in figures
    assert "+1.7312%" in manuscript and "+6.4880%" in manuscript
    assert "100/199" in manuscript and "5/199" in manuscript
    assert "−2.9233%" in manuscript and "+8.1461%" in manuscript
    assert "1.4_percent_as_biological_realization_rate" in ledger["claims_excluded_from_manuscript"]
    assert "14_percentage_point_island_effect" in ledger["claims_excluded_from_manuscript"]
    assert "13 of 920" in manuscript and "1.413%" in manuscript

    failures, exceptions, checked = [], [], 0
    for entry in ledger["evidence"]:
        source = BASE / entry["source"]
        if not source.is_file():
            failures.append((entry["id"], "missing source", entry["source"]))
            continue
        nums = load_numbers(source)
        auxiliary = entry.get("supplemental_source", "")
        if auxiliary and not auxiliary.startswith("GitHub Actions"):
            other = BASE / auxiliary
            if not other.is_file():
                failures.append((entry["id"], "missing auxiliary source", auxiliary))
                continue
            nums.extend(load_numbers(other))
        for key, value in entry["values"].items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                continue
            if entry["id"] == "marine_current_map" and key in EXTERNAL_CURRENT_VALUES:
                exceptions.append((entry["id"], key, "separate successful workflow artifact"))
                continue
            if not matched(float(value), nums):
                failures.append((entry["id"], key, value))
            else:
                checked += 1

    report = {
        "status": "PASS" if not failures else "FAIL",
        "evidence_entries": len(ledger["evidence"]),
        "numeric_assertions_found_in_committed_sources": checked,
        "numeric_assertions_in_external_artifacts_not_rechecked_here": exceptions,
        "failures": failures,
        "claim_status": "source/numeric consistency only; neither ecological refit nor journal submission"
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if failures:
        raise SystemExit(1)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=BASE,
                        help="Unused compatibility argument: only the checkout rooted at this script is audited")
    parser.parse_args()
    audit()


if __name__ == "__main__":
    main()

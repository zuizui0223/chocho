#!/usr/bin/env python3
"""Frozen country-panel overlap audit; no island/port effect is estimated.

Replay from repo root: python scripts/analyze_island_host_identifiability_20261029.py
The geographic class is a conservative country screen, not measured landmass
identity. Zero indexed detections are retained as source recording outcomes,
never interpreted as biological absences or establishment failures.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_gbif_country_record_turnover_20261020 import read_sources, states_for

ROOT = Path(__file__).resolve().parents[1]
FROZEN = "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json"
COUNTS = "data/exploratory/gbif_country_20261020/gbif_country_counts.csv"
MANIFEST = "data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json"
ISLAND_COUNTRIES = frozenset("""
AG BS BB BN CU CV CY DM DO FJ GB GD HT ID IE IS JM JP KM KN LC LK
MG MT MU MV NC NZ PG PH SB SC SG ST TL TO TT TW VC VU WS
""".split())
CLASSES = ("wholly_island_country", "continental_or_mixed_country")


def country_class(code):
    return CLASSES[0] if code in ISLAND_COUNTRIES else CLASSES[1]


def directional_rows(rows, direction):
    """Use original temporal risk-set definitions, without changing exposures."""
    if direction not in ("late_only", "early_only"):
        raise ValueError("Expected late_only or early_only")
    prior, response = ("early", "late") if direction == "late_only" else ("late", "early")
    return [{**r, "cells": [{**c, "positive": c[response]} for c in r["cells"]
                              if c[prior] == 0]} for r in rows]


def analyze_support(rows, cohort):
    """Count exposure/outcome overlap within preserved species×Level1 groups."""
    margins = {k: Counter() for k in CLASSES}
    support = {k: Counter() for k in CLASSES}
    four = Counter()
    country = {}
    stratum_rows = []
    all_cells = []
    for row in rows:
        cells = row["cells"]
        all_cells.extend(cells)
        if any(c["added"] not in (0, 1) or c["positive"] not in (0, 1) for c in cells):
            raise ValueError("Exposure and source recording outcome must be binary")
        counts = {}
        for cls in CLASSES:
            selected = [c for c in cells if country_class(c["country"]) == cls]
            n = len(selected)
            a = sum(c["added"] for c in selected)
            k = sum(c["positive"] for c in selected)
            contrast = 0 < a < n
            informative = contrast and 0 < k < n
            margins[cls].update(cells=n, added_cells=a, positive_cells=k,
                                added_positive_cells=sum(c["added"] * c["positive"] for c in selected))
            support[cls].update(host_contrast_strata=int(contrast),
                                cells_in_host_contrast_strata=n if contrast else 0,
                                positive_cells_in_host_contrast_strata=k if contrast else 0,
                                informative_strata=int(informative))
            for c in selected:
                ct = country.setdefault(c["country"], Counter())
                ct.update(cells=1, added_cells=c["added"], positive_cells=c["positive"],
                          added_positive_cells=c["added"] * c["positive"],
                          cells_in_within_class_host_contrast_strata=int(contrast))
            counts[cls] = (n, a, k)
            stratum_rows.append({"cohort": cohort, "species": row["species"],
                                 "genus": row["genus"], "level1": row["stratum"],
                                 "country_class": cls, "cells": n, "added_cells": a,
                                 "positive_cells": k, "host_contrast": int(contrast),
                                 "informative": int(informative)})
        full = all(0 < a < n for n, a, _ in counts.values())
        n = len(cells)
        k = sum(c["positive"] for c in cells)
        if full:
            four.update(strata=1, cells=n, positive_cells=k,
                        island_positive_cells=counts[CLASSES[0]][2],
                        informative_strata=int(0 < k < n))
    # Supply explicit zeros rather than omitted keys when no support exists.
    for cls in CLASSES:
        for key in ("cells", "added_cells", "positive_cells", "added_positive_cells"):
            margins[cls].setdefault(key, 0)
        for key in ("host_contrast_strata", "cells_in_host_contrast_strata",
                    "positive_cells_in_host_contrast_strata", "informative_strata"):
            support[cls].setdefault(key, 0)
    for key in ("strata", "cells", "positive_cells", "island_positive_cells", "informative_strata"):
        four.setdefault(key, 0)
    country_rows = [{"cohort": cohort, "country": ct,
                     "country_class": country_class(ct), **dict(values)}
                    for ct, values in sorted(country.items())]
    summary = {"strata_preserved": len(rows), "candidate_cells": len(all_cells),
               "eligible_countries": sorted(country),
               "descriptive_margins_confounded": {k: dict(v) for k, v in margins.items()},
               "within_class_host_support": {k: dict(v) for k, v in support.items()},
               "all_four_groups_support": dict(four)}
    return summary, country_rows, stratum_rows


def run(root, output):
    strata, counts, background, receipts = read_sources(root / COUNTS, root / MANIFEST, root / FROZEN)
    if any(c["other_species_richness"] < 3 for r in strata for c in r["cells"]):
        raise ValueError("Frozen primary >=3 other-taxa pooled effort gate violated")
    paired, paired_info = states_for(strata, counts, background, "all_accepted_bases", 3, 1)
    paired_pooled = [{**r, "cells": [{**c, "positive": int(c["early"] or c["late"])}
                                      for c in r["cells"]]} for r in paired]
    cohorts = {"frozen_primary_pooled": strata,
               "paired_other3_pooled": paired_pooled,
               "paired_other3_late_only_risk_set": directional_rows(paired, "late_only"),
               "paired_other3_early_only_risk_set": directional_rows(paired, "early_only")}
    results, country_rows, stratum_rows = {}, [], []
    for name, rows in cohorts.items():
        result, cr, sr = analyze_support(rows, name)
        results[name] = result
        country_rows.extend(cr)
        stratum_rows.extend(sr)
    state_counts = {cls: dict(sorted(Counter(c["state"] for r in paired for c in r["cells"]
                                         if country_class(c["country"]) == cls).items()))
                    for cls in CLASSES}
    unsupported = all(
        r["within_class_host_support"][CLASSES[0]]["positive_cells_in_host_contrast_strata"] == 0
        and r["all_four_groups_support"]["positive_cells"] == 0
        for r in results.values())
    result = {
        "schema": "chocho_frozen_country_island_host_identifiability_20261029_v1",
        "status": ("NO_POSITIVE_OUTCOME_SUPPORT_IN_ISLAND_HOST_CONTRAST_OR_FOUR_GROUP_STRATA"
                   if unsupported else "POSITIVE_SUPPORT_PRESENT_REQUIRES_FURTHER_IDENTIFICATION_AUDIT"),
        "source_files": {"frozen_primary_strata": FROZEN, "country_counts": COUNTS,
                         "original_taxon_manifest": MANIFEST},
        "source_receipts": receipts,
        "country_classification": {
            "method": "Conservative fixed political-country screen; not a landmass geometry join",
            "wholly_island_country_screen": sorted(ISLAND_COUNTRIES),
            "eligible_wholly_island_countries": sorted(set(results["frozen_primary_pooled"]["eligible_countries"]) & ISLAND_COUNTRIES),
            "complement_label": CLASSES[1],
            "limitations": "Countries may share islands or include multiple islands. Classification is not island independence, mainland connectivity, port exposure, isolation or a definitive global country ontology."},
        "definitions": {
            "primary": "Preserve frozen species×Level1 strata, native-host exclusion, >=90% climate analogue, >=3 other source taxa pooled, and >=1 focal indexed record outcome.",
            "temporal": "Same original candidates; >=3 OTHER TAXA detected in EACH original window, not visits. Windows 2010–2017 and 2018–2025. Late-only risk set excludes early detections; early-only excludes late detections.",
            "informative_stratum": "Both host classes within selected geography and 0<K<N indexed-record outcomes.",
            "four_groups": "Both host classes in BOTH country geography classes within the same species×Level1 stratum; informative additionally requires 0<K<N across those cells.",
            "zero": "Zero indexed detections under complete source queries; never biological absence, invasion failure or survey nondetection."},
        "cohorts": results,
        "paired_temporal_source_gate": paired_info,
        "paired_temporal_state_counts_by_country_class": state_counts,
        "interpretation": "Marginal geography tables are composition-confounded. No pooled, causal, colonization or port effect estimated. Source overlap audit only; changing model assumptions does not create missing within-stratum positive support.",
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "island_host_support_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for name, rows in (("island_country_support.csv", country_rows), ("island_stratum_support.csv", stratum_rows)):
        with (output / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps({"status": result["status"], "primary": results["frozen_primary_pooled"]}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/exploratory/island_port_20261029")
    args = parser.parse_args()
    run(args.root, args.output_dir)

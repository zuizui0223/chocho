#!/usr/bin/env python3
"""Country-grain host-climate niche-entry opportunity vs GBIF recording, source-only.

This is an ecological *opportunity* comparison, not evidence of butterfly
colonisation, range dynamics, realised climatic niche, or local host use.
See frozen protocol BUTTERFLY_NICHE_ENTRY_COUNTRY_20261026.md.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from analyze_gbif_country_host_opportunity_20261020 import load_observations

SEED = 20261026
BLOBS = {
    "areas": "6c104dfeeb4796d8a4e149772750952b68d2f763",
    "design": "452e85f4adfeb264054076a2b5a57c7963cad332",
    "counts": "10fe5b3f89d510b9c4d19dd94124ed51ee3b2c31",
    "manifest": "c6869cb14bc221d88c5c88781f1f18671f608492",
    "receipt": "1855e4332fb5eae5afee65524bf4f32dcc117bdb",
}
PRIMARY = {
    "basis": "all_accepted_bases", "window": "pooled", "effort_gate": 3,
    "min_area": .05, "purity": .80, "min_records": 1,
}
MIN_INFORMATIVE = 50
MIN_SPECIES = 25
MIN_POSITIVE = 25


def blob_digest(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\\0" + data).hexdigest()


def check_source(path, name):
    data = Path(path).read_bytes()
    digest = blob_digest(data)
    if digest != BLOBS[name]:
        raise ValueError(f"Changed original {name} Git blob: {digest}")
    return data


def parse_bool(value):
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"Unexpected source boolean {value!r}")


def opportunity_class(added, aligned, *, min_area=.05, purity=.80):
    """Classify source potential area, not butterfly habitat or presence."""
    if not all(math.isfinite(z) for z in (added, aligned)):
        raise ValueError("Nonfinite country host area")
    if not 0 <= aligned <= added + 1e-8 or not 0 <= added <= 1 + 1e-8:
        raise ValueError("Impossible source country coverage fractions")
    if not (0 < min_area <= 1) or not (.5 < purity <= 1):
        raise ValueError("Invalid fixed area/purity gate")
    if added < min_area - 1e-10:
        return "low_or_zero_added"
    analogous_share = min(1., max(0., aligned / added))
    if analogous_share >= purity - 1e-10:
        return "analogue"
    if 1 - analogous_share >= purity - 1e-10:
        return "novel"
    return "mixed"


def load_areas(areas_path, design_path):
    design = json.loads(check_source(design_path, "design").decode("utf-8"))
    raw = check_source(areas_path, "areas")
    rows = []
    seen = set()
    with gzip.open(areas_path, "rt", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {
            "species", "country", "Level1", "country_area_km2",
            "native_guard_present", "added_fraction", "climate_fraction",
            "added_climate_fraction", "added_country",
        }
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Source country exposure schema changed")
        for r in reader:
            key = (r["species"], r["country"])
            if key in seen:
                raise ValueError("Repeated species/country area in frozen source")
            seen.add(key)
            area = float(r["country_area_km2"])
            added = float(r["added_fraction"])
            climate = float(r["climate_fraction"])
            aligned = float(r["added_climate_fraction"])
            if not area > 0 or not 0 <= climate <= 1 + 1e-8:
                raise ValueError("Invalid country area or source climate coverage")
            opportunity_class(added, aligned)
            if aligned > climate + 1e-8:
                raise ValueError("Intersection exceeds native climate area")
            native = parse_bool(r["native_guard_present"])
            added_country = parse_bool(r["added_country"])
            if native and added_country:
                raise ValueError("Country labelled both native-guard and new host")
            if not native and (added > 1e-10) != added_country:
                raise ValueError("Introduced-host flag disagrees with mapped added area")
            rows.append({
                "species": r["species"], "country": r["country"],
                "level1": str(r["Level1"]), "country_km2": area,
                "native": native, "added": added, "analog": aligned,
                "novel": max(0., added - aligned),
            })
    countries = set(design["geography"]["eligible_countries"])
    species = {x["species"] for x in rows}
    if len(countries) != 140 or len(species) != 235:
        raise ValueError("Source country or 5-95% butterfly cohort changed")
    if len(rows) != len(species) * len(countries):
        raise ValueError("Source species/country grid is not complete")
    if {x["country"] for x in rows} != countries:
        raise ValueError("Country source roster changed")
    if design["host_climate_eligible_species"] != len(species):
        raise ValueError("Country design source cohort mismatch")
    return rows, {
        "source_area_gzip_bytes": len(raw),
        "source_species": len(species), "source_countries": len(countries),
        "source_species_country_cells": len(rows),
        "source_country_area_definition": "full country area, disjoint geographic atoms",
    }


def geography_summary(rows):
    """Area-weighted available host-resource geography, not insect density."""
    out = Counter()
    by_region = defaultdict(Counter)
    for r in rows:
        out["source_cells_all"] += 1
        if r["native"]:
            out["source_cells_native_guard_excluded"] += 1
            continue
        out["source_cells_outside_native"] += 1
        out["added_species_km2"] += r["country_km2"] * r["added"]
        out["analog_species_km2"] += r["country_km2"] * r["analog"]
        out["novel_species_km2"] += r["country_km2"] * r["novel"]
        kind = opportunity_class(r["added"], r["analog"])
        out["cells_" + kind] += 1
        by_region[r["level1"]]["analog_species_km2"] += r["country_km2"] * r["analog"]
        by_region[r["level1"]]["novel_species_km2"] += r["country_km2"] * r["novel"]
        by_region[r["level1"]]["cells_" + kind] += 1
    if out["added_species_km2"] <= 0:
        raise ValueError("No source added host resource area")
    if not math.isclose(
        out["added_species_km2"],
        out["analog_species_km2"] + out["novel_species_km2"],
        rel_tol=1e-9, abs_tol=1e-3,
    ):
        raise ValueError("Potential niche opportunity area fails conservation")
    return {
        **dict(out),
        "novel_share_area_weighted": out["novel_species_km2"] / out["added_species_km2"],
        "level1": {k: dict(v) for k, v in sorted(by_region.items())},
    }


def exact_convolution(strata, observed):
    """Artificial fixed-margin hypergeometric geography reference, exact DP."""
    pmf = {0: 1.0}
    for row in strata:
        n, a, k = row["N"], row["A"], row["K"]
        lo, hi = max(0, k - (n - a)), min(a, k)
        denom = math.comb(n, k)
        dist = {
            j: math.comb(a, j) * math.comb(n - a, k - j) / denom
            for j in range(lo, hi + 1)
        }
        next_pmf = defaultdict(float)
        for x, px in pmf.items():
            for j, pj in dist.items():
                next_pmf[x + j] += px * pj
        pmf = dict(next_pmf)
    if abs(sum(pmf.values()) - 1) > 1e-9:
        raise ValueError("Exact artificial geography null normalization failed")
    return {
        "right_tail_not_biological_p": sum(v for k, v in pmf.items() if k >= observed),
        "left_tail_not_biological_p": sum(v for k, v in pmf.items() if k <= observed),
        "conditional_support": [min(pmf), max(pmf)] if pmf else [0, 0],
    }


def genus_bootstrap(strata, *, draws=1999):
    groups = defaultdict(list)
    for r in strata:
        groups[r["genus"]].append(r)
    names = sorted(groups)
    if not names:
        return None
    matrix = np.array([
        [sum(x[z] for x in groups[g]) for z in ("O", "expected", "K")]
        for g in names
    ])
    rng = np.random.default_rng(SEED)
    picks = rng.integers(0, len(names), size=(draws, len(names)))
    totals = matrix[picks].sum(axis=1)
    valid = totals[:, 2] > 0
    if valid.sum() < .9 * draws:
        return {"status": "NOT_ENOUGH_POSITIVE_GENUS_BOOTSTRAP_DRAWS", "valid": int(valid.sum())}
    delta = 100 * (totals[valid, 0] - totals[valid, 1]) / totals[valid, 2]
    return {
        "status": "DESCRIPTIVE_DEPENDENT_SOURCE_BOOTSTRAP",
        "distinct_genera": len(names), "valid_draws": int(valid.sum()),
        "95pct_percentile_pp": [float(np.quantile(delta, q)) for q in (.025, .975)],
    }


def compare(rows, counts, complete, *, basis, window, effort_gate, min_area,
            purity, min_records):
    known = set(complete[basis])
    observed = counts[basis, window]
    effort = defaultdict(set)
    for (species, country), n in observed.items():
        if species in known and n > 0:
            effort[country].add(species)
    disposition = Counter()
    candidates = defaultdict(list)
    for r in rows:
        sp, country = r["species"], r["country"]
        if sp not in known:
            disposition["unresolved_or_unqualified_butterfly"] += 1
            continue
        if r["native"]:
            disposition["native_host_country_excluded"] += 1
            continue
        kind = opportunity_class(r["added"], r["analog"], min_area=min_area, purity=purity)
        disposition["source_" + kind] += 1
        if kind not in {"novel", "analogue"}:
            continue
        effort_others = len(effort[country] - {sp})
        if effort_others < effort_gate:
            disposition["excluded_insufficient_other_butterfly_records"] += 1
            continue
        positive = int(observed.get((sp, country), 0) >= min_records)
        disposition["matched_candidate_" + kind] += 1
        disposition["matched_positive_" + kind] += positive
        candidates[(sp, r["level1"])].append((kind, positive))
    strata = []
    for (sp, level1), cells in sorted(candidates.items()):
        n = len(cells)
        a = sum(kind == "novel" for kind, _ in cells)
        k = sum(pos for _, pos in cells)
        o = sum(kind == "novel" and pos for kind, pos in cells)
        if not 0 < a < n:
            disposition["discarded_strata_not_both_climates"] += 1
            continue
        strata.append({
            "species": sp, "genus": sp.split()[0], "level1": level1,
            "N": n, "A": a, "K": k, "O": o, "expected": a * k / n,
        })
    N = sum(x["N"] for x in strata)
    A = sum(x["A"] for x in strata)
    K = sum(x["K"] for x in strata)
    O = sum(x["O"] for x in strata)
    expectation = sum(x["expected"] for x in strata)
    informative = sum(0 < x["K"] < x["N"] for x in strata)
    species_positive = {x["species"] for x in strata if x["K"] > 0}
    passed = (
        informative >= MIN_INFORMATIVE
        and len(species_positive) >= MIN_SPECIES
        and K >= MIN_POSITIVE
    )
    return {
        "settings": dict(
            basis=basis, window=window, effort_gate=effort_gate, min_area=min_area,
            purity=purity, min_records=min_records
        ),
        "status": "EXPLORATORY_INFORMATION_GATE_MET" if passed else "INSUFFICIENT_NICHE_ENTRY_INFORMATION",
        "information_thresholds": {
            "informative_species_level1_strata_min": MIN_INFORMATIVE,
            "positive_species_min": MIN_SPECIES,
            "positive_country_species_cells_min": MIN_POSITIVE,
        },
        "source_diagnostics": dict(sorted(disposition.items())),
        "strata_with_both_resource_classes": len(strata),
        "informative_strata": informative,
        "positive_species": len(species_positive),
        "candidate_country_species_cells": N, "novel_candidate_cells": A,
        "analog_candidate_cells": N - A, "positive_country_species_cells": K,
        "positive_novel_cells": O, "positive_analog_cells": K - O,
        "expected_positive_novel_cells": expectation,
        "observed_novel_share_among_positive": O / K if K else None,
        "expected_novel_share_among_positive": expectation / K if K else None,
        "observed_minus_expected_novel_pp": 100 * (O - expectation) / K if K else None,
        "exact_conditional_reference": exact_convolution(strata, O),
        "genus_cluster_bootstrap": genus_bootstrap(strata),
        "interpretation_gate": (
            "MATHEMATICAL_COUNTRY_INDEX_ASSOCIATION_ONLY; no actual butterfly "
            "climatic niche entry, successful establishment, host use or colonisation"
        ),
    }


def run(args):
    files = {
        "areas": args.areas, "design": args.design, "counts": args.counts,
        "manifest": args.original_manifest, "receipt": args.country_receipt,
    }
    for name, path in files.items():
        check_source(path, name)
    rows, source_area = load_areas(args.areas, args.design)
    geo = geography_summary(rows)
    with tempfile.TemporaryDirectory(prefix="butterfly-niche-entry-") as temp:
        receipt_path = Path(temp) / "original_country_receipt.json"
        with gzip.open(args.country_receipt, "rt", encoding="utf-8") as src:
            receipt_path.write_text(src.read(), encoding="utf-8")
        counts, complete, old, receipt, aliases = load_observations(
            args.counts, receipt_path, args.original_manifest
        )
    if aliases != ["Oeneis bore", "Oeneis polixenes"]:
        raise ValueError("Accepted GBIF taxon alias set changed")
    if len(complete["all_accepted_bases"]) != 229 or len(complete["core"]) != 229:
        raise ValueError("Source country butterfly completeness gate changed")
    settings = {
        "primary": PRIMARY,
        "core": {**PRIMARY, "basis": "core"},
        "early_2010_2017": {**PRIMARY, "window": "2010_2017"},
        "late_2018_2025": {**PRIMARY, "window": "2018_2025"},
        "purity_90": {**PRIMARY, "purity": .90},
        "min_area_20pct": {**PRIMARY, "min_area": .20},
        "effort_1": {**PRIMARY, "effort_gate": 1},
        "effort_5": {**PRIMARY, "effort_gate": 5},
        "minimum_5_indexed_records": {**PRIMARY, "min_records": 5},
    }
    outputs = {name: compare(rows, counts, complete, **parameters)
               for name, parameters in settings.items()}
    out = {
        "schema": "chocho_host_climate_niche_entry_country_source_index_v01",
        "status": outputs["primary"]["status"],
        "source": {
            **source_area, "complete_species": 229,
            "original_taxon_aliases_excluded": aliases,
            "source_git_blob_sha1": BLOBS,
            "retained_source_index": "exact same GBIF full-index country facets as PR51",
        },
        "opportunity_geography": geo,
        "primary": outputs.pop("primary"),
        "one_factor_at_a_time_sensitivity": outputs,
        "inferential_red_lines": [
            "Source climate rectangle is native HOST-PLANT, not butterfly climatic tolerance.",
            "Species-country plant area and GBIF country presence do not locate butterfly in the plant subarea.",
            "No host introduction dates, surveyed nondetection, breeding or independent populations.",
            "Exact geography tails are artificial, not biological p values; genus bootstrap shares host and geography dependence.",
            "Original GEB paper remains unchanged; all claims require an information gate.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
    print(json.dumps({
        "status": out["status"],
        "opportunity_novel_share_country_area": geo["novel_share_area_weighted"],
        "source_outside_native_country_cells": geo["source_cells_outside_native"],
        "primary": {
            k: v for k, v in out["primary"].items()
            if k in (
                "status", "informative_strata", "positive_species",
                "positive_country_species_cells", "novel_candidate_cells",
                "analog_candidate_cells", "positive_novel_cells",
                "positive_analog_cells", "observed_minus_expected_novel_pp",
            )
        },
        "sensitivity_status": {k: v["status"] for k, v in outputs.items()},
    }, sort_keys=True))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--areas", type=Path, required=True)
    p.add_argument("--design", type=Path, required=True)
    p.add_argument("--counts", type=Path, required=True)
    p.add_argument("--original-manifest", type=Path, required=True)
    p.add_argument("--country-receipt", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    run(p.parse_args())

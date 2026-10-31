#!/usr/bin/env python3
"""Source-pinned UNCTAD historical *country-pair* Liner Shipping Bilateral Connectivity.

Reads the independently published UNCTAD 2026-06 full public country-pair CSV;
does not turn countries into ports, nor old plant observations into ship import.
The actual source is downloaded by the CI acquisition workflow only and is
NEVER added to this repository. Output is a small island-pair derived table and
a zero/unknown-preserving evidence audit.

This is observational geography. Only 2010 index values are primary. The
2006-2010 mean is a sensitivity where available, not a causal exposure.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pycountry

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "provenance/reviewer_defenses/results"
INPUTS = {
    "port_island": BASE / "island_dated_host_match_20261030/island_port_early_plant_availability.csv",
    "port_country": BASE / "port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv",
    "eligible_islands": BASE / "island_host_network_20261030/island_host_network_islands.csv",
    "island_pairs": BASE / "island_host_network_20261030/island_host_network_pairs.csv",
}
PINNED_BLOBS = {
    "port_island": "0a01831056c2cd59a4bd41bd1dbd4c84d77cf922",
    "port_country": "fdf73b11ddcb1644f46730123350ae52cc61cae6",
    "eligible_islands": "18cecde5bf3066f276380167f451f1fa5cc1caa1",
    "island_pairs": "5663d2ffe42419f5c4496b92c9f3bbf76175907b",
}
RAW_URL = "https://raw.githubusercontent.com/DeckardShaw31/Eureka-2026/main/data/raw/unctad/US_LSBCI.csv"
RAW_SHA256 = "edab41de33b9a6b4fc568bf454544dedfcbf32117e2bad7ccc77203e1ee6c33e"
PILOT_URL = "https://raw.githubusercontent.com/DeckardShaw31/Eureka-2026/main/data/processed/panel_bilateral_year.csv"
PILOT_BLOB = "b24aeae93acf7f987bea313b4eb99c99e68bdb65"

ALIASES = {
    "United Kingdom": "GBR", "UK": "GBR",
    "United Kingdom of Great Britain and Northern Ireland": "GBR",
    "China, Taiwan Province of": "TWN", "Taiwan Province of China": "TWN",
    "Taiwan, Province of China": "TWN",
    "China, Hong Kong SAR": "HKG", "Hong Kong, China": "HKG",
    "China, Macao SAR": "MAC", "Macao, China": "MAC",
    "Hong Kong Special Administrative Region of China": "HKG",
    "Macao Special Administrative Region of China": "MAC",
    "Russia": "RUS", "Republic of Korea": "KOR",
    "Democratic People's Republic of Korea": "PRK",
    "Viet Nam": "VNM", "Brunei Darussalam": "BRN",
    "Czechia": "CZE", "Czech Republic": "CZE",
    "Réunion": "REU", "Reunion": "REU",
    "United States": "USA", "United States of America": "USA",
    "New Caledonia": "NCL", "Christmas Island": "CXR",
    "Guam": "GUM",
}


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as file:
        return list(csv.DictReader(file))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def iso3(name):
    if not name:
        return None
    s = str(name).strip()
    if s in ALIASES:
        return ALIASES[s]
    try:
        country = pycountry.countries.lookup(s)
        return country.alpha_3
    except LookupError:
        return None


def validate_source(raw):
    actual = sha256(raw)
    if actual != RAW_SHA256:
        raise ValueError(
            "Full UNCTAD source SHA-256 changed: expected " + RAW_SHA256 +
            "; observed " + actual + ". Stop rather than silently use revisions."
        )
    return actual


def fixed_islands():
    eligible = {x["landmass_id"] for x in read_csv(INPUTS["eligible_islands"])
                if x["eligible"] == "True"}
    assert len(eligible) == 51
    cbyport = {x["port_id"]: x["country"] for x in read_csv(INPUTS["port_country"])}
    members = defaultdict(set)
    mapped_port_ids = defaultdict(set)
    for r in read_csv(INPUTS["port_island"]):
        island = r["landmass_id"]
        if island not in eligible:
            continue
        code = cbyport.get(r["port_id"])
        if code is None:
            raise ValueError("Missing source WPI port-country: " + r["port_id"])
        country = pycountry.countries.get(alpha_2=code)
        if not country:
            raise ValueError("Unknown ISO-2 in WPI catalogue: " + code)
        members[island].add(country.alpha_3)
        mapped_port_ids[island].add(r["port_id"])
    assert set(members) == eligible
    stable = {i: next(iter(v)) for i, v in members.items() if len(v) == 1}
    mixed = {i: sorted(v) for i, v in members.items() if len(v) > 1}
    assert len(stable) == 50 and len(mixed) == 1
    assert sum(map(len, mapped_port_ids.values())) == 194
    return stable, mixed


def source_bilateral(raw):
    """Average listed quarters/directions per country-pair/year, retaining zero."""
    values = defaultdict(list)
    unresolved = Counter()
    inspected = 0
    matched_rows = 0
    seen_years = Counter()
    with Path(raw).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        needed = {"Quarter", "Economy Label", "Partner Label", "Index"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError("UNCTAD columns changed: " + repr(reader.fieldnames))
        for r in reader:
            inspected += 1
            q = r["Quarter"].strip()
            if len(q) < 4 or not q[:4].isdigit():
                continue
            year = int(q[:4])
            if year < 2006 or year > 2010:
                continue
            seen_years[year] += 1
            a, b = iso3(r["Economy Label"]), iso3(r["Partner Label"])
            if a is None or b is None:
                unresolved[(r["Economy Label"] if a is None else r["Partner Label"])] += 1
                continue
            if a == b or not r["Index"].strip():
                continue
            try:
                value = float(r["Index"])
            except ValueError:
                continue
            if not (math.isfinite(value) and 0 <= value <= 1):
                raise ValueError("Bad LSBCI score " + str(value))
            key = tuple(sorted((a, b)))
            values[(key, year)].append(value)
            matched_rows += 1
    annual = {k: sum(v) / len(v) for k, v in values.items()}
    bypair = defaultdict(dict)
    for (key, year), val in annual.items():
        bypair[key][year] = val
    # Source was not expected to omit all European or non-Asian countries.
    if len(bypair) < 200:
        raise ValueError("Not a global bilateral panel: only " + str(len(bypair)) + " country pairs")
    return bypair, {
        "source_csv_rows": inspected, "year_counts_including_missing": dict(seen_years),
        "valid_2006_2010_source_rows": matched_rows,
        "source_country_pairs": len(bypair),
        "unknown_country_labels_top20": unresolved.most_common(20),
        "unknown_country_label_rows": sum(unresolved.values()),
    }


def aggregate(raw, pilot=None):
    stable, mixed = fixed_islands()
    bypair, audit = source_bilateral(raw)
    prior = read_csv(INPUTS["island_pairs"])
    assert len(prior) == 1275
    out = []
    for r in prior:
        a, b = r["island_a"], r["island_b"]
        ac, bc = stable.get(a), stable.get(b)
        pair = tuple(sorted((ac, bc))) if ac and bc else None
        years = bypair.get(pair, {}) if pair is not None else {}
        if not ac or not bc:
            status = "MULTICOUNTRY_ISLAND_UNRESOLVED"
        elif ac == bc:
            status = "SAME_COUNTRY_ISLAND_PAIR_UNRESOLVED"
        elif 2010 not in years:
            status = "MISSING_2010_BILATERAL_COUNTRY_INDEX"
        else:
            status = "INDEX_2010_OBSERVED"
        yr_2006_2010 = [years[y] for y in range(2006, 2011) if y in years]
        out.append({
            "island_a": a, "island_b": b,
            "country_a": ac or "", "country_b": bc or "",
            "status": status,
            "lsbci_2010": (years[2010] if 2010 in years and status == "INDEX_2010_OBSERVED" else ""),
            "lsbci_2006_2010_mean": (
                sum(yr_2006_2010) / len(yr_2006_2010)
                if len(yr_2006_2010) >= 3 and status == "INDEX_2010_OBSERVED" else ""),
            "lsbci_pre2011_observed_years": len(yr_2006_2010),
            "host_jaccard": r["host_jaccard"],
            "new_butterfly_jaccard": r["new_butterfly_jaccard"],
            "earlier_butterfly_jaccard": r["earlier_butterfly_jaccard"],
            "port_distance_km": r["port_distance_km"],
        })
    counts = Counter(r["status"] for r in out)
    if pilot:
        ref = defaultdict(list)
        for p in read_csv(pilot):
            if p["year"] != "2010" or p["lsbci"] == "":
                continue
            ref[tuple(sorted((p["exporter_iso3"], p["importer_iso3"])))].append(float(p["lsbci"]))
        checks = []
        for key, vals in ref.items():
            yrs = bypair.get(key, {})
            if 2010 in yrs:
                checks.append(abs(sum(vals) / len(vals) - yrs[2010]))
        audit["pilot_matched_country_pairs"] = len(checks)
        audit["pilot_max_absolute_difference"] = max(checks) if checks else None
        if len(checks) < 25 or max(checks, default=0) > 2e-4:
            raise ValueError("Full global source does not reproduce 2010 subset pilot")
    audit.update({
        "status": "SOURCE_VERIFIED_COUNTRY_BILATERAL_COHORT_NOT_PORT_OD",
        "source_url": RAW_URL,
        "source_sha256": RAW_SHA256,
        "pilot_subset_url": PILOT_URL if pilot else None,
        "pilot_blob_sha": PILOT_BLOB if pilot else None,
        "source_port_count": 194,
        "eligible_islands": 51,
        "single_country_islands": len(stable),
        "multicountry_islands": mixed,
        "candidate_pairs": len(out),
        "pair_status_counts": dict(sorted(counts.items())),
        "matched_2010_islands": len({i for r in out
                                     if r["status"] == "INDEX_2010_OBSERVED"
                                     for i in (r["island_a"], r["island_b"])}),
        "matched_2010_country_pairs": len({
            tuple(sorted((r["country_a"], r["country_b"])))
            for r in out if r["status"] == "INDEX_2010_OBSERVED"}),
        "source_git_blobs": PINNED_BLOBS,
        "limits": [
            "The 2010 LSBCI is between countries, not measured voyages between matched WPI ports.",
            "Same-country island pairs have no bilateral score and are NOT zero-ship pairs.",
            "2010 index can precede 2011-2022 butterflies but not all 2001-2010 plants.",
            "The 2001-2024 candidate host catalogue is retrospective and may postdate butterflies.",
            "No shipping/host mediation, direct vessel transport or population establishment is identified.",
            "The 51-island host-eligible cohort is selected on plant and previous butterfly observation.",
        ],
    })
    return out, audit


def ranks(x):
    ordered = sorted(range(len(x)), key=lambda k: x[k])
    result = np.empty(len(x), dtype=float)
    start = 0
    while start < len(x):
        end = start + 1
        while end < len(x) and x[ordered[end]] == x[ordered[start]]:
            end += 1
        result[ordered[start:end]] = (start + 1 + end) / 2
        start = end
    return result


def partial_rank(rows, left, right, controls=()):
    x = ranks([float(r[left]) for r in rows])
    y = ranks([float(r[right]) for r in rows])
    if controls:
        z = np.column_stack(
            [np.ones(len(rows))] + [ranks([float(r[c]) for r in rows]) for c in controls])
        x = x - z @ np.linalg.lstsq(z, x, rcond=None)[0]
        y = y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
    if np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def measure(rows, audit):
    selected = [r for r in rows if r["status"] == "INDEX_2010_OBSERVED"]
    if len(selected) < 30:
        raise ValueError("Insufficient historical bilateral country index support")
    for r in selected:
        r["log_port_distance"] = math.log1p(float(r["port_distance_km"]))
    # 1275 pairs are dependent; do not report conventional pair-independent p.
    comparisons = {
        "ship_host_rho": ("lsbci_2010", "host_jaccard", ()),
        "ship_butterfly_rho": ("lsbci_2010", "new_butterfly_jaccard", ()),
        "host_butterfly_rho": ("host_jaccard", "new_butterfly_jaccard", ()),
        "ship_host_adjust_distance_prior": (
            "lsbci_2010", "host_jaccard", ("log_port_distance", "earlier_butterfly_jaccard")),
        "ship_butterfly_adjust_distance_prior": (
            "lsbci_2010", "new_butterfly_jaccard", ("log_port_distance", "earlier_butterfly_jaccard")),
        "ship_butterfly_adjust_distance_prior_host": (
            "lsbci_2010", "new_butterfly_jaccard",
            ("log_port_distance", "earlier_butterfly_jaccard", "host_jaccard")),
        "host_butterfly_adjust_distance_prior_ship": (
            "host_jaccard", "new_butterfly_jaccard",
            ("log_port_distance", "earlier_butterfly_jaccard", "lsbci_2010")),
    }
    vals = {name: partial_rank(selected, left, right, controls)
            for name, (left, right, controls) in comparisons.items()}
    island_ids = sorted({i for r in selected for i in (r["island_a"], r["island_b"])})
    country_ids = sorted({i for r in selected for i in (r["country_a"], r["country_b"])})
    ranges = {}
    for name in ("ship_host_adjust_distance_prior",
                 "ship_butterfly_adjust_distance_prior_host",
                 "host_butterfly_adjust_distance_prior_ship"):
        l, r, controls = comparisons[name]
        v = [partial_rank([z for z in selected if id_ not in
                          (z["island_a"], z["island_b"])], l, r, controls)
             for id_ in island_ids]
        good = [x for x in v if x is not None]
        ranges[name] = [min(good), max(good)] if good else None
    audit.update({
        "sample_size": {"island_pairs": len(selected),
                        "island_nodes": len(island_ids),
                        "country_nodes": len(country_ids),
                        "country_pairs": audit["matched_2010_country_pairs"]},
        "rank_associations": vals,
        "leave_one_island_ranges": ranges,
        "note": "Post-hoc descriptive correlations; country-pair LSBCI values repeat over multiple island pairs.",
    })
    return audit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-csv", type=Path, required=True)
    ap.add_argument("--pilot-csv", type=Path)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()
    validate_source(args.raw_csv)
    out, audit = aggregate(args.raw_csv, args.pilot_csv)
    audit = measure(out, audit)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(out[0]))
        wr.writeheader()
        wr.writerows(out)
    args.output_json.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n",
                                encoding="utf8")
    print(json.dumps({k: v for k, v in audit.items()
                      if k not in ("source_git_blobs", "unknown_country_labels_top20")},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

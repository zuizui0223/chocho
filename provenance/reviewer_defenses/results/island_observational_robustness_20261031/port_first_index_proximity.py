#!/usr/bin/env python3
"""Screen whether newly GBIF-indexed island butterflies are nearer ports.

This is strictly retrospective detectability evidence, not true arrival.
Requires ORIGINAL per-record coordinates: an island-year aggregate cannot be
substituted. Within the same island and year, other butterfly records are the
minimum observation-location control. No unmatched island is port-free.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

RADIUS_KM = 6371.0088


def read_csv(path):
    with Path(path).open(encoding="utf8", newline="") as f:
        yield from csv.DictReader(f)


def haversine(lat1, lon1, lat2, lon2):
    a, b = math.radians(lat1), math.radians(lat2)
    dlat, dlon = b - a, math.radians(lon2 - lon1)
    h = math.sin(dlat / 2) ** 2 + math.cos(a) * math.cos(b) * math.sin(dlon / 2) ** 2
    return 2 * RADIUS_KM * math.asin(math.sqrt(min(1., max(0., h))))


def build_ports(port_rows):
    ports = defaultdict(list)
    for x in port_rows:
        if not all(k in x for k in ("port_id", "landmass_id", "latitude", "longitude")):
            raise ValueError("port CSV needs port_id,landmass_id,latitude,longitude")
        if not x["port_id"] or not x["landmass_id"]:
            raise ValueError("port ID and landmass ID are mandatory")
        lat, lon = float(x["latitude"]), float(x["longitude"])
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            raise ValueError("invalid port coordinates")
        ports[x["landmass_id"]].append((lat, lon))
    if not ports:
        raise ValueError("empty 2019-WPI matched port mapping")
    return ports


def collect(occurrences, ports):
    """Memory grows with species × island × year, not record count."""
    years = defaultdict(lambda: [0, 0.])
    all_island_year = defaultdict(lambda: [0, 0.])
    species_per_island_year = defaultdict(set)
    count = defaultdict(int)
    for r in occurrences:
        count["seen"] += 1
        necessary = ("speciesKey", "landmass_id", "year",
                     "decimalLatitude", "decimalLongitude")
        if any(k not in r for k in necessary):
            raise ValueError("occurrence CSV requires " + ",".join(necessary))
        island = r["landmass_id"]
        if island not in ports:
            count["no_geocoded_port_support"] += 1
            continue
        try:
            year = int(r["year"])
            lat = float(r["decimalLatitude"])
            lon = float(r["decimalLongitude"])
        except (ValueError, TypeError):
            count["missing_location_or_year"] += 1
            continue
        if not (2001 <= year <= 2022 and -90 <= lat <= 90 and -180 <= lon <= 180):
            count["outside_period_or_invalid_location"] += 1
            continue
        species = r["speciesKey"]
        if not species:
            count["missing_taxon_key"] += 1
            continue
        d = min(haversine(lat, lon, a, b) for a, b in ports[island])
        key = (island, species, year)
        years[key][0] += 1
        years[key][1] += math.log1p(d)
        all_island_year[(island, year)][0] += 1
        all_island_year[(island, year)][1] += math.log1p(d)
        species_per_island_year[(island, year)].add(species)
        count["eligible_records"] += 1
    return years, all_island_year, species_per_island_year, dict(count)


def audit(occurrences, port_rows, min_other_records=10,
          min_earlier_years=2, min_earlier_species=5):
    ports = build_ports(port_rows)
    obs, background, detected, count = collect(occurrences, ports)
    keys = defaultdict(list)
    for island, species, year in obs:
        keys[(island, species)].append(year)
    findings = []
    exclusions = defaultdict(int)
    for (island, species), seen_years in keys.items():
        first_year = min(seen_years)
        if not 2011 <= first_year <= 2022:
            exclusions["not_new_index_2011_2022"] += 1
            continue
        prior = [y for y in range(max(2001, first_year - 5), first_year)
                 if len(detected.get((island, y), set()) - {species}) >= min_earlier_species]
        if len(prior) < min_earlier_years:
            exclusions["insufficient_prior_butterfly_recording_opportunity"] += 1
            continue
        n, value = obs[(island, species, first_year)]
        bg_n, bg_value = background[(island, first_year)]
        other_n, other_value = bg_n - n, bg_value - value
        if other_n < min_other_records:
            exclusions["insufficient_same_island_year_background"] += 1
            continue
        first_mean = value / n
        background_mean = other_value / other_n
        later = []
        for y in sorted(seen_years):
            if first_year < y <= min(2022, first_year + 5):
                q, v = obs[(island, species, y)]
                bn, bv = background[(island, y)]
                if bn - q >= min_other_records:
                    later.append((v / q) - ((bv - v) / (bn - q)))
        findings.append({
            "landmass_id": island, "speciesKey": species,
            "first_index_year": first_year, "first_records": n,
            "prior_opportunity_years": len(prior),
            "same_year_other_butterfly_records": other_n,
            "first_relative_log1p_nearest_port_km": first_mean - background_mean,
            "mean_relative_log1p_nearest_port_km_later5y":
                (sum(later) / len(later)) if later else None,
            "difference_first_minus_later5y": (
                first_mean - background_mean - sum(later) / len(later)
                if later else None),
            "later_eligible_years": len(later),
        })
    island_differences = defaultdict(list)
    island_first = defaultdict(list)
    for row in findings:
        island_first[row["landmass_id"]].append(row["first_relative_log1p_nearest_port_km"])
        if row["difference_first_minus_later5y"] is not None:
            island_differences[row["landmass_id"]].append(
                row["difference_first_minus_later5y"])
    means = [sum(v)/len(v) for v in island_first.values()]
    later_means = [sum(v)/len(v) for v in island_differences.values()]
    return {
        "status": "OBSERVED_FIRST_INDEX_LOCATION_NOT_ARRIVAL",
        "record_audit": count, "screened_species_island_cases": len(keys),
        "eligible_first_year_species_island_cases": len(findings),
        "eligible_islands": len(island_first),
        "island_equal_first_minus_other_mean_log_distance":
            (sum(means)/len(means) if means else None),
        "eligible_first_and_later_islands": len(later_means),
        "island_equal_first_minus_later_relative_log_distance":
            (sum(later_means)/len(later_means) if later_means else None),
        "reasons_unmatched": dict(exclusions),
        "cases": findings,
        "limits": [
            "GBIF earliest indexed year does not establish biological arrival year.",
            "2019-WPI ports are anachronistic proxies for the full 2001-2022 records.",
            "Same-island-year background is a basic detection-location control, not a provider matched experiment.",
            "Nearest port is a geographic distance, not vessel connection, direct maritime route or host transport.",
            "No butterfly recording is never coded as confirmed species absence.",
            "Inference should be assessed by island-level clusters and independent data source matching.",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--occurrence-csv", type=Path, required=True)
    p.add_argument("--island-wpi-ports-csv", type=Path, required=True)
    p.add_argument("--output-json", type=Path, required=True)
    args = p.parse_args()
    result = audit(read_csv(args.occurrence_csv),
                   list(read_csv(args.island_wpi_ports_csv)))
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2)
                                + "\n", encoding="utf8")
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()

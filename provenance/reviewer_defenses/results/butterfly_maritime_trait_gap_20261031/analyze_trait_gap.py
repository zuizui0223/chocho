#!/usr/bin/env python3
"""Reproduce a species-level, non-causal screen of butterfly island record gaps.

All input CSVs are already public source aggregates in the chocho repository.
Ports contain pooled 2001-2024 plant records and MUST NOT be treated as
prior-year larval-host supply or ship routes. Wing size is NOT actual migration.
The outcome is a minimum previous-GBIF-record gap, NOT a flight distance.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
SOURCES = {
    "traits": "data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "species": "provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv",
    "gaps": "provenance/reviewer_defenses/results/original_island_continuity_20261030/global_original_continuity_species_counts.csv",
    "ports": "provenance/reviewer_defenses/results/port_connectivity_hosts_20261030/port_connectivity_host_all225.csv",
}
GIT_BLOB_PINS = {
    "traits": "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
    "species": "39fe84981cab2b0e47063399f0a5c0c9731e800d",
    "gaps": "17bf86ad6bd3622d36cea738a363f3e76774cf04",
    "ports": "239dc82992fa3374e4862a1e2202663e1aeef7da",
}


def read_csv(kind: str) -> list[dict[str, str]]:
    with (ROOT / SOURCES[kind]).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def make_screen() -> tuple[list[dict], list[dict]]:
    traits = {r["species"]: r for r in read_csv("traits")}
    species_rows = read_csv("species")
    ports = {r["taxon_key"]: r for r in read_csv("ports")}
    assert len(species_rows) == 225 and len(ports) == 225
    by_key: dict[str, dict] = {}
    for r in read_csv("gaps"):
        if r["cohort"] != "both_period_opportunity":
            continue
        key = r["speciesKey"]
        n, k = int(r["total_events"]), int(r["events"])
        x = by_key.setdefault(key, {"n": n, "far500": 0, "far100": 0, "sum": 0})
        assert x["n"] == n
        x["sum"] += k
        if r["distance_band"] == "500 km以上":
            x["far500"] += k
            x["far100"] += k
        elif r["distance_band"] == "100–500 km":
            x["far100"] += k
        else:
            assert r["distance_band"] in {"10 km未満", "10–100 km"}
    assert all(x["sum"] == x["n"] for x in by_key.values())

    included, screen = [], []
    for r in species_rows:
        key, name = r["taxon_key"], r["species"]
        t, p, g = traits.get(name), ports.get(key), by_key.get(key)
        status = "MATCHED" if t and p and g and g["n"] > 0 else "NO_BOTH_PERIOD_OPPORTUNITY_OR_TRAIT"
        screen.append({"species": name, "taxon_key": key, "disposition": status,
                       "both_period_index_events": g["n"] if g else 0})
        if status != "MATCHED":
            continue
        x = {
            "species": name,
            "speciesKey": key,
            "n": g["n"],
            "far500": g["far500"],
            "far100": g["far100"],
            "wing": float(t["wing_size_proxy"]),
            "host": float(t["introduced_share_of_contemporary"]),
            "port": float(p["island_ports_with_matching_alien_host"]),
            "host_species": float(t["resolved_host_species"]),
            "native": float(t["native_resource_units"]),
            "old": float(r["earlier_records_on_eligible_islands"]),
            "news": float(r["new_records_on_eligible_islands"]),
        }
        assert all(math.isfinite(x[k]) for k in ("wing", "host", "port", "host_species",
                                                "native", "old", "news"))
        included.append(x)
    assert len(included) == 193 and len(screen) == 225
    assert sum(x["n"] for x in included) == 3136
    assert sum(x["far500"] for x in included) == 164
    included.sort(key=lambda x: int(x["speciesKey"]))
    return included, screen


def rank(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda p: p[1])
    out = [0.0] * len(values)
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][1] == ordered[i][1]:
            j += 1
        av = (i + j + 1) / 2.0
        for k in range(i, j):
            out[ordered[k][0]] = av
        i = j
    return out


def pearson(a: list[float], b: list[float]) -> float | None:
    n = len(a)
    if n < 3:
        return None
    am, bm = sum(a) / n, sum(b) / n
    aa = sum((x-am)**2 for x in a)
    bb = sum((x-bm)**2 for x in b)
    ab = sum((a[i]-am)*(b[i]-bm) for i in range(n))
    return ab / math.sqrt(aa * bb) if aa > 0 and bb > 0 else None


def spearman(a: list[float], b: list[float]) -> float | None:
    return pearson(rank(a), rank(b))


def solve(matrix: list[list[float]], vector: list[float]) -> list[float] | None:
    n = len(vector)
    aug = [matrix[i][:] + [vector[i]] for i in range(n)]
    for k in range(n):
        best = max(range(k, n), key=lambda j: abs(aug[j][k]))
        if abs(aug[best][k]) < 1e-12:
            return None
        aug[k], aug[best] = aug[best], aug[k]
        v = aug[k][k]
        for j in range(k, n + 1):
            aug[k][j] /= v
        for i in range(n):
            if i == k:
                continue
            a = aug[i][k]
            for j in range(k, n + 1):
                aug[i][j] -= a * aug[k][j]
    return [row[-1] for row in aug]


def residual(y: list[float], controls: list[list[float]]) -> list[float] | None:
    xx = [[1.0] + [col[i] for col in controls] for i in range(len(y))]
    k = len(xx[0])
    aa = [[sum(x[i]*x[j] for x in xx) for j in range(k)] for i in range(k)]
    bb = [sum(x[i]*z for x, z in zip(xx, y)) for i in range(k)]
    coefs = solve(aa, bb)
    if coefs is None:
        return None
    return [v - sum(xx[i][j]*coefs[j] for j in range(k))
            for i, v in enumerate(y)]


def partial_rank(x: list[float], y: list[float],
                 controls: list[list[float]]) -> float | None:
    a = residual(rank(x), [rank(c) for c in controls])
    b = residual(rank(y), [rank(c) for c in controls])
    return pearson(a, b) if a is not None and b is not None else None


def statistics(rows: list[dict]) -> dict[str, float | None]:
    wing = [x["wing"] for x in rows]
    host = [x["host"] for x in rows]
    port = [x["port"] for x in rows]
    y = [x["far500"]/x["n"] for x in rows]
    y100 = [x["far100"]/x["n"] for x in rows]
    l_n = [math.log1p(x["n"]) for x in rows]
    l_h = [math.log1p(x["host_species"]) for x in rows]
    l_nat = [math.log1p(x["native"]) for x in rows]
    return {
        "wing_gap": spearman(wing, y),
        "host_gap": spearman(host, y),
        "ports_gap": spearman(port, y),
        "wing_100": spearman(wing, y100),
        "ports_prior": spearman(port, [x["old"] for x in rows]),
        "ports_new": spearman(port, [x["news"] for x in rows]),
        "partial_wing": partial_rank(wing, y, [l_n, l_h, l_nat, host]),
        "partial_host": partial_rank(host, y, [l_n, l_h, l_nat, wing]),
        "partial_ports": partial_rank(port, y, [l_n, l_h, l_nat, wing, host]),
    }


def percentile(v: list[float], p: float) -> float:
    vals = sorted(x for x in v if x is not None and math.isfinite(x))
    i = (len(vals)-1)*p
    lo, hi = math.floor(i), math.ceil(i)
    return vals[lo] + (vals[hi] - vals[lo]) * (i - lo)


def report(rows: list[dict], n_boot: int = 800) -> dict:
    points = statistics(rows)
    boot = {k: [] for k in points}
    state = 20261031
    for _ in range(n_boot):
        sample = []
        for _j in range(len(rows)):
            state = (state * 1664525 + 1013904223) & 0xffffffff
            sample.append(rows[state % len(rows)])
        for k, v in statistics(sample).items():
            if v is not None and math.isfinite(v):
                boot[k].append(v)
    return {
        "schema": "chocho_butterfly_maritime_trait_record_gap_v1",
        "status": "EXPLORATORY_NO_RELIABLE_CAUSAL_HOST_OR_PORT_RESULT",
        "source_paths": SOURCES,
        "source_git_blob_pins": GIT_BLOB_PINS,
        "cohort": "both_period_opportunity",
        "n_candidate_species": 225,
        "n_matched_species": len(rows),
        "n_index_events": sum(x["n"] for x in rows),
        "n_record_gaps_500km": sum(x["far500"] for x in rows),
        "n_species_with_gap_500km": sum(x["far500"] > 0 for x in rows),
        "n_species_with_port_candidate_host_index": sum(x["port"] > 0 for x in rows),
        "species_resample_seed": 20261031,
        "bootstrap_draws": n_boot,
        "correlations": {
            k: {"point": v, "species_resampled_95": [
                percentile(boot[k], .025), percentile(boot[k], .975)
            ] if boot[k] else None,
                "bootstrap_n": len(boot[k])}
            for k, v in points.items()
        },
        "minimum_event_sensitivity": {
            str(n): {"species": len(s), "correlations": statistics(s)}
            for n in (5, 10, 20)
            for s in [[x for x in rows if x["n"] >= n]]
        },
        "limits": [
            "Distances link GBIF indices, not observed butterfly flights.",
            "2001-2024 port flora can postdate butterfly records and is NOT a historical exposure.",
            "Wing size is a limited dispersal proxy, NOT observed migration.",
            "Global introduced-host share is not island-specific or independently dated.",
            "Species bootstrap fixes observation/host/geometry errors and ignores shared islands.",
            "Historical and new-island index counts are reporting-sensitive, NOT negative causal controls.",
            "Missing plant records cannot be coded as biological absence.",
            "No within-island port versus non-port plant control in these four source CSVs.",
        ],
    }


def main() -> None:
    rows, screened = make_screen()
    result = report(rows)
    with (OUT / "matched_species_screen.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "all225_disposition.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(screened[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(screened)
    (OUT / "trait_gap_results.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: result[k] for k in
                      ("n_candidate_species", "n_matched_species", "n_index_events",
                       "n_record_gaps_500km")}, indent=2))


if __name__ == "__main__":
    main()

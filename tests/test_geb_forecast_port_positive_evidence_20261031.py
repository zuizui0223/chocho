"""Guard positive forecast, island and port correlations in GEB gateway synthesis.

Numbers here are source-derived and descriptive. This file does not establish
port-caused dispersal, plant delivery by ships, larval use or colonisation.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "provenance/reviewer_defenses/results"
SYNTHESIS = R / "geb_port_island_realization_20261031/GEB_port_island_realization_results_ja.md"
LEDGER = ROOT / "provenance/reviewer_defenses/EVIDENCE_STATUS_20261031.md"
FORECAST = R / "host_forecast_validation_20261029"
ISLANDS = R / "island_host_network_20261030/island_host_network_pairs.csv"
SHARED = R / "island_shared_host_match_20261030/island_shared_host_match_audit.json"
NULL = R / "island_host_null_20261030/island_host_family_null_audit.json"
PORT = R / "geb_port_island_realization_20261031/GEB_to_port_species_bridge.json"


def blob(path: Path) -> str:
    raw = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\x00" + raw).hexdigest()


def ranks(seq):
    order = sorted(range(len(seq)), key=lambda i: seq[i])
    result = [0.0] * len(seq)
    k = 0
    while k < len(seq):
        j = k + 1
        while j < len(seq) and seq[order[k]] == seq[order[j]]:
            j += 1
        for ind in order[k:j]:
            result[ind] = (k + j + 1) / 2
        k = j
    return result


def pearson(a, b):
    n = len(a)
    av = sum(a) / n
    bv = sum(b) / n
    aa = sum((x - av) ** 2 for x in a)
    bb = sum((x - bv) ** 2 for x in b)
    return sum((x - av) * (y - bv) for x, y in zip(a, b)) / math.sqrt(aa * bb)


def solve(m, v):
    n = len(v)
    a = [[*m[i], v[i]] for i in range(n)]
    for i in range(n):
        pivot = max(range(i, n), key=lambda k: abs(a[k][i]))
        a[i], a[pivot] = a[pivot], a[i]
        den = a[i][i]
        if abs(den) < 1e-10:
            raise AssertionError("A source changed and made the rank covariates singular")
        for k in range(i, n + 1):
            a[i][k] /= den
        for j in range(n):
            if j == i:
                continue
            z = a[j][i]
            for k in range(i, n + 1):
                a[j][k] -= z * a[i][k]
    return [r[n] for r in a]


def residuals(y, z):
    n = len(y)
    p = len(z[0])
    gram = [[sum(z[i][j] * z[i][k] for i in range(n)) for k in range(p)]
            for j in range(p)]
    rhs = [sum(z[i][j] * y[i] for i in range(n)) for j in range(p)]
    beta = solve(gram, rhs)
    return [y[i] - sum(z[i][j] * beta[j] for j in range(p))
            for i in range(n)]


def test_forecast_improvement_verified_for_three_windows():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island-port evidence not included in independent blinded GEB paper")
    sources = [
        ("forecast_temporal_results.json", "77c3663b8d6b0e1dab732becb8a2a4f03706c857",
         208482, 0.041, "4.12%"),
        ("forecast_spatial_results.json", "ee541ca8208def59aa13abbd4d1a3cb36bbb9fb0",
         208482, 0.030, "3.00%"),
        ("forecast_spatial_buffer100_results.json", "f432b78c628fc46a198b51d6ddf5bf88bf605010",
         102728, 0.0274, "2.74%"),
    ]
    doc = SYNTHESIS.read_text(encoding="utf8")
    ledger = LEDGER.read_text(encoding="utf8")
    for file, sha, n, expected, label in sources:
        path = FORECAST / file
        assert blob(path) == sha
        d = json.loads(path.read_text(encoding="utf8"))
        m3 = d["results"]["M3_both"]["species_equal_mse"]
        m5 = d["results"]["M5_host_composition"]["species_equal_mse"]
        assert d["results"]["M3_both"]["rows"] == n
        assert d["results"]["M5_host_composition"]["rows"] == n
        assert m5 < m3
        assert abs((m3 - m5) / m3 - expected) < 0.0005
        assert label in doc
    assert all(s in doc for s in ("155", "208,482", "119/155",
                                  "2.98%", "3.76%", "5.55%", "5.34%",
                                  "未学習の蝶種への転用は支持されない"))
    assert "119/155種" in ledger
    assert "2.74%" in ledger


def test_island_similarity_and_plant_network_overlap_recalculated():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island botanical source deliberately excluded from GEB review")
    assert blob(ISLANDS) == "5663d2ffe42419f5c4496b92c9f3bbf76175907b"
    with ISLANDS.open(newline="", encoding="utf8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1275
    h = ranks([float(r["host_jaccard"]) for r in rows])
    newer = ranks([float(r["new_butterfly_jaccard"]) for r in rows])
    prior = ranks([float(r["earlier_butterfly_jaccard"]) for r in rows])
    distance = ranks([float(r["port_distance_km"]) for r in rows])
    assert math.isclose(pearson(h, newer), 0.29980576483377025, abs_tol=1e-11)
    assert math.isclose(pearson(h, prior), 0.4178255833828515, abs_tol=1e-11)
    z = [[1.0, distance[i], prior[i]] for i in range(len(rows))]
    assert math.isclose(pearson(residuals(h, z), residuals(newer, z)),
                        0.2153535842700985, abs_tol=1e-11)
    original = json.loads(SHARED.read_text(encoding="utf8"))
    reference = json.loads(NULL.read_text(encoding="utf8"))
    assert blob(SHARED) == "230b5601a6bc292eeb7996d2d3aebdc2ebbfaf5d"
    assert blob(NULL) == "91402de9b98b17838b8e1e9d29fde310d5bc2bb4"
    assert original["shared_new_species_pairs"] == reference["both_new_species_pairs"] == 1163
    assert original["matched_shared_new_species_pairs"] == 208
    assert reference["draws"] == 199
    assert math.isclose(original["fraction_matched"], 208/1163)
    assert math.isclose(reference["null_mean_fraction"], 0.11896109956489237)
    doc = SYNTHESIS.read_text(encoding="utf8")
    assert all(v in doc for v in ("ρ = +0.300", "ρ = +0.418",
                                 "+0.215", "17.9%", "11.9%",
                                 "5/1,163", "1,275", "51島"))


def test_gEB_to_port_record_correlation_not_transport_causality():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory port source omitted from independent GEB manuscript")
    assert blob(PORT) == "e2d942dbba4844096eb4fe11677ff471d3e0237d"
    port = json.loads(PORT.read_text(encoding="utf8"))
    cohort = port["aligned_port_candidate_comparison"]
    assert cohort["species"] == 225
    assert math.isclose(
        cohort["spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"],
        0.6524282493600642, abs_tol=1e-12)
    assert math.isclose(
        cohort["partial_rank_assoc_after_adjusting_log1p_host_count_native_units_and_prior_recorded_islands"],
        0.5495755354184105, abs_tol=1e-12)
    assert "2001–2024" in cohort["temporal_leak"]
    doc = SYNTHESIS.read_text(encoding="utf8")
    ledger = LEDGER.read_text(encoding="utf8")
    assert all(s in doc for s in ("+0.652", "+0.550", "39港でρ=+0.023",
                                  "予測改善4.1%・島間ρ0.300・港種別ρ0.652"))
    assert "港を通じて実現した割合" in ledger
    assert "島間" in ledger and "+0.652" not in ledger or "0.652" in ledger

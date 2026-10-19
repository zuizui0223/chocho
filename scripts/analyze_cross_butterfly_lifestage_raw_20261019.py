#!/usr/bin/env python3
"""Original-data, stage-specific butterfly host-performance audit (2026-10-19).

Exploratory, independent reanalysis of CC0 primary data already published by
Ferreira & Rodrigues, Saldivar & Wilson-Rankin, and Verspagen et al.
NEVER interpret these studies as 239-species global recruitment measurements.
Python standard library only. No network, no synthetic butterfly outcomes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "exploratory" / "realized_host_20261019"
FILES = {
    "vanessa/larv.plants.csv": ("0af2cbb249494e11e8dc81a68dcf7bae", 69769, 1154),
    "vanessa/larv.switched.csv": ("af98d718442e42d56a89d728478ff21e", 33476, 530),
    "danaus/field_performance.txt": ("5437505ae948d0133927e44d18acc617", 2695, 54),
    "danaus/lab_performance.txt": ("4e0c3fd3866a9668bab975398a01d708", 1852, 40),
    "melitaea/Survival_data.csv": ("373acf78433faddd218bc9d9c427fc9a", 48291, 1680),
}
DOI = {
    "vanessa": "10.5061/dryad.wdbrv15v5",
    "danaus": "10.5061/dryad.3bk3j9kk2",
    "melitaea": "10.5061/dryad.kwh70rz2m",
}


def exact_bytes(name: str) -> str:
    expected_md5, expected_bytes, _ = FILES[name]
    b = (BASE / name).read_bytes()
    assert len(b) == expected_bytes, f"{name}: byte size {len(b)} != {expected_bytes}"
    actual = hashlib.md5(b).hexdigest()
    assert actual == expected_md5, f"{name}: MD5 {actual} != source {expected_md5}"
    return b.decode("utf-8-sig")


def csv_records(name: str) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(exact_bytes(name), newline="")))


def _txt_rows(name: str, maxsplit: int, required_header: str):
    lines = exact_bytes(name).splitlines()
    assert lines[1].strip() == required_header
    data = [line.split(maxsplit=maxsplit) for line in lines[2:] if line.strip()]
    assert all(len(row) == maxsplit + 1 for row in data), name
    return data


def vanessa():
    outputs = {}
    combined = []
    for name in ("vanessa/larv.plants.csv", "vanessa/larv.switched.csv"):
        rows = csv_records(name)
        assert len(rows) == FILES[name][2]
        assert set(x["status"] for x in rows) == {"2"}, "censor status must not be silently recoded"
        assert all(float(x["time"]) >= 0 for x in rows)
        # Source publication reports FOUR reached pupation; the tables do not
        # contain a success identifier and all 1,684 rows are coded 'dead'.
        per_treatment = {}
        for sp in sorted(set(x["species"] for x in rows)):
            xs = [x for x in rows if x["species"] == sp]
            age = [float(x["time"]) for x in xs]
            per_treatment[sp] = {
                "n": len(xs),
                "event_code_2": len(xs),
                "recorded_time_mean": round(statistics.mean(age), 5),
                "recorded_time_median": round(statistics.median(age), 5),
                "recorded_time_max": max(age),
                "provenance_counts": dict(sorted(Counter(x["source"] for x in xs).items())),
                "distinct_species_trial_groups": len(set(x["plant.trial"] for x in xs)),
            }
        outputs[name.split("/")[-1]] = {
            "n": len(rows), "statuses": dict(Counter(x["status"] for x in rows)),
            "species": per_treatment,
            "suspect_avg_temp_ge_40C": sum(float(x["avg.temp"]) >= 40 for x in rows),
        }
        combined.extend(rows)
    assert len(combined) == 1684
    return {
        "doi": DOI["vanessa"],
        "observations": len(combined),
        "all_status_2": True,
        "pupation_endpoint_identifiable_from_CSV": False,
        "source_paper_reports_4_reaching_pupation": True,
        "source_paper_reports_2_completed_chrysalis": True,
        "published_total_includes_all_1684": True,
        "files": outputs,
        "interpretation": "Only death-event time (as coded) is available; all status=2 despite published pupation; no rederived pupation success, origin effect, or adult production.",
    }


def danaus():
    field = _txt_rows(
        "danaus/field_performance.txt", 5,
        "obs id host.plant treatment survival info")
    lab = _txt_rows(
        "danaus/lab_performance.txt", 7,
        "obs id host.plant danaini.sp survival dev.time wing.len info")
    assert len(field) == 54 and len(lab) == 40
    fd = {}
    for host in ("Asclepias", "Calotropis"):
        for exposure in ("bagged", "non.bagged"):
            subset = [x for x in field if x[2] == host and x[3] == exposure]
            y = Counter(x[4] for x in subset)
            assert set(y).issubset({"0", "1", "NA"})
            n = y["0"] + y["1"]
            fd[f"{host}|{exposure}"] = {
                "n_assigned": len(subset), "n_observed": n,
                "positive_stage1_2": y["1"], "negative": y["0"],
                "unobserved_cattle_or_other": y["NA"],
                "positive_rate_observed": (round(y["1"] / n, 8) if n else None),
            }
    assert sum(z["n_assigned"] for z in fd.values()) == 54
    assert sum(z["unobserved_cattle_or_other"] for z in fd.values()) == 4
    for host in ("Asclepias", "Calotropis"):
        fd[f"{host}|bagging_minus_exposed_pp"] = round(
            100 * (fd[f"{host}|bagged"]["positive_rate_observed"] -
                   fd[f"{host}|non.bagged"]["positive_rate_observed"]), 5
        )
    ld = {}
    for species in ("D.erippus", "D.gilippus"):
        for host in ("Asclepias", "Calotropis"):
            subset = [x for x in lab if x[3] == species and x[2] == host]
            successes = [x for x in subset if x[4] == "1"]
            assert len(subset) == 10
            assert set(x[4] for x in subset).issubset({"0", "1"})
            assert len(successes) == sum(x[6] != "NA" for x in subset)
            ld[f"{species}|{host}"] = {
                "n_assigned": len(subset), "survived": len(successes),
                "survival_rate": round(len(successes) / 10, 8),
                "development_days_survivors_mean": round(
                    statistics.mean(float(x[5]) for x in successes), 8),
                "development_days_survivors_n": len(successes),
            }
        ld[f"{species}|exotic_minus_native_days_survivors"] = round(
            ld[f"{species}|Calotropis"]["development_days_survivors_mean"] -
            ld[f"{species}|Asclepias"]["development_days_survivors_mean"], 8
        )
    return {
        "doi": DOI["danaus"], "lab": ld, "field": fd,
        "lab_n": 40, "field_n": 54, "field_missing_n": 4,
        "field_site_confounding": "Native Asclepias and exotic Calotropis were at different field sites; between-host field contrast is not causal.",
        "field_stage": "Survival to first/second instar status, NOT later adult emergence.",
        "lab_stage": "Source survival (adult wing measurement on survivors); development time calculated among survivors only.",
        "no_linked_egg_allocation": True,
    }


def percentile_sorted(s: list[float], q: float) -> float:
    return s[int((len(s) - 1) * q)]


def family_ci(values, nboot=9999, seed=20261019):
    rng = random.Random(seed)
    n = len(values)
    mu = statistics.mean(values)
    boots = []
    for _ in range(nboot):
        boots.append(sum(values[rng.randrange(n)] for _ in range(n)) / n)
    boots.sort()
    return {
        "estimate": round(mu, 8),
        "cluster_bootstrap_95ci": [
            round(percentile_sorted(boots, 0.025), 8),
            round(percentile_sorted(boots, 0.975), 8)
        ],
        "families_positive": sum(v > 0 for v in values),
        "families_negative": sum(v < 0 for v in values),
        "families_zero": sum(v == 0 for v in values),
    }


def melitaea():
    rows = csv_records("melitaea/Survival_data.csv")
    assert len(rows) == 1680
    assert set(x["Dead"] for x in rows) == {"0", "3"}
    assert all(x["Group"] == f'{x["Family"]}_{x["Temperature"]}_{x["Plant"]}' for x in rows)
    temps = ("28", "30", "32", "34")
    plants = ("Pl", "Vs")
    families = sorted(set(x["Family"] for x in rows), key=int)
    assert len(families) == 14
    groups = defaultdict(list)
    for row in rows:
        groups[(row["Family"], row["Temperature"], row["Plant"])].append(row)
    assert len(groups) == 14 * 4 * 2
    assert all(len(x) == 15 for x in groups.values())
    anomalous = [
        {k: x[k] for k in ("Clutch.ID", "Family", "Temperature", "Plant", "Time", "Time2", "Dead")}
        for x in rows if (x["Dead"] == "0" and x["Time2"] != "NA")
        or (x["Dead"] == "3" and (x["Time2"] == "NA" or float(x["Time2"]) < float(x["Time"])))
    ]
    assert len(anomalous) == 4
    def cell_survival(f, temp, plant):
        return sum(x["Dead"] == "0" for x in groups[f, temp, plant]) / 15
    records = {}
    mother_diff = {t: [] for t in temps}
    family_rows = []
    for fam in families:
        values = {}
        for t in temps:
            p = cell_survival(fam, t, "Pl")
            v = cell_survival(fam, t, "Vs")
            values[t] = v - p
            mother_diff[t].append(values[t])
        family_rows.append({
            "family_id": fam,
            "Vs_minus_Pl_pp_by_temperature": {t: round(100 * values[t], 5) for t in temps},
            "mean_Vs_minus_Pl_pp": round(25 * sum(values.values()), 5),
        })
    for t in temps:
        r = {plant: sum(x["Dead"] == "0" for x in rows
                         if x["Temperature"] == t and x["Plant"] == plant)
             for plant in plants}
        records[t] = {
            "n_per_host": 210,
            "reached_diapause_Pl": r["Pl"], "reached_diapause_Vs": r["Vs"],
            "Pl_pct": round(100 * r["Pl"] / 210, 5),
            "Vs_pct": round(100 * r["Vs"] / 210, 5),
            "paired_family_Vs_minus_Pl_pp": round(100 * statistics.mean(mother_diff[t]), 5),
            "mother_cluster_ci": family_ci(mother_diff[t], seed=20261019 + int(t)),
        }
    mean_diff = [sum(mother_diff[t][i] for t in temps) / 4 for i in range(14)]
    interaction = [mother_diff["34"][i] - mother_diff["28"][i] for i in range(14)]
    mean_ci = family_ci(mean_diff)
    inter_ci = family_ci(interaction, seed=20261020)
    leave_one_family_out = [statistics.mean([x for j,x in enumerate(mean_diff) if j!=i]) for i in range(14)]
    return {
        "doi": DOI["melitaea"], "n": len(rows),
        "family_n": len(families),
        "clutch_n": len(set(x["Clutch.ID"] for x in rows)),
        "independent_treatment_group_n": len(groups),
        "each_family_temperature_host_group_n": 15,
        "endpoint": "Dead=0 alive at diapause; Dead=3 died before diapause; NOT adult emergence.",
        "host_codes": {"Pl": "Plantago lanceolata", "Vs": "Veronica spicata"},
        "per_temperature": records,
        "overall_paired_family_Vs_minus_Pl": mean_ci,
        "interaction_Vs_minus_Pl_34_vs_28": inter_ci,
        "leave_one_family_out_overall_range": [
            round(min(leave_one_family_out), 8), round(max(leave_one_family_out), 8)
        ],
        "families_with_positive_and_negative_host_differences_across_temperatures":
            sum(any(mother_diff[t][i] > 0 for t in temps) and
                any(mother_diff[t][i] < 0 for t in temps) for i in range(14)),
        "discordant_Dead_Time2_rows": anomalous,
        "four_discordant_rows_adjudicated": False,
        "data_quality_limit": "Four Dead=0 rows carry non-NA Time2=4; primary uses stated Dead event code. Carry as flagged sensitivity; do not silently relabel.",
        "family_descriptors": family_rows,
    }


def build():
    return {
        "schema": "chocho_20261019_multispecies_original_lifestage_receipt_v1",
        "status": "ORIGINAL_SOURCE_VERIFIED_WITH_ENDPOINT_LIMITS",
        "distinct_butterfly_species": 4,
        "original_larval_records_total": 3458,
        "source_checksums_md5_verified": {k: md5 for k, (md5, _, _) in FILES.items()},
        "source_sizes_bytes": {k: size for k, (_, size, _) in FILES.items()},
        "vanessa": vanessa(), "danaus": danaus(), "melitaea": melitaea(),
        "macro_geography_crosswalk": {
            "Danaus gilippus": {"known_hosts":17, "native_cells":154,"contemporary_cells":170,"introduced_added_cells":16},
            "Vanessa cardui": {"known_hosts":104, "native_cells":282,"contemporary_cells":360,"introduced_added_cells":78},
            "source": "Frozen chocho data/frozen/figure_sources/anthropogenic_species_metrics.csv",
            "not_proof_of_local_host_use": True,
        },
        "source_independence": {
            "n_independent_public_studies":3,
            "n_distinct_butterflies":4,
            "incompatible_lifestage_endpoints": True,
            "total_recruitment_sign_reversal_identified": False,
            "cross_species_effect_meta_analysis_permitted": False,
            "novel_field_observations":0,
            "manuscript_submitted_or_modified":False
        },
        "biological_boundary": [
            "Reported host and temperature associations in primary papers are PRIOR ART; this script is an independent source audit/recomputation, not a new ecological phenomenon.",
            "Melitaea family-cluster bootstrap reflects between-family resampling in a single experiment; it does not make temperature placements, host identity, or sites randomized and does not generalize to global geography.",
            "The original Melitaea host×temperature literature effect may involve developmental traits; this survival-only extreme-temperature interaction is an exploratory diagnostic, not a refutation of published models.",
            "Danaus field experiment assigns different host taxa to two different sites and its endpoint is early instar survival; no direct field exotic causal test.",
            "Vanessa CSVs contain status=2 for all 1684 even though source authors report four reached pupation; no recoverable adult recruitment endpoint.",
            "No matched host-availability, within-cohort egg-diversion, same-cohort adult recruitment or source-sink movement information allows testing nonmonotone total butterfly production."
        ]
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", help="write JSON receipt at path")
    args = p.parse_args()
    ans = build()
    rendered = json.dumps(ans, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(rendered, encoding="utf-8")
    print(json.dumps({
        "status": ans["status"],
        "raw_rows": ans["original_larval_records_total"],
        "vanessa_event_code_only": ans["vanessa"]["all_status_2"],
        "danaus_lab": ans["danaus"]["lab"],
        "melitaea_overall": ans["melitaea"]["overall_paired_family_Vs_minus_Pl"],
        "melitaea_34_vs_28_interaction": ans["melitaea"]["interaction_Vs_minus_Pl_34_vs_28"],
        "melitaea_discordant_rows": len(ans["melitaea"]["discordant_Dead_Time2_rows"]),
        "production_sign_reversal_identified": False
    }, sort_keys=True))

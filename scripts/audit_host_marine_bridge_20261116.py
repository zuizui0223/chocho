#!/usr/bin/env python3
"""Recheck whether hosts and marine predictors supply complementary or
synergistic retrospective information for butterfly-community similarity.

Source: frozen 51-held-island model-fold MSEs, not species colonization data.
No re-fitting, re-querying, or future-first-arrival inference is performed.
Island bootstrap intervals describe resampling of correlated saved folds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

SOURCE = "results/island51_host_circuit_joint_20261106/model_fold_losses.csv"
MODEL_SUMMARY = "results/island51_host_circuit_joint_20261106/model_summary.csv"
COMPARISON_SUMMARY = "results/island51_host_circuit_joint_20261106/comparison_summary.csv"
EXPECTED_SOURCE_GIT_BLOB = "d2df6457895091a3652c6cdad3fd081238dfdcf4"
MODELS = ("G", "H", "C", "HC", "I")
TARGETS = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
CONTRASTS = (("G", "H"), ("G", "C"), ("C", "HC"), ("H", "HC"), ("HC", "I"))
DRAWS = 9999
SEED = 20261116


class LCG:
    """Unsigned 32-bit generator shared with the initial independent JS audit."""
    def __init__(self, seed: int) -> None:
        self.state = seed

    def index(self, n: int) -> int:
        self.state = (self.state * 1664525 + 1013904223) & 0xFFFFFFFF
        return (self.state * n) // 4294967296


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def quantile_type7(values: list[float], p: float) -> float:
    values = sorted(values)
    z = (len(values) - 1) * p
    low = int(z)
    weight = z - low
    return values[low] + (values[min(low + 1, len(values) - 1)] - values[low]) * weight


def summarize(means: dict[str, float]) -> dict[str, float]:
    out = {f"{a}->{b}": 100 * (means[a] - means[b]) / means[a]
           for a, b in CONTRASTS}
    # Positive nonadditivity means the *predicted MSE* benefit of combining
    # main-effect predictors falls short of summing their separate gains.
    # This is NOT a statistical coefficient for biological interaction.
    out["nonadditive_MSE_fraction_G_pct"] = (
        100 * (means["HC"] - means["H"] - means["C"] + means["G"]) / means["G"]
    )
    out["interaction_absolute_MSE"] = means["HC"] - means["I"]
    return out


def aggregate(by_island: dict[int, dict[str, float]], selected: list[int]) -> dict[str, float]:
    return {model: sum(by_island[k][model] for k in selected) / len(selected)
            for model in MODELS}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path,
                        default=Path("results/host_marine_bridge_audit_20261116/bridge_bootstrap.json"))
    args = parser.parse_args()

    source = args.repo_root / SOURCE
    raw = source.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if blob != EXPECTED_SOURCE_GIT_BLOB:
        raise ValueError(f"Original Git blob changed: {blob}")

    data: dict[str, dict[int, dict[str, float]]] = defaultdict(dict)
    for row in read_csv(source):
        target, island, model = row["target"], int(row["heldout_island_id"]), row["model"]
        if target not in TARGETS or model not in MODELS:
            raise ValueError(f"Unexpected source key: {target}/{model}")
        if int(row["test_pairs"]) != 50 or int(row["train_pairs"]) != 1225:
            raise ValueError("Changed held-island fold design")
        record = data[target].setdefault(island, {})
        if model in record:
            raise ValueError("Duplicate island/model fold")
        record[model] = float(row["MSE"])

    expected = {(r["target"], r["model"]): float(r["MSE"])
                for r in read_csv(args.repo_root / MODEL_SUMMARY)}
    if len(expected) != 10 or set(data) != set(TARGETS):
        raise ValueError("Unexpected model summary or target set")
    for target in TARGETS:
        if len(data[target]) != 51 or any(set(x) != set(MODELS) for x in data[target].values()):
            raise ValueError("Each response needs the same complete 51-island, five-model panel")
        av = aggregate(data[target], sorted(data[target]))
        for model in MODELS:
            if abs(av[model] - expected[target, model]) > 1e-12:
                raise ValueError("Saved fold mean disagrees with fixed model summary")

    comparisons = {(r["target"], r["reference_model"], r["augmented_model"]):
                   float(r["relative_MSE_reduction_percent"])
                   for r in read_csv(args.repo_root / COMPARISON_SUMMARY)}
    if len(comparisons) != 10:
        raise ValueError("Missing original prespecified comparisons")

    rng = LCG(SEED)
    reports = {}
    for target in TARGETS:
        tab = data[target]
        ids = sorted(tab)
        n = len(ids)
        point_means = aggregate(tab, ids)
        point = summarize(point_means)
        for a, b in CONTRASTS:
            if abs(point[f"{a}->{b}"] - comparisons[target, a, b]) > 1e-9:
                raise ValueError("Recomputed comparison disagrees with original frozen result")

        bootstrap = {k: [] for k in point}
        for _ in range(DRAWS):
            chosen = [ids[rng.index(n)] for _ in range(n)]
            sample = summarize(aggregate(tab, chosen))
            for key, value in sample.items():
                bootstrap[key].append(value)

        contrasts = {}
        for key, value in point.items():
            leave_one_out = [
                summarize(aggregate(tab, [i for i in ids if i != held]))[key]
                for held in ids
            ]
            contrasts[key] = {
                "point": value,
                "island_fold_bootstrap_95pct_descriptive": [
                    quantile_type7(bootstrap[key], .025),
                    quantile_type7(bootstrap[key], .975),
                ],
                "leave_one_heldout_island_range": [
                    min(leave_one_out), max(leave_one_out)
                ],
                "bootstrap_fraction_positive": sum(x > 0 for x in bootstrap[key]) / DRAWS,
            }
        reports[target] = {
            "held_islands": n,
            "models": point_means,
            "contrasts": contrasts,
            "islands_with_lower_fold_mse": {
                f"{a}->{b}": sum(tab[i][b] < tab[i][a] for i in ids)
                for a, b in CONTRASTS
            },
        }

    result = {
        "analysis": "host_marine_bridge_audit_20261116",
        "source_path": SOURCE,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "source_git_blob": blob,
        "draws": DRAWS,
        "generator": "32-bit LCG, seed 20261116, sorted numeric island IDs, type-7 percentiles",
        "estimand": "Paired retrospective MSE reduction across saved 51 held-island folds",
        "models": {
            "G": "uniform-ocean geographic reference",
            "H": "G plus pooled candidate introduced-host similarity",
            "C": "G plus later-period AIS and wind-derived marine fields",
            "HC": "G plus both host and marine-field main effects",
            "I": "HC plus two host-by-signed-connectivity terms",
        },
        "results": reports,
        "limits": [
            "Bootstrap resamples fixed model folds, not independently re-trained models",
            "The 51 held islands share overlapping pair responses and training sets",
            "Plant host indices pooled 2001–2024 may postdate target butterfly records",
            "Native AIS and wind fields are not measured historical transport",
            "The newer-butterfly response means first GBIF index in selected windows, not actual island colonization",
            "Positive conditional predictive value is not ecological synergy or a causal mechanism",
            "Model I is an exploratory two-interaction predictive extension",
        ],
    }
    output = args.output if args.output.is_absolute() else args.repo_root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({
        "status": "SOURCE_VERIFIED",
        "source_git_blob": blob,
        "output": str(output),
        "later_MSE_reductions_pct": {
            key: reports[TARGETS[0]]["contrasts"][key]["point"]
            for key in ("C->HC", "H->HC", "HC->I")
        },
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()

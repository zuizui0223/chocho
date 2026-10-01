#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save(fig, outdir, stem):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def occurrence_figure(rows, null_result, outdir):
    informative = [
        row for row in rows if int(row["observed_outside_native"]) > 0
    ]
    informative.sort(
        key=lambda row: (
            float(row["fraction_outside_native_explained_by_introduced"] or 0),
            int(row["observed_outside_native"]),
        )
    )
    names = [row["species"] for row in informative]
    fractions = np.asarray([
        float(row["fraction_outside_native_explained_by_introduced"])
        for row in informative
    ])
    totals = np.asarray([int(row["observed_outside_native"]) for row in informative])
    rescued = np.asarray([
        int(row["outside_native_units_explained_by_introduced_host_ranges"])
        for row in informative
    ])

    fig, ax = plt.subplots(figsize=(8.2, 7.8))
    y = np.arange(len(names))
    ax.barh(y, fractions)
    ax.set_yticks(y, labels=names, fontsize=8)
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Fraction of outside-native occurrence units recovered")
    ax.set_title("Introduced host geography recovers independent butterfly occurrences")
    ax.axvline(
        null_result["observed_fraction_recovered"],
        linestyle="--",
        linewidth=1,
        label="Aggregate observed recovery = 57.4%",
    )
    for yi, frac, k, n in zip(y, fractions, rescued, totals):
        ax.text(min(frac + 0.02, 0.98), yi, f"{k}/{n}", va="center", fontsize=7)
    ax.text(
        0.02,
        0.02,
        "Aggregate recovered units: 66/115\n"
        "Regional-composition null: median 49 (95% 43–56), p=5×10⁻⁶",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.9),
    )
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    save(fig, outdir, "Figure2_occurrence_validation")


def climate_figure(rows, outdir):
    columns = [
        ("original_score_recomputed", "Original"),
        ("same_level1_score", "Same\nregion"),
        ("distance_matched_score_250km", "250 km"),
        ("distance_matched_score_500km", "500 km"),
        ("distance_matched_score_1000km", "1000 km"),
    ]
    values = [
        np.asarray([float(row[key]) for row in rows if row[key] != ""], dtype=float)
        for key, _ in columns
    ]
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    positions = np.arange(1, len(columns) + 1)
    ax.boxplot(values, positions=positions, widths=0.55, showfliers=False)
    for pos, vals in zip(positions, values):
        x = np.full(len(vals), pos, dtype=float)
        offsets = np.linspace(-0.12, 0.12, len(vals))
        ax.scatter(x + offsets, vals, s=16, alpha=0.65)
        ax.text(
            pos,
            min(1.03, np.max(vals) + 0.035),
            f"med={np.median(vals):.3f}",
            ha="center",
            fontsize=8,
        )
    ax.axhline(0.5, linestyle="--", linewidth=1)
    ax.set_ylim(0.2, 1.08)
    ax.set_xticks(positions, [label for _, label in columns])
    ax.set_ylabel("Climate-filtering score")
    ax.set_title("Climate-associated filtering persists after geographic matching")
    ax.text(
        0.02,
        0.03,
        "n=24 species for every sensitivity\nNeutral score = 0.5",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
    )
    fig.tight_layout()
    save(fig, outdir, "Figure4_climate_distance_sensitivity")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--occurrence-csv", type=Path, required=True)
    parser.add_argument("--occurrence-null-json", type=Path, required=True)
    parser.add_argument("--climate-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args=parser.parse_args()
    occurrence_figure(
        read_csv(args.occurrence_csv),
        json.loads(args.occurrence_null_json.read_text(encoding="utf-8")),
        args.output_dir,
    )
    climate_figure(read_csv(args.climate_csv), args.output_dir)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STRATA = ["1_family", "2_families", "3_to_5_families", "6plus_families"]
STRATA_LABELS = {
    "1_family": "1 family",
    "2_families": "2 families",
    "3_to_5_families": "3–5 families",
    "6plus_families": "6+ families",
}


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def average_ranks(values):
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def spearman(a, b):
    return float(np.corrcoef(average_ranks(a), average_ranks(b))[0, 1])


def jitter(label, width=0.12):
    import hashlib
    value = int(hashlib.sha256(label.encode("utf-8")).hexdigest()[:8], 16)
    return ((value / 0xFFFFFFFF) - 0.5) * 2 * width


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def panel_labels(axes):
    for label, ax in zip("abcde", np.ravel(np.atleast_1d(axes))):
        ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=12)


def boxplot_by_stratum(ax, rows, field):
    values = [
        [float(row[field]) for row in rows if row["host_breadth_stratum"] == stratum]
        for stratum in STRATA
    ]
    labels = [
        f"{STRATA_LABELS[stratum]}\n(n={len(values[i])})"
        for i, stratum in enumerate(STRATA)
    ]
    ax.boxplot(values, tick_labels=labels, showfliers=False)
    for i, stratum in enumerate(STRATA, start=1):
        for row in rows:
            if row["host_breadth_stratum"] != stratum:
                continue
            ax.scatter(
                i + jitter(row["species"], 0.17),
                float(row[field]),
                s=13,
                alpha=0.38,
                edgecolors="none",
            )


def figure1(anth, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    panel_labels(axes)
    native = np.asarray([int(r["native_resource_units"]) for r in anth], dtype=float)
    contemporary = np.asarray([int(r["contemporary_resource_units"]) for r in anth], dtype=float)

    axes[0].scatter(native, contemporary, s=17, alpha=0.42, edgecolors="none")
    lo = max(1, min(native.min(), contemporary.min()))
    hi = max(native.max(), contemporary.max())
    axes[0].plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Native host-resource breadth (WGSRPD3 units)")
    axes[0].set_ylabel("Contemporary host-resource breadth (WGSRPD3 units)")
    expanded = sum(int(r["introduced_added_units"]) > 0 for r in anth)
    total_native = int(native.sum())
    total_contemporary = int(contemporary.sum())
    gain = total_contemporary - total_native
    axes[0].set_title(f"{expanded}/{len(anth)} species gained resource opportunity")
    axes[0].text(
        0.04, 0.96,
        f"Species × region units\n{total_native:,} → {total_contemporary:,}\n"
        f"+{gain:,} (+{100*gain/total_native:.1f}%)",
        transform=axes[0].transAxes, va="top", fontsize=9,
        bbox={"boxstyle":"round,pad=0.3","facecolor":"white","alpha":0.85},
    )

    boxplot_by_stratum(axes[1], anth, "log_resource_expansion")
    axes[1].axhline(0, linestyle="--", linewidth=1)
    rho = spearman(
        [float(r["host_family_count"]) for r in anth],
        [float(r["log_resource_expansion"]) for r in anth],
    )
    axes[1].set_ylabel("Log proportional resource expansion")
    axes[1].set_title(f"Little proportional gradient across diet breadth\nSpearman ρ = {rho:.3f}")
    axes[1].tick_params(axis="x", rotation=18)
    fig.suptitle("Anthropogenic host redistribution expands butterfly resource geography", fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.94))
    save(fig, outdir / "Figure1_resource_expansion_v02")


def figure2(validation, null_result, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.6))
    panel_labels(axes)
    total = int(validation["outside_native_species_x_units"])
    recovered = int(validation["outside_native_species_x_units_recovered_by_introduced_hosts"])
    unrecovered = total - recovered

    axes[0].barh([0], [recovered], label="Recovered by introduced hosts")
    axes[0].barh([0], [unrecovered], left=[recovered], label="Still outside contemporary envelope")
    axes[0].set_xlim(0, total)
    axes[0].set_yticks([])
    axes[0].set_xlabel("Outside-native butterfly species × WGSRPD3 units")
    axes[0].set_title(f"{recovered}/{total} units recovered ({100*recovered/total:.1f}%)")
    axes[0].legend(frameon=False, fontsize=8, loc="lower center")

    names = ["Envelope-size\nnull", "Region-matched\nnull", "Observed"]
    means = [
        null_result["uniform_global_null"]["mean"],
        null_result["level1_composition_matched_null"]["mean"],
        recovered,
    ]
    lower = [
        null_result["uniform_global_null"]["q025"],
        null_result["level1_composition_matched_null"]["q025"],
        recovered,
    ]
    upper = [
        null_result["uniform_global_null"]["q975"],
        null_result["level1_composition_matched_null"]["q975"],
        recovered,
    ]
    x=np.arange(3)
    axes[1].scatter(x, means, s=55, zorder=3)
    for i in range(2):
        axes[1].vlines(i, lower[i], upper[i], linewidth=4, alpha=0.55)
    axes[1].axhline(recovered, linestyle="--", linewidth=1)
    axes[1].set_xticks(x, names)
    axes[1].set_ylabel("Recovered species × region units")
    axes[1].set_title("Overlap exceeds structural expectations\nregion-matched p = 5 × 10⁻⁶")
    axes[1].text(
        0.03, 0.97,
        "Null intervals: 95%",
        transform=axes[1].transAxes, va="top", fontsize=8
    )
    fig.suptitle("Introduced host geography aligns with independent butterfly occurrences", fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.94))
    save(fig, outdir / "Figure2_occurrence_validation_null_v02")


def figure3(ceiling, regional, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    panel_labels(axes)

    thresholds = []
    rhos = []
    for key, value in ceiling["restricted_native_breadth_thresholds"].items():
        thresholds.append(int(key))
        rhos.append(float(value["rho_family_vs_log_expansion"]))
    order=np.argsort(thresholds)
    thresholds=np.asarray(thresholds)[order]
    rhos=np.asarray(rhos)[order]
    axes[0].plot(thresholds, rhos, marker="o")
    axes[0].axhline(0, linestyle="--", linewidth=1)
    axes[0].scatter([369], [ceiling["rho_host_family_vs_log_expansion"]], marker="s", s=55)
    axes[0].text(369, ceiling["rho_host_family_vs_log_expansion"]+0.012, "Full panel", ha="right", fontsize=8)
    axes[0].set_xlim(85, 380)
    axes[0].set_ylim(-0.12, 0.15)
    axes[0].set_xlabel("Maximum native resource breadth retained")
    axes[0].set_ylabel("Spearman ρ: host-family breadth vs expansion")
    axes[0].set_title("Near-zero association persists below ceiling")

    regions = []
    values = []
    ns = []
    for region, item in regional["dominant_native_resource_region"].items():
        if "dominant_region" in item:
            d = item["dominant_region"]
            name = item.get("name", region)
            rho = d["rho_host_family_vs_log_expansion"]
            n = d["species"]
        else:
            d = item
            name = region
            rho = d.get("rho", d.get("rho_host_family_vs_log_expansion"))
            n = d["species"]
        if n < 8 or rho is None:
            continue
        regions.append(name)
        values.append(float(rho))
        ns.append(int(n))
    y=np.arange(len(regions))
    axes[1].barh(y, values)
    axes[1].axvline(0, linestyle="--", linewidth=1)
    axes[1].set_yticks(y, [f"{name} (n={n})" for name,n in zip(regions,ns)])
    axes[1].set_xlim(-0.30,0.30)
    axes[1].set_xlabel("Within-region Spearman ρ")
    axes[1].set_title("Most major geographic strata remain near zero")
    fig.suptitle("The absence of a broad-generalist advantage is robust to map support and geography", fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.94))
    save(fig, outdir / "Figure3_resource_expansion_robustness_v02")


def figure4(primary, same_region, distance, effect, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.7))
    panel_labels(axes)

    rows = list(primary["species"])
    x=[float(row["host_family_count"]) + jitter(row["species"],0.10) for row in rows]
    y=[float(row["climate_filtering_score"]) for row in rows]
    axes[0].scatter(x,y,s=27,alpha=0.65,edgecolors="none")
    axes[0].axhline(0.5,linestyle="--",linewidth=1)
    axes[0].set_xlabel("Larval host-plant families")
    axes[0].set_ylabel("Climate-filtering score")
    axes[0].set_title("Original independent cross-fit\nmedian 0.801; 23/24 > 0.5")

    labels=["Original","Same\nregion","250 km","500 km","1000 km"]
    medians=[
        0.8013787349014623,
        same_region["same_level1_median_score"],
        distance["distance_matched_250km"]["median_score"],
        distance["distance_matched_500km"]["median_score"],
        distance["distance_matched_1000km"]["median_score"],
    ]
    above=[
        23,
        same_region["same_level1_species_above_0_5"],
        distance["distance_matched_250km"]["species_above_0_5"],
        distance["distance_matched_500km"]["species_above_0_5"],
        distance["distance_matched_1000km"]["species_above_0_5"],
    ]
    axes[1].plot(np.arange(5),medians,marker="o")
    axes[1].axhline(0.5,linestyle="--",linewidth=1)
    axes[1].set_xticks(np.arange(5),labels)
    axes[1].set_ylim(0.45,0.9)
    axes[1].set_ylabel("Median filtering score")
    axes[1].set_title("Filtering persists after spatial controls")
    for i,(m,n) in enumerate(zip(medians,above)):
        axes[1].text(i,m+0.018,f"{n}/24",ha="center",fontsize=8)

    obs=float(effect["observed_partial_spearman"])
    lo,hi=map(float,effect["bootstrap_percentile_95_ci"])
    axes[2].errorbar([obs],[0],xerr=[[obs-lo],[hi-obs]],fmt="o",capsize=5)
    axes[2].axvline(0,linestyle="--",linewidth=1)
    mde=float(effect["approximate_minimum_detectable_absolute_partial_correlation"]["fisher_z_approximation"])
    axes[2].axvline(-mde,linestyle=":",linewidth=1)
    axes[2].axvline(mde,linestyle=":",linewidth=1)
    axes[2].set_xlim(-0.75,0.55)
    axes[2].set_yticks([])
    axes[2].set_xlabel("Partial Spearman ρ")
    axes[2].set_title("Predicted negative effect unsupported\n95% bootstrap interval")
    axes[2].text(
        0.04,0.08,
        f"Observed = {obs:.3f}\n95% CI [{lo:.3f}, {hi:.3f}]\n"
        f"Approx. 80% MDE |ρ| = {mde:.3f}",
        transform=axes[2].transAxes,fontsize=8,va="bottom"
    )
    fig.suptitle("Climate-associated filtering remains after resource and distance controls", fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.94))
    save(fig, outdir / "Figure4_climate_filtering_robustness_v02")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--anthropogenic-csv",type=Path,required=True)
    ap.add_argument("--climate-primary-json",type=Path,required=True)
    ap.add_argument("--occurrence-validation-json",type=Path,required=True)
    ap.add_argument("--occurrence-null-json",type=Path,required=True)
    ap.add_argument("--ceiling-json",type=Path,required=True)
    ap.add_argument("--regional-json",type=Path,required=True)
    ap.add_argument("--same-level1-json",type=Path,required=True)
    ap.add_argument("--distance-json",type=Path,required=True)
    ap.add_argument("--climate-effect-json",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    args=ap.parse_args()

    anth=read_csv(args.anthropogenic_csv)
    figure1(anth,args.outdir)
    figure2(
        read_json(args.occurrence_validation_json),
        read_json(args.occurrence_null_json),
        args.outdir,
    )
    figure3(read_json(args.ceiling_json),read_json(args.regional_json),args.outdir)
    figure4(
        read_json(args.climate_primary_json),
        read_json(args.same_level1_json),
        read_json(args.distance_json),
        read_json(args.climate_effect_json),
        args.outdir,
    )


if __name__=="__main__":
    main()

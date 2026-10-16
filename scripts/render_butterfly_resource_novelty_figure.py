#!/usr/bin/env python3
"""Reconstruct an exploratory butterfly macroecology Figure from frozen receipts.

This DOES NOT replace the submitted GEB figures. It distinguishes potential
native-resource overlap from unique newly opened larval host-resource regions.
Bars denote means of TWO simulation-chain medians, not confidence intervals.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def extract_results(no_climate: dict, climate: dict, local: dict) -> dict:
    """Check the original estimands and frozen reconstruction before plotting."""
    abc = no_climate["accounting"]
    H = int(abc["fixed_candidate_incidences"])
    A = int(abc["observed_incidences_absorbed_by_native_resource"])
    N = int(abc["observed_incidences_supporting_new_resource_regions"])
    R = int(abc["observed_redundant_incidences_within_new_resource_regions"])
    U = int(abc["observed_unique_added_butterfly_region_cells"])
    if (H, A, N, R, U) != (75346, 46387, 28959, 14406, 14553):
        raise ValueError("observed potential-resource identity has changed")
    if H != A + N or U != N - R:
        raise ValueError("host-resource incidence accounting does not balance")
    if (
        climate["baseline"]["introduced_host_butterfly_region_incidences"] != H
        or climate["baseline"]["introduced_incidences_in_native_supported_cells"] != A
        or climate["baseline"]["added_units"] != U
    ):
        raise ValueError("climate-null source identity mismatch")
    if (
        local["status"] != "QUALIFIED_EXPLORATORY_MCMC_SENSITIVITY"
        or local["qualification"]["status"] != local["status"]
        or local["study"]["spatial_caliper_km"] != 1500
        or local["study"]["fixed_candidate_incidences"] != H
        or local["observed"]["introduced_incidences_in_native_supported_cells"] != A
        or local["observed"]["added_units"] != U
        or not all(x["pass"] for x in local["qualification"]["checks"])
        or not local["qualification"]["between_chain"]["pass"]
    ):
        raise ValueError("long-chain spatial null is not qualified")
    if (
        climate["models"]["temperature_precipitation_2x2"]["status"]
        != "EVALUABLE"
    ):
        raise ValueError("climate-matched null is not evaluable")

    base = no_climate["chains"]
    matched = climate["models"]["temperature_precipitation_2x2"]["chains"]
    regional = local["chains"]
    if not (len(base) == len(matched) == len(regional) == 2):
        raise ValueError("exactly two independent runs required per null")
    if sorted(x["seed"] for x in base) != [20261015, 20261016]:
        raise ValueError("unconstrained seed identity mismatch")
    if sorted(x["seed"] for x in matched) != [20261015, 20261016]:
        raise ValueError("climate seed identity mismatch")
    if sorted(x["seed"] for x in regional) != [20261017, 20261018]:
        raise ValueError("long-chain seed identity mismatch")
    if any(
        x["status"] != "EXPLORATORY_CLIMATE_CONDITIONED_NULL_EVALUABLE"
        for x in matched
    ):
        raise ValueError("some climate chain failed support gate")
    controls = [
        {
            "name": "Plant family\n+ continent",
            "a_excess": [
                x["null_summary"]["introduced_incidences_in_native_supported_cells"][
                    "observed_minus_null_median"
                ]
                for x in base
            ],
            "u_deficit": [
                -x["null_summary"]["added_units"]["observed_minus_null_median"]
                for x in base
            ],
            "draws_per_chain": base[0]["replicates"],
        },
        {
            "name": "+ Climate\n(2×2 bins)",
            "a_excess": [
                x["absorbed"]["observed_minus_null_median"] for x in matched
            ],
            "u_deficit": [
                -x["unique_added"]["observed_minus_null_median"] for x in matched
            ],
            "draws_per_chain": 299,
        },
        {
            "name": "+ Local geography\n(1,500 km)",
            "a_excess": [
                x["observed_absorption_minus_null"] for x in regional
            ],
            "u_deficit": [
                -x["observed_unique_added_minus_null"] for x in regional
            ],
            "draws_per_chain": 99,
        },
    ]
    for control in controls:
        if not all(v > 0 for k in ("a_excess", "u_deficit") for v in control[k]):
            raise ValueError("control effects unexpectedly change sign")
    return {
        "focal_butterflies": 239,
        "resource_regions": 365,
        "total_candidate_incidences": H,
        "accounting": [
            {"label": "Within native-supported\nresource geography", "n": A},
            {"label": "Extra introduced hosts\nwithin novel resource cells", "n": R},
            {"label": "Unique new butterfly ×\nregion resource cells", "n": U},
        ],
        "controls": controls,
        "definition": (
            "Potential larval-resource geography reconstructed from known "
            "butterfly–host links; not local larval use, host abundance, "
            "butterfly colonization, competitive pressure, or fitness."
        ),
    }


def render(figure_data: dict, out_png: Path, out_svg: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12.5, 8.0), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2, ncols=2, height_ratios=[1.0, 1.15],
        left=0.105, right=0.965, top=0.835, bottom=0.15,
        wspace=0.28, hspace=0.50,
    )
    first = fig.add_subplot(gs[0, :])
    ax_a = fig.add_subplot(gs[1, 0])
    ax_u = fig.add_subplot(gs[1, 1])

    total = figure_data["total_candidate_incidences"]
    fig.suptitle(
        "Where introduced larval hosts add — or overlap — butterfly resources",
        fontsize=16, weight="semibold", y=0.955,
    )
    fig.text(
        0.105, 0.898,
        f"{figure_data['focal_butterflies']} butterflies  •  "
        f"{figure_data['resource_regions']} WGSRPD3 regions  •  "
        f"{total:,} introduced host × butterfly × region incidences",
        fontsize=10.5,
    )

    # A. Incidence units in each mutually exclusive accounting class.
    parts = figure_data["accounting"]
    xpos = 0
    bars = []
    for item in parts:
        pct = 100 * item["n"] / total
        bar = first.barh(0, pct, left=xpos, height=0.48, linewidth=0)
        bars.append(bar)
        first.text(
            xpos + pct/2, 0.0, f"{item['n']:,}\n({pct:.1f}%)",
            ha="center", va="center", fontsize=10, weight="semibold",
        )
        xpos += pct
    first.set_ylim(-0.58, 0.57)
    first.set_xlim(0, 100)
    first.set_yticks([])
    first.set_xlabel("Percent of introduced host × butterfly × region incidences")
    first.set_title(
        "A  |  Most mapped host exposure does not open a distinct new resource region",
        loc="left", fontsize=11, pad=13,
    )
    first.set_xticks([0, 20, 40, 60, 80, 100])
    first.spines[["top", "right", "left"]].set_visible(False)
    fig.legend(
        [b[0] for b in bars],
        [x["label"].replace("\n", " ") for x in parts],
        loc="upper center", bbox_to_anchor=(0.53, 0.64),
        ncol=3, frameon=False, fontsize=8.3,
    )

    # B-C. Changes in observed-minus-null median across the nested MAIN
    # set of structural controls. Each marker is a different RNG chain.
    models = figure_data["controls"]
    label = [x["name"] for x in models]
    for ax, metric, title, ylabel in (
        (
            ax_a, "a_excess",
            "B  |  Native-resource absorption above matched expectation",
            "Observed minus null (host incidences)",
        ),
        (
            ax_u, "u_deficit",
            "C  |  Fewer distinct novel resource regions than expected",
            "Null minus observed (butterfly × region cells)",
        ),
    ):
        means = [sum(m[metric])/2 for m in models]
        ax.bar(list(range(3)), means, width=0.55, alpha=0.52)
        for i, row in enumerate(models):
            ax.scatter([i-0.10, i+0.10], row[metric], s=33, zorder=5)
            ax.annotate(
                f"{min(row[metric]):.0f}–{max(row[metric]):.0f}",
                (i, max(row[metric])), xytext=(0, 9),
                textcoords="offset points", ha="center", fontsize=9,
            )
        ax.set_xticks(range(3), label)
        ax.set_title(title, fontsize=10.6, loc="left", pad=14)
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, max(means)*1.25)
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
        ax.set_axisbelow(True)

    fig.text(
        0.105, 0.070,
        "Accounting: one introduced-host incidence supplies the first support "
        "in each newly added butterfly×region cell; any further hosts count as extra supports.\n"
        "Bars: mean of two simulation-chain medians; points: individual chain medians "
        "(499 / 299 / 99 sampled worlds per chain across B and C).\n"
        "All comparisons are exploratory, post-hoc and conditional on HOSTS/WCVP; "
        "model differences are not causal fractions explained by climate or distance.",
        fontsize=8.2, va="center",
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=210, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--plant-null-json", type=Path, required=True)
    p.add_argument("--climate-null-json", type=Path, required=True)
    p.add_argument("--local-null-json", type=Path, required=True)
    p.add_argument("--output-png", type=Path, required=True)
    p.add_argument("--output-svg", type=Path, required=True)
    p.add_argument("--output-values-json", type=Path, required=True)
    args = p.parse_args()
    data = extract_results(
        load(args.plant_null_json),
        load(args.climate_null_json),
        load(args.local_null_json),
    )
    render(data, args.output_png, args.output_svg)
    args.output_values_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_values_json.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": "EXPLORATORY_FIGURE_RENDERED",
        "png": str(args.output_png),
        "svg": str(args.output_svg),
        "data": data,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

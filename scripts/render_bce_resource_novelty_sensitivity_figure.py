#!/usr/bin/env python3
"""Plot source dependence of introduced larval-resource novelty (not population).

One composition per source inventory. All numbers come from a frozen,
source-verified receipt, never from a manually rewritten list of plant hosts.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_validate(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("status") != "SUCCESS_POSTHOC_EUROPE_BIASED_INVENTORY_SENSITIVITY_NOT_SPATIAL_NULL":
        raise ValueError("unverified BCE source-sensitivity status")
    a = data["original_HOSTS"]
    b = data["BCE_Clarke_partial_augmented"]
    d = data["paired_source_change"]
    if any(c["focal_butterflies"] != 239 for c in (a, b)):
        raise ValueError("focal butterfly panel drift")
    expected = [
        (a, (75346, 46387, 14406, 14553)),
        (b, (103349, 67740, 20642, 14967)),
    ]
    for o, (H,A,R,U) in expected:
        if (
            o["candidate_incidences"], o["native_absorbed_incidences"],
            o["redundant_new_cell_incidences"], o["distinct_added_cells"]
        ) != (H,A,R,U) or H != A+R+U:
            raise ValueError("resource novelty accounting doesn't match frozen source")
    if (
        d["previously_unadded_regions_newly_added"] != 1817
        or d["old_added_cells_reclassified_native"] != 1403
        or d["delta_U"] != 414
        or 1817 - 1403 != 414
    ):
        raise ValueError("BCE-source reclassification identity drift")
    if data["source"]["butterflies_with_candidate_hosts"] != 81 or (
        data["source"]["candidate_additional_host_links"] != 1027
    ):
        raise ValueError("BCE candidate host sample unexpectedly changed")
    return {
        "schema": "chocho_butterfly_bce_inventory_novelty_figure_v0.1",
        "sample": {"focal_butterflies":239, "BCE_augmented_butterflies":81,
                   "candidate_links":1027},
        "scenarios": [
            {"label":"Original HOSTS inventory", "H":a["candidate_incidences"],
             "A":a["native_absorbed_incidences"],
             "R":a["redundant_new_cell_incidences"],
             "U":a["distinct_added_cells"]},
            {"label":"HOSTS + European BCE/Clarke candidate links",
             "H":b["candidate_incidences"],
             "A":b["native_absorbed_incidences"],
             "R":b["redundant_new_cell_incidences"],
             "U":b["distinct_added_cells"]},
        ],
        "newly_added_regions":1817,
        "reclassified_as_native":1403,
        "net_new_regions":414,
        "denominator_warning":"H differs by inventory; percentage shares are not true butterfly utilization probabilities, nor matched null significance",
    }


def render(source: dict, png: Path, svg: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    categories = [
        ("A", "Already native-supported host×butterfly×region incidence"),
        ("R", "Extra introduced supports in a new resource cell"),
        ("U", "First support of a uniquely new butterfly×region cell"),
    ]
    fig, ax = plt.subplots(figsize=(12.4, 5.9))
    fig.subplots_adjust(left=0.13, right=0.94, top=0.73, bottom=0.19)
    handles = []
    for idx, case in enumerate(source["scenarios"]):
        y = 1 - idx
        cursor = 0.0
        H = case["H"]
        for key, label in categories:
            percent = 100*case[key]/H
            bar=ax.barh(y, percent, left=cursor, height=0.41)
            if not handles or len(handles) < 3 and idx==0:
                handles.append(bar[0])
            ax.text(
                cursor + percent/2, y,
                f'{percent:.1f}%\n({case[key]:,})',
                ha="center", va="center", fontsize=10.4,
                weight="semibold",
            )
            cursor += percent
        if abs(cursor-100)>1e-8:
            raise ValueError("percentages do not sum to 100")
        ax.text(
            101.1, y, f'H = {H:,}',
            ha="left", va="center", fontsize=11, weight="medium"
        )
    ax.set_xlim(0, 113)
    ax.set_ylim(-0.65, 1.63)
    ax.set_yticks([1, 0], [x["label"] for x in source["scenarios"]])
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("Percentage of introduced host × butterfly × region candidate incidences")
    ax.grid(axis="x", alpha=0.10)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.suptitle(
        "How host-list completeness changes inferred novelty of butterfly resources",
        fontsize=17, fontweight="semibold", x=0.10, ha="left", y=0.97,
    )
    fig.text(
        0.10,0.898,
        "239 butterflies  |  81 source-augmented  |  1,027 candidate accepted-ID host links",
        fontsize=11, ha="left",
    )
    fig.legend(
        handles, [label for _, label in categories],
        loc="upper left", bbox_to_anchor=(0.10, 0.84),
        frameon=False, ncol=1, fontsize=8.5,
    )
    fig.text(
        0.12,0.085,
        "Source revision reclassified 1,403 previously introduced-added cells as natively supported; "
        "1,817 other cells became newly supplied (net +414).",
        fontsize=10.3,
    )
    fig.text(
        0.12,0.046,
        "European BCE/Clarke partial-source sensitivity, not a representative global correction or an independent spatial null.",
        fontsize=8.9,
    )
    png.parent.mkdir(parents=True, exist_ok=True)
    svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--receipt-json",type=Path,required=True)
    p.add_argument("--output-png",type=Path,required=True)
    p.add_argument("--output-svg",type=Path,required=True)
    p.add_argument("--output-data-json",type=Path,required=True)
    a=p.parse_args()
    validated=load_validate(a.receipt_json)
    render(validated,a.output_png,a.output_svg)
    a.output_data_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_data_json.write_text(json.dumps(validated,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"schema":validated["schema"],"source_sample":validated["sample"],
      "png":str(a.output_png),"svg":str(a.output_svg)},indent=2))


if __name__=="__main__":
    main()

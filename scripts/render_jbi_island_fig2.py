#!/usr/bin/env python3
"""Render JBI Figure 2: original island resources and matched-stratum contrast.

Reads ONLY three Git-blob-pinned published JSON source receipts. No model fitting.
The archival butterfly panel is a selected observation sample, NOT an estimate
of biological realization; no record is NOT an absence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = {
    "opportunity": (
        "results/island-mainland-resource-q1/island_mainland_summary.json",
        "9a7c4b2e434002d7125e328ba7f2e7b0b21ef2b5",
    ),
    "matched": (
        "results/island-mainland-resource-area-v02/area_adjusted_summary.json",
        "d74c651280d9790941b6bdfab14dc7ae3b200912",
    ),
    "gbif": (
        "results/q3-observed-butterfly-island-opportunity-v02/q3_archival_gbif_overlap_audit.json",
        "9240587123e05e5c7b0456a5d3404583bd105065",
    ),
}


def blob_sha(raw: bytes):
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load():
    data, checks = {}, {}
    for label, (path, sha) in SOURCE.items():
        raw = (ROOT / path).read_bytes()
        actual = blob_sha(raw)
        if sha != actual:
            raise ValueError(f"Source Git blob hash changed: {path}: {actual}")
        data[label] = json.loads(raw)
        checks[label] = {"path": path, "git_blob_sha": actual,
                         "sha256": hashlib.sha256(raw).hexdigest()}
    island = data["opportunity"]["groups"]["island"]
    mainland = data["opportunity"]["groups"]["mainland"]
    match = data["matched"]["stratum_test"]
    gbif = data["gbif"]["validation_panel"]

    groups = [island, mainland]
    for x in groups:
        assert x["native_butterfly_region_units"] + x["added_butterfly_region_units"] == x["contemporary_butterfly_region_units"]
        calc = 100 * x["added_butterfly_region_units"] / x["native_butterfly_region_units"]
        assert math.isclose(calc, x["pooled_increase_pct_if_native_positive"], abs_tol=1e-9)
        calc2 = 100 * x["added_butterfly_region_units"] / x["contemporary_butterfly_region_units"]
        assert math.isclose(calc2, x["pooled_added_share_pct_if_contemporary_positive"], abs_tol=1e-9)
    assert island["all_regions"] == 126 and mainland["all_regions"] == 237
    assert island["regions_new_from_zero_baseline"] == 10
    assert math.isclose(match["weighted_island_minus_mainland_share"], 0.10702159757115921, abs_tol=1e-12)
    assert match["within_stratum_permutation_n"] == 4999
    assert match["within_stratum_descriptive_p_two_sided"] == 0.0002
    n = gbif["sampled_species_added_island_opportunity_cells"]
    found = gbif["recorded_butterfly_in_new_opportunity_cells"]
    assert gbif["gbif_species_with_original_records"] == 31
    assert n == 920 and found == 13
    assert found < n
    assert math.isclose(found/n, gbif["fraction_sampled_added_opportunity_cells_with_any_gbif_record"], abs_tol=1e-12)
    return data, checks


def render(data, out):
    op = data["opportunity"]["groups"]
    isl, land = op["island"], op["mainland"]
    match = data["matched"]["stratum_test"]
    val = data["gbif"]["validation_panel"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.9),
                            gridspec_kw={"width_ratios": [1.43, 1.0, 1.05]})
    ax, bx, cx = axes
    fig.suptitle("Introduced plants expand potential butterfly host resources particularly on islands",
                 fontsize=15, fontweight="bold", y=0.985)
    ax.set_title("A  Reconstructed host opportunities\n239 butterflies; WGSRPD3 regions", loc="left")
    yy = [1, 0]
    native = np.array([isl["native_butterfly_region_units"], land["native_butterfly_region_units"]])
    added = np.array([isl["added_butterfly_region_units"], land["added_butterfly_region_units"]])
    ax.barh(yy, native, height=0.48, label="Native-only", color="#537492")
    ax.barh(yy, added, left=native, height=0.48, label="Added by introduced hosts", color="#d99b4e")
    for y, nat, ad, group in zip(yy, native, added, (isl, land)):
        ax.text(nat+ad+350, y, f"{nat+ad:,.0f}  (+{group['pooled_increase_pct_if_native_positive']:.1f}%)",
                va="center", fontsize=10.3)
        ax.text(nat/2, y, f"{nat:,.0f}", va="center", ha="center",
                color="white", fontsize=10.5, fontweight="bold")
        ax.text(nat+ad/2, y, f"+{ad:,.0f}", va="center", ha="center",
                fontsize=9.8, color="#241b10")
    ax.set_yticks(yy, ["Islands\n126 regions", "Mainland\n237 regions"])
    ax.set_xlim(0, 38000)
    ax.set_xlabel("Potential butterfly × botanical-region combinations")
    ax.legend(loc="upper left", bbox_to_anchor=(0,-0.20), frameon=False,
              fontsize=9, ncol=2)
    ax.grid(axis="x", alpha=0.14)
    ax.set_axisbelow(True)

    bx.set_title("B  Introduced-added share\nof contemporary opportunities", loc="left")
    share = np.array([isl["pooled_added_share_pct_if_contemporary_positive"],
                      land["pooled_added_share_pct_if_contemporary_positive"]])
    bx.barh(yy, share, height=0.5, color=["#537492", "#aeb7be"])
    for y,v in zip(yy,share):
        bx.text(v+1.1, y, f"{v:.1f}%", fontsize=11, va="center")
    bx.set_yticks(yy, ["Islands", "Mainland"])
    bx.set_xlim(0,65)
    bx.set_xlabel("Added / contemporary resource cells (%)")
    bx.grid(axis="x", alpha=0.14)
    bx.set_axisbelow(True)
    bx.text(0.04, -0.37,
            "Pooled added-share contrast is descriptive;\n"
            "the area-and-biogeography-matched contrast is in C.",
            transform=bx.transAxes, ha="left", va="top", fontsize=9.2)

    # Visualise the source-fixed overlap strata directly rather than
    # elevating the very limited GBIF index-overlap screen to a main panel.
    strata = match["strata_with_both_classes"]
    assert len(strata) == 15
    weights = np.array([z["weight"] for z in strata],dtype=float)
    deviations = 100*np.array([z["difference"] for z in strata],dtype=float)
    weighted = np.average(deviations,weights=weights)
    assert math.isclose(weighted,100*match["weighted_island_minus_mainland_share"],abs_tol=1e-10)
    yjit = np.linspace(-0.21,0.21,len(strata))
    cx.axvline(0,color="#505861",linewidth=1,linestyle="--",zorder=1)
    cx.axvline(weighted,color="#d59d55",linewidth=2.4,zorder=2)
    cx.scatter(deviations,yjit,s=38+28*weights,color="#537492",
               alpha=.84,edgecolors="white",linewidths=.7,zorder=3)
    cx.set_title("C  Island excess in matched strata\n15 Level 1 × area overlap strata",
                 loc="left")
    cx.set_xlim(-4,30); cx.set_ylim(-.34,.34); cx.set_yticks([])
    cx.set_xticks([0,10,20,30])
    cx.set_xlabel("Island − mainland added-share difference (pp)")
    cx.grid(axis="x",alpha=.13); cx.set_axisbelow(True)
    cx.text(weighted+.85,.295,f"Weighted +{weighted:.1f} pp",
            color="#805521",fontsize=10.2,va="top")
    cx.text(0.03,-0.37,
            "Each dot is one matched stratum; size ∝ minimum\n"
            "of island and mainland region counts. Descriptive\n"
            f"within-stratum permutation p = {match['within_stratum_descriptive_p_two_sided']:.4f}.",
            transform=cx.transAxes,fontsize=9.0,va="top")

    fig.text(0.5,0.052,
             "Resource cells are reconstructed host opportunities, not observed feeding or butterfly populations. "
             "Matched strata are descriptive, not independent causal replicates.",
             fontsize=9.4, ha="center")
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.40, top=0.79, wspace=0.44)
    out.mkdir(parents=True, exist_ok=True)
    hashes={}
    for ext in ("png","svg","pdf"):
        p=out/f"fig2_island_resource_expansion.{ext}"
        fig.savefig(p,dpi=180,bbox_inches="tight")
        hashes[ext]=hashlib.sha256(p.read_bytes()).hexdigest()
    plt.close(fig)
    return hashes


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument("--out-dir", type=Path, default=ROOT/"_jbi_fig2")
    args=cli.parse_args()
    data, checks=load()
    hashes=render(data,args.out_dir)
    report={
        "schema":"chocho_JBI_Fig2_blob_locked_source_render_v1",
        "status":"PASS_SOURCE_AND_NUMERIC_AUDIT",
        "source_files":checks,
        "output_sha256":hashes,
        "limits":["Not observed larval feeding or butterfly population growth.",
                  "The 13/920 sample describes indexed GBIF occurrence, not realisation of the full 239-butterfly resource network.",
                  "Matched +10.7 percentage points describes area-and-broad-biogeography overlap, distinct from unadjusted pooled shares."]
    }
    (args.out_dir/"fig2_receipt.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()

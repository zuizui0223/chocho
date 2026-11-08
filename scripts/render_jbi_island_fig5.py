#!/usr/bin/env python3
"""Render the source-locked Fig. 5 maritime-proxy specificity panels.

No regression, null shuffling, or route simulation is performed here.
The script independently checks original CSV values before drawing summaries.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript/island_resources_marine_access_20261108/FIG5_FROZEN_VALUES_V01.json"


def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def read_committed_csv(spec):
    path = ROOT / spec["path"]
    raw = path.read_bytes()
    assert git_blob_sha(raw) == spec["blob_sha"], f"changed original source: {path}"
    return list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))


def only(rows, **filters):
    matches = [row for row in rows if all(row.get(k) == v for k, v in filters.items())]
    assert len(matches) == 1, f"Expected one source row for {filters}; found {len(matches)}"
    return matches[0]


def eq(label, got, expected):
    assert math.isclose(float(got), float(expected), rel_tol=0, abs_tol=1e-9), (
        label, got, expected
    )


def audit_source(data):
    source = {k: read_committed_csv(v) for k, v in data["source_files"].items()}
    for item in data["latitude_controlled"]["samples"]:
        row = only(source["lat"], target="new_butterfly_jaccard", contrast=item["contrast"],
                   pair_class="all_pairs")
        for key, field in (("mse_gain", "MSE_reduction_percent"),
                           ("mae_gain", "MAE_reduction_percent"),
                           ("n_pairs", "unique_evaluated_pairs")):
            eq(f"latitude:{item['contrast']}:{key}", row[field], item[key])
    for spec, key in zip(data["spatial_placebo"]["sets"], ("latNull", "localNull")):
        row = only(source[key], target="new_butterfly_jaccard",
                   contrast="AIS_without_wind", pair_class="all_pairs")
        for source_key, target_key in [
            ("MSE_observed_gain_percent", "actual_mse_gain"),
            ("MSE_null_median_gain_percent", "median_mse_gain"),
            ("MSE_null_025_gain_percent", "interval025"),
            ("MSE_null_975_gain_percent", "interval975"),
            ("MSE_null_ge_observed_count", "placebos_at_least_as_good"),
            ("null_replicates", "n_placebos"),
            ("unique_test_pairs", "n_pairs")
        ]:
            eq(f"{key}:{source_key}", row[source_key], spec[target_key])
    hist = data["historic"]["rows"]
    h1 = only(source["hist"], scope="all17", target="new_butterfly_jaccard",
              contrast="BASE_to_MED_SUM")
    h2 = only(source["histDiag"], evaluation_scope="FULL17_TRAINED_NO9_TEST120",
              target="new_butterfly_jaccard", contrast="BASE_to_MED_SUM")
    h3 = only(source["hist"], scope="leave_original_landmass9_out",
              target="new_butterfly_jaccard", contrast="BASE_to_MED_SUM")
    for label, row, entry in zip(("full17", "same120-full17-trained", "same120-refit-no9"),
                                 (h1, h2, h3), hist):
        for name, field in (("mse_gain", "MSE_reduction_percent"),
                            ("mae_gain", "MAE_reduction_percent")):
            eq(f"{label}:{name}", row[field], entry[name])
        n_field = "original_source_unique_test_pairs" if label.startswith("same120-full17") else "test_unique_island_pairs"
        eq(f"{label}:pairs", row[n_field], entry["n_pairs"])
    eq("actual AIS source agreement", data["spatial_placebo"]["sets"][0]["actual_mse_gain"],
       data["spatial_placebo"]["sets"][1]["actual_mse_gain"])
    eq("actual AIS model agreement", data["spatial_placebo"]["sets"][0]["actual_mse_gain"],
       data["latitude_controlled"]["samples"][0]["mse_gain"])
    return {"source_blobs_verified": len(source), "latitude_scores": 3,
            "placebo_summaries": 2, "historical_scores": 3, "status": "PASS"}


def render(data, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12,
                         "axes.labelsize": 10, "savefig.dpi": 180,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.8),
                            gridspec_kw={"width_ratios": [1.04, 1.14, 1.38]})
    a, b, c = axes
    fig.suptitle("Can maritime accessibility be distinguished from geographic structure?",
                 fontweight="bold", fontsize=15, y=0.98)
    a.set_title("A  Modern marine proxies\nLatitude-adjusted, 1,275 pairs", loc="left")
    axdata = data["latitude_controlled"]["samples"]
    for y, item in zip((2, 1, 0), axdata):
        a.barh(y, item["mse_gain"], height=0.54)
        a.text(item["mse_gain"] + (0.14 if item["mse_gain"] >= 0 else -0.14),
               y, f"{item['mse_gain']:+.2f}%", ha="left" if item["mse_gain"] >= 0 else "right",
               va="center", fontsize=10)
    a.set_yticks([2, 1, 0], [x["label"] for x in axdata])
    a.axvline(0, color="black", lw=0.8)
    a.set_xlim(-6.9, 2.6)
    a.set_xlabel("Held-out MSE improvement (%)")
    a.grid(axis="x", alpha=0.2)
    a.text(0.02, -0.29, "Positive = improved prediction.\nAll comparisons use the same 51-island source cohort.",
           transform=a.transAxes, fontsize=8.4, va="top")

    b.set_title("B  AIS graph identity\n199 placebos per definition", loc="left")
    po = data["spatial_placebo"]["sets"]
    obs = po[0]["actual_mse_gain"]
    b.axvline(obs, linestyle="--", linewidth=1.3, label=f"Observed AIS {obs:+.2f}%")
    b.axvline(0, color="black", lw=0.7)
    for y, item in zip((1, 0), po):
        b.plot([item["interval025"], item["interval975"]], [y, y], lw=4,
               solid_capstyle="round", alpha=0.66)
        b.plot(item["median_mse_gain"], y, marker="D", markersize=7)
        b.text(-3.85, y + 0.22,
               f"{item['placebos_at_least_as_good']}/{item['n_placebos']} at least as good",
               fontsize=8.5)
    b.set_yticks([1, 0], ["Latitude-band\nrelabeling", "Local ≤750 km\nnode swapping"])
    b.set_xlim(-4.05, 4.6)
    b.set_ylim(-0.65, 1.56)
    b.set_xlabel("AIS MSE improvement (%)")
    b.grid(axis="x", alpha=0.2)
    b.legend(loc="lower right", frameon=False, fontsize=8)
    b.text(0.02, -0.29, "Line = 2.5–97.5% surrogate range;\ndiamond = surrogate median. Not a causal null.",
           transform=b.transAxes, fontsize=8.4, va="top")

    c.set_title("C  Historical liner connectedness\n2006–2010 index → 2011–2022 records", loc="left")
    cc = data["historic"]["rows"]
    short = ["All 17; 136 test pairs", "Full-17 trained;\n120 no-group-9 tests",
             "Refit without group 9;\nsame 120 tests"]
    for y, item in zip((2, 1, 0), cc):
        c.barh(y, item["mse_gain"], height=0.56)
        c.text(item["mse_gain"] + (0.22 if item["mse_gain"] >= 0 else -0.22), y,
               f"{item['mse_gain']:+.2f}%", ha="left" if item["mse_gain"] >= 0 else "right",
               va="center", fontsize=10)
    c.set_yticks([2, 1, 0], short)
    c.axvline(0, color="black", lw=0.8)
    c.set_xlim(-5.2, 11)
    c.set_xlabel("Held-out MSE improvement (%)")
    c.grid(axis="x", alpha=0.2)
    c.text(0.02, -0.29, "Train-set membership changes the sign.\nPLSCIQ is an index, not observed ship voyages.",
           transform=c.transAxes, fontsize=8.4, va="top")

    fig.text(0.5, 0.055, "Distinct source windows and baselines; percentages are not additive. "
             "No recorded plant carriage, butterfly arrival or breeding is identified.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=0.11, right=0.975, top=0.79, bottom=0.27, wspace=0.48)
    outputs = {}
    for extension in ("png", "svg", "pdf"):
        path = out / f"fig5_maritime_proxy_specificity.{extension}"
        fig.savefig(path, bbox_inches="tight", dpi=180)
        outputs[extension] = hashlib.sha256(path.read_bytes()).hexdigest()
    plt.close(fig)
    return outputs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", type=Path, default=ROOT / "_jbi_fig5")
    args = p.parse_args()
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert data["schema"] == "chocho_JBI_fig5_frozen_observational_proxy_comparison_v1"
    audited = audit_source(data)
    outputs = render(data, args.out_dir)
    receipt = {
        "schema": "chocho_JBI_fig5_source_render_receipt_v1",
        "frozen_data_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_checks": audited, "output_sha256": outputs,
        "limits": data["cautions"], "status": "SOURCE_AUDITED_FIGURE_RENDERED_NOT_CAUSAL_TRANSPORT"
    }
    (args.out_dir / "fig5_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "output": str(args.out_dir),
                      "audit": audited, "figure_sha256": outputs["png"]}, indent=2))


if __name__ == "__main__":
    main()

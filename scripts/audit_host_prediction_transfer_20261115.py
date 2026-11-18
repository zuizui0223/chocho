#!/usr/bin/env python3
"""Verify ecological breadth and transfer limits of the frozen 155-butterfly host forecast.

The source already contains retrospectively fitted models. This script checks
immutable reported results and derives comparable loss reductions without
re-fitting models or claiming past plant movements, occupancy or establishment.
All fixed-prediction intervals are descriptive (species × spatial resampling).
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

SOURCE = Path("provenance/reviewer_defenses/results/host_forecast_validation_20261029")
OUTPUT = Path("results/host_prediction_transfer_bridge_20261115")
FILES = {
    "metrics": "forecast_validation_all_metrics.csv",
    "pairs": "forecast_paired_comparisons.csv",
    "initially_unrecorded": "forecast_initially_unrecorded_comparison.csv",
    "scope": "forecast_initially_unrecorded_scope.json",
    "model_lock": "forecast_model_lock.json",
    "validator": "forecast_metric_validation.json",
}
MODE_NAMES = {
    "temporal": "known species, held 2019–2022 years",
    "spatial": "known species, held space and years",
    "spatial_buffer100": "spatial-fold subset >=100 km from all training cells",
    "species": "unseen butterfly species and years",
}


def csv_rows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as file:
        return list(csv.DictReader(file))


def exact(rows, **need):
    selected = [r for r in rows if all(r[k] == v for k, v in need.items())]
    if len(selected) != 1:
        raise ValueError(f"Expected one source row with {need}, found {len(selected)}")
    return selected[0]


def losses(metrics, pairs, mode):
    a = exact(metrics, mode=mode, model="M3_both")
    h = exact(metrics, mode=mode, model="M4_host_resources")
    c = exact(metrics, mode=mode, model="M5_host_composition")
    compare = exact(pairs, mode=mode, added="M5_host_composition",
                    reference="M3_both", metric="mse")
    component = exact(pairs, mode=mode, added="M5_host_composition",
                      reference="M4_host_resources", metric="mse")
    vals = {k: float(row["species_equal_mse"]) for k, row in
            (("base", a), ("host_total", h), ("host_family", c))}
    if abs((vals["host_family"]-vals["base"]) - float(compare["difference"])) > 1e-11:
        raise ValueError("M3 vs M5 model-score comparison changed")
    if abs((vals["host_family"]-vals["host_total"]) - float(component["difference"])) > 1e-11:
        raise ValueError("Host composition beyond host total disagreement")
    if (int(compare["species"]) != 155 or
        int(a["species"]) != 155 or int(compare["rows"]) != int(a["rows"])):
        raise ValueError("155-butterfly sample or test row identity changed")
    reduction = 100 * (vals["base"]-vals["host_family"])/vals["base"]
    ci = [-100*float(compare["ci_high"])/vals["base"],
          -100*float(compare["ci_low"])/vals["base"]]
    return {
        "interpretation": MODE_NAMES[mode],
        "n_butterflies": int(compare["species"]),
        "test_records": int(compare["rows"]),
        "species_improved_M3_to_M5": int(compare["species_improved"]),
        "M3_background_climate_landcover_species_equal_MSE": vals["base"],
        "M4_plus_host_resource_totals_MSE": vals["host_total"],
        "M5_plus_host_family_composition_MSE": vals["host_family"],
        "total_host_information_MSE_reduction_pct": reduction,
        "fixed_prediction_bootstrap_95pct_reduction_pct": ci,
        "host_resource_totals_MSE_reduction_pct": 100*(vals["base"]-vals["host_total"])/vals["base"],
        "additional_host_family_composition_MSE_reduction_pct_vs_M4":
            100*(vals["host_total"]-vals["host_family"])/vals["host_total"],
        "host_composition_beyond_total_95pct_difference_absolute_MSE":
            [float(component["ci_low"]), float(component["ci_high"])],
        "host_information_positive_95pct_fixed_prediction":
            float(compare["ci_high"]) < 0,
        "host_composition_beyond_total_positive_95pct_fixed_prediction":
            float(component["ci_high"]) < 0,
    }


def initially_unrecorded(source):
    assert source["status"] == "COMPLETE_FOR_AVAILABLE_MODES"
    assert "NOT never-recorded historically" in source["definition"]
    return {
        "temporal": "recording probability exactly 0 in 2015–2018, evaluated 2019–2022",
        "subsets": {},
    }


def run(root: Path, out_dir: Path):
    folder = root / SOURCE
    inp = {k: folder / path for k, path in FILES.items()}
    hashes = {k: hashlib.sha256(f.read_bytes()).hexdigest() for k, f in inp.items()}
    metrics = csv_rows(inp["metrics"])
    pairs = csv_rows(inp["pairs"])
    early = csv_rows(inp["initially_unrecorded"])
    scope = json.loads(inp["scope"].read_text())
    model_lock = json.loads(inp["model_lock"].read_text())
    validator = json.loads(inp["validator"].read_text())
    if model_lock["status"] != "LOCKED_BEFORE_MODEL_FIT" or validator["status"] != "PASS":
        raise ValueError("Frozen forecast evaluation lock or original validator missing")
    if "static host catalogue" not in model_lock["historical_knowledge_limit"]:
        raise ValueError("Historical plant map provenance changed")
    result = {mode: losses(metrics, pairs, mode) for mode in MODE_NAMES}
    assert {mode: result[mode]["test_records"] for mode in MODE_NAMES} == {
        "temporal": 208482, "spatial": 208482,
        "spatial_buffer100": 102728, "species": 208482
    }
    # Three independently locked viewpoints on the SAME 133,557 test records:
    early_report = initially_unrecorded(scope)
    for mode in ("temporal", "spatial", "species"):
        a = exact(early, mode=mode, metric="mse")
        b = exact(early, mode=mode, metric="mae")
        if (int(a["species"]), int(a["rows"])) != (154, 133557):
            raise ValueError("Initially-unrecorded fixed source changed")
        if int(b["rows"]) != 133557:
            raise ValueError("Initially-unrecorded MAE denominator changed")
        if abs(float(a["reference_M3"]) - float(exact(metrics, mode=mode, model="M3_both")[
            "initially_unrecorded_mse"])) > 1e-12:
            raise ValueError("Reused initially-unrecorded MSE source changed")
        early_report["subsets"][mode] = {
            "n_species": int(a["species"]), "rows": int(a["rows"]),
            "M3_MSE": float(a["reference_M3"]),
            "M5_MSE": float(a["full_M5"]),
            "MSE_relative_reduction_pct": float(a["relative_improvement_percent"]),
            "fixed_prediction_95pct_MSE_difference_CI": [
                float(a["ci_low"]), float(a["ci_high"])
            ],
            "M3_MAE": float(b["reference_M3"]),
            "M5_MAE": float(b["full_M5"]),
            "MAE_relative_reduction_pct": float(b["relative_improvement_percent"]),
            "fixed_prediction_95pct_MAE_difference_CI": [
                float(b["ci_low"]), float(b["ci_high"])
            ],
        }
    assert abs(result["temporal"]["total_host_information_MSE_reduction_pct"]-4.12017)<0.01
    assert abs(result["spatial"]["total_host_information_MSE_reduction_pct"]-2.995)<0.01
    assert abs(result["spatial_buffer100"]["total_host_information_MSE_reduction_pct"]-2.735)<0.02
    assert abs(early_report["subsets"]["temporal"]["MSE_relative_reduction_pct"]-5.2317)<0.001

    report = {
        "analysis": "155 known-butterfly host identity predicts recorded geography with a species-transfer boundary",
        "original_data_doi": "10.15468/dl.umqmsf",
        "inputs_sha256": hashes,
        "frozen_model_lock_date": model_lock["date_jst"],
        "primary_response": model_lock["response"],
        "n_eligible_butterflies": 155,
        "four_main_transfer_modes": result,
        "earlier_unrecorded_cell_test": early_report,
        "summary": {
            "positive_known_species_time_and_spatial_transfer": True,
            "no_positive_unseen_species_transfer": True,
            "spatial_host_family_composition_adds_beyond_total_host_resource": result[
                "spatial"]["host_composition_beyond_total_positive_95pct_fixed_prediction"],
            "observed_recordings_not_breeding_or_occupancy": True,
            "do_not_infer_temporal_introduced_host_colonization": True,
            "cross_study_239_resource_and_155_forecast_populations_distinct": True,
        },
        "inferential_boundary": [
            "Temporal evaluation uses contemporary static host geography and 2026 GBIF data extraction; NOT 2018 contemporaneously issued forecasts",
            "Original response is rarefied probability of indexed record in a species × cell, NOT occupancy or actual colonization",
            "Absence of a 2015–2018 focal record is not a historic absence nor a first island arrival",
            "The 100-km distance check is a post-hoc nested spatial subset, not an independently refitted model or replication",
            "M5 adds both host resources and family composition; association cannot be assigned exclusively to introduced host plants",
            "Bootstrap intervals condition on fixed predictions, not new-model or new-recording-survey uncertainty",
            "The same 155 butterflies across modes do not form independent replications",
            "Unseen-species prediction not improved by host geography, limiting transfer beyond studied taxa",
            "155 forecasting taxa and 239 resource reconstruction taxa use overlapping host catalogue but different data supports",
            "No host-history x wind x shipping x successful-recruitment pathway is directly tested",
        ],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "summary.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    with (out_dir/"transfer_contrasts.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["mode", "species", "test_records", "base_M3_mse",
                  "host_M4_mse", "host_family_M5_mse",
                  "M5_mse_reduction_pct", "low_pct", "high_pct",
                  "species_improved", "M5_vs_M4_pct"]
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for mode, val in result.items():
            writer.writerow({
                "mode":mode, "species":val["n_butterflies"],"test_records":val["test_records"],
                "base_M3_mse":val["M3_background_climate_landcover_species_equal_MSE"],
                "host_M4_mse":val["M4_plus_host_resource_totals_MSE"],
                "host_family_M5_mse":val["M5_plus_host_family_composition_MSE"],
                "M5_mse_reduction_pct":val["total_host_information_MSE_reduction_pct"],
                "low_pct":val["fixed_prediction_bootstrap_95pct_reduction_pct"][0],
                "high_pct":val["fixed_prediction_bootstrap_95pct_reduction_pct"][1],
                "species_improved":val["species_improved_M3_to_M5"],
                "M5_vs_M4_pct":val["additional_host_family_composition_MSE_reduction_pct_vs_M4"],
            })
    print(json.dumps({
        "source_population":155,
        "transfer":[
            {k:round(v["total_host_information_MSE_reduction_pct"],4),
             "improved_species":v["species_improved_M3_to_M5"],
             "composition_beyond_total":v[
                "host_composition_beyond_total_positive_95pct_fixed_prediction"]}
            for k,v in result.items()],
        "initial_unrecorded":{k: v["MSE_relative_reduction_pct"]
                              for k,v in early_report["subsets"].items()},
        "source_lock":"PASS"
    }, ensure_ascii=False))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--out-dir", type=Path, default=OUTPUT)
    args = p.parse_args()
    run(args.repo_root, args.out_dir)

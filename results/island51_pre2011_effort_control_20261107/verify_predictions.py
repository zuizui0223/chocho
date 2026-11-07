#!/usr/bin/env python3
"""Independent standard-library audit of pre-2011 botanical-record effort sensitivity predictions."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

MODELS = ("CR","HR","HCR","CRE","HRE","HCRE")
TARGETS = ("new_butterfly_jaccard","earlier_butterfly_jaccard")
COMPARISONS = {
    "CR_to_HCR":("CR","HCR"),"HR_to_HCR":("HR","HCR"),
    "CR_to_CRE":("CR","CRE"),"CRE_to_HCRE":("CRE","HCRE"),
    "HR_to_HRE":("HR","HRE"),"HRE_to_HCRE":("HRE","HCRE"),
    "HCR_to_HCRE":("HCR","HCRE"),
}

def audit(directory):
    directory = Path(directory)
    fold_values = defaultdict(list)
    pair_count = Counter()
    observed_pair = {}
    count = 0
    with (directory/"heldout_predictions.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            target, model, island = row["target"], row["model"], row["heldout_island_id"]
            if target not in TARGETS or model not in MODELS:
                raise ValueError("Unexpected target/model")
            idx = int(row["source_pair_index"])
            if idx < 0 or idx >= 1275:
                raise ValueError("Source pair index out of range")
            obs, pred = float(row["observed"]), float(row["predicted"])
            sq, ae = float(row["squared_error"]), float(row["absolute_error"])
            if not all(map(math.isfinite, (obs, pred, sq, ae))):
                raise ValueError("Nonfinite prediction")
            if not (0 <= obs <= 1 and 0 <= pred <= 1):
                raise ValueError("Jaccard range invalid")
            if not math.isclose(sq, (obs-pred)**2, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Incorrect saved squared error")
            if not math.isclose(ae, abs(obs-pred), rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Incorrect saved absolute error")
            key = (target, idx)
            if key in observed_pair and not math.isclose(observed_pair[key], obs, rel_tol=0, abs_tol=1e-12):
                raise ValueError("Original observation differs between endpoints")
            observed_pair[key] = obs
            pair_count[(target, model, idx)] += 1
            fold_values[(target, model, island)].append((sq, ae))
            count += 1
    if (count != 30600 or len(pair_count) != 2*6*1275 or
            any(v != 2 for v in pair_count.values()) or len(observed_pair) != 2*1275):
        raise ValueError("Incomplete endpoint-prediction support")
    fold = {}
    with (directory/"model_fold_losses.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            key = (row["target"], row["model"], row["heldout_island_id"])
            if key in fold or len(fold_values.get(key, ())) != 50:
                raise ValueError("Duplicated or incomplete island fold")
            if (int(row["train_pairs"]), int(row["test_pairs"])) != (1225, 50):
                raise ValueError("Fold train/test sizes changed")
            mse = math.fsum(v[0] for v in fold_values[key])/50
            mae = math.fsum(v[1] for v in fold_values[key])/50
            if not math.isclose(mse, float(row["MSE"]), abs_tol=1e-12, rel_tol=1e-12):
                raise ValueError("Fold MSE reconstruction mismatch")
            if not math.isclose(mae, float(row["MAE"]), abs_tol=1e-12, rel_tol=1e-12):
                raise ValueError("Fold MAE reconstruction mismatch")
            fold[key] = mse, mae
    if len(fold) != 612 or len(fold_values) != 612:
        raise ValueError("Expected 612 complete model folds")
    expected = {}
    for target in TARGETS:
        for name, (ref, aug) in COMPARISONS.items():
            ids = sorted({k[2] for k in fold if k[:2] == (target, ref)}, key=int)
            if len(ids) != 51 or set(ids) != {k[2] for k in fold if k[:2] == (target, aug)}:
                raise ValueError("Missing or unmatched island folds")
            a, b = [fold[(target, ref, i)] for i in ids], [fold[(target, aug, i)] for i in ids]
            mse_a, mse_b = math.fsum(x[0] for x in a)/51, math.fsum(x[0] for x in b)/51
            mae_a, mae_b = math.fsum(x[1] for x in a)/51, math.fsum(x[1] for x in b)/51
            expected[(target, name)] = {
                "MSE_reference": mse_a, "MSE_augmented": mse_b,
                "MSE_reduction_percent": 100*(mse_a-mse_b)/mse_a,
                "MAE_reference": mae_a, "MAE_augmented": mae_b,
                "MAE_reduction_percent": 100*(mae_a-mae_b)/mae_a,
                "islands_mse_improved": sum(x[0] > y[0] for x, y in zip(a,b)),
                "islands_mse_worse": sum(x[0] < y[0] for x, y in zip(a,b)),
                "islands_mse_tied": sum(x[0] == y[0] for x, y in zip(a,b)),
                "islands_mae_improved": sum(x[1] > y[1] for x, y in zip(a,b)),
                "islands_mae_worse": sum(x[1] < y[1] for x, y in zip(a,b)),
            }
    checked = set()
    with (directory/"comparison_summary.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            key = row["target"], row["comparison"]
            if key not in expected or key in checked:
                raise ValueError("Unexpected or duplicate comparison")
            checked.add(key)
            for column, value in expected[key].items():
                recorded = float(row[column]) if isinstance(value, float) else int(row[column])
                if not math.isclose(value, recorded, rel_tol=1e-11, abs_tol=1e-11):
                    raise ValueError("Comparison value differs: " + str((key, column)))
    if len(checked) != 14:
        raise ValueError("Missing comparisons")
    receipt = {"schema": "chocho_island51_pre2011_effort_independent_audit_v1",
               "status": "ALL_INDEPENDENT_ENDPOINT_FOLD_AND_SUMMARY_CHECKS_PASS",
               "predictions": count, "model_target_pair_endpoints_checked": len(pair_count),
               "source_target_observations_checked": len(observed_pair),
               "folds_checked": len(fold), "comparisons_checked": len(checked)}
    (directory/"independent_audit.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("PRE2011_EFFORT_AUDIT_JSON="+json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, required=True)
    audit(parser.parse_args().results_dir)

#!/usr/bin/env python3
"""Audited prospective prediction of E. editha oviposition, NOT evolved hysteresis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

EXPECTED_MD5 = "6d9770cce9cacfccd98f7d19c67d0ec5"
URLS = [
    "https://zenodo.org/records/4318182/files/oviposition.csv?download=1",
    "https://zenodo.org/api/records/4318182/files/oviposition.csv/content",
]
HOSTS = ["CAHI", "CALE", "PLLA"]


def get_original():
    problems = []
    for url in URLS:
        try:
            req = Request(url, headers={"User-Agent": "Mozilla/5.0 (academic source integrity audit)"})
            with urlopen(req, timeout=25) as r:
                data = r.read(500000)
            got = hashlib.md5(data).hexdigest()
            if got != EXPECTED_MD5:
                problems.append({"url": url, "error": "source MD5 mismatch", "got_md5": got})
                continue
            return data, url, problems
        except (HTTPError, URLError, TimeoutError, OSError) as e:
            problems.append({"url": url, "error": str(e)[:200]})
    return None, None, problems


def parse_rows(raw: bytes):
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    expected = {"matriline", "date", "choice.number", "CAHI.eggs", "CALE.eggs", "PLLA.eggs", "other.eggs",
                "CAHI.cl", "CALE.cl", "PLLA.cl"}
    if not rows or not expected.issubset(rows[0]):
        raise ValueError("source columns not as expected")
    output = []
    for i, row in enumerate(rows):
        try:
            rec = {"source_row": i + 2, "matriline": row["matriline"].strip(),
                   "date": datetime.strptime(row["date"].strip(), "%m/%d/%Y").date(),
                   "choice_num": int(row["choice.number"])}
            for h in HOSTS + ["other"]:
                rec[h] = int(row[f"{h}.eggs"])
                if rec[h] < 0:
                    raise ValueError("negative egg count")
            for h in HOSTS:
                rec[f"{h}_cl"] = int(row[f"{h}.cl"])
                if rec[f"{h}_cl"] < 0:
                    raise ValueError("negative egg-clutch count")
            if not rec["matriline"]:
                raise ValueError("missing matriline")
            rec["host_eggs"] = sum(rec[h] for h in HOSTS)
            rec["introduced_share"] = rec["PLLA"] / rec["host_eggs"] if rec["host_eggs"] else None
            output.append(rec)
        except (ValueError, TypeError) as e:
            raise ValueError(f"bad original row {i+2}: {e}") from e
    # Every matriline/choice index represents a single exposure; preserve all source rows,
    # but repeated index would invalidate our within-female chronology.
    pairs = [(r["matriline"], r["choice_num"]) for r in output]
    if len(set(pairs)) != len(pairs):
        raise ValueError("repeated matriline × choice index")
    return output


def sequence_dataset(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[r["matriline"]].append(r)
    min_date = min(r["date"] for r in rows)
    trials = []
    future_trials = []
    for key, rr in groups.items():
        ordered = sorted(rr, key=lambda r: (r["choice_num"], r["date"]))
        # The previous egg-producing trial must occur strictly earlier in trial
        # sequence; empty trials are not "previous choices".
        for j, r in enumerate(ordered):
            if not r["host_eggs"]:
                continue
            preceding = [z for z in ordered[:j] if z["host_eggs"]]
            following = [z for z in ordered[j+1:] if z["host_eggs"]]
            if not preceding:
                continue
            prior = preceding[-1]
            x0 = [(r["date"] - min_date).days, r["choice_num"]]
            record = {"group": key, "row": r["source_row"],
                      "y": int(r["PLLA"] > r["CAHI"] + r["CALE"]),
                      "base": x0, "lag": prior["introduced_share"],
                      "other_eggs": r["other"], "n_host_eggs": r["host_eggs"]}
            trials.append(record)
            if following:
                future_trials.append({**record, "future": following[0]["introduced_share"]})
    return trials, future_trials


def predictions(items, feature):
    groups = sorted(set(r["group"] for r in items))
    y = np.array([r["y"] for r in items], dtype=int)
    x = np.array([r["base"] + ([r[feature]] if feature else []) for r in items], dtype=float)
    preds = np.full(len(items), np.nan)
    for group in groups:
        train = np.array([r["group"] != group for r in items])
        test = ~train
        if not train.any() or not test.any():
            raise ValueError("invalid leave-one-matriline-out partition")
        if len(set(y[train])) == 1:
            preds[test] = (sum(y[train]) + 0.5) / (len(y[train]) + 1)
        else:
            model = make_pipeline(StandardScaler(), LogisticRegression(C=1, max_iter=2000, solver="lbfgs"))
            model.fit(x[train], y[train])
            preds[test] = model.predict_proba(x[test])[:, 1]
    return np.clip(preds, 1e-8, 1-1e-8)


def score(items, a, b):
    ys = np.array([r["y"] for r in items])
    la = -(ys*np.log(a)+(1-ys)*np.log(1-a))
    lb = -(ys*np.log(b)+(1-ys)*np.log(1-b))
    return {"baseline_logloss": float(la.mean()), "augmented_logloss": float(lb.mean()),
            "logloss_improvement": float((la-lb).mean()),
            "brier_improvement": float((np.square(a-ys)-np.square(b-ys)).mean())}, (la-lb)


def cluster_ci(items, differences, draws=2000, seed=20261010):
    groups = sorted(set(r["group"] for r in items))
    idx = {g: np.array([i for i,r in enumerate(items) if r["group"] == g]) for g in groups}
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(draws):
        sample = rng.choice(groups, size=len(groups), replace=True)
        ii = np.concatenate([idx[g] for g in sample])
        values.append(float(np.mean(differences[ii])))
    return list(map(float, np.quantile(values, [0.025, 0.975])))


def run(rows):
    trials, fut = sequence_dataset(rows)
    groups = set(r["group"] for r in trials)
    positives = sum(r["y"] for r in trials)
    result = {"n_original_trials": len(rows), "n_original_matrilines": len({r["matriline"] for r in rows}),
              "n_eggs_named_hosts": sum(r["host_eggs"] for r in rows),
              "n_other_surface_eggs": sum(r["other"] for r in rows),
              "n_sequential_trials": len(trials), "n_sequential_matrilines": len(groups),
              "sequential_positive": positives,
              "same_trial_clutch_count_used_as_predictor": False,
              "outcome":"PLLA eggs strictly exceed CAHI + CALE eggs",
              "causal_evolutionary_hysteresis_estimated": False}
    if len(trials) < 25 or len(groups) < 8 or min(positives, len(trials)-positives) < 5:
        result["decision"] = "NOT_ESTIMABLE_FROZEN_FEASIBILITY_GATE"
        return result
    a = predictions(trials, None)
    b = predictions(trials, "lag")
    metric, dd = score(trials,a,b)
    result.update(metric)
    result["matriline_bootstrap_95ci"] = cluster_ci(trials,dd)
    result["n_forward_and_reverse_trials"] = len(fut)
    if len(fut) >= 25 and len(set(r["group"] for r in fut)) >= 8 and len({r["y"] for r in fut}) == 2:
        baseline = predictions(fut,None)
        past = predictions(fut,"lag")
        future = predictions(fut,"future")
        pm, _ = score(fut,baseline,past)
        fm, _ = score(fut,baseline,future)
        result["past_skill_same_subset"] = pm["logloss_improvement"]
        result["future_placebo_skill_same_subset"] = fm["logloss_improvement"]
        result["directional_signal"] = bool(pm["logloss_improvement"] > fm["logloss_improvement"])
    else:
        result["directional_signal"] = None
    result["sensitivity_other_eggs"] = {
        "n_trials_without_other_eggs": sum(r["other_eggs"] == 0 for r in trials),
        "note": "Report support; do not select fitted sample after outcome."
    }
    ci=result["matriline_bootstrap_95ci"]
    result["decision"] = ("PREDICTIVE_CARRYOVER_PILOT_PASSES_NONCAUSAL_GATE"
                          if ci[0] > 0 and result["directional_signal"] is True
                          else "NO_ROBUST_PREDICTIVE_CARRYOVER_GATE_PASS")
    result["claim_boundary"] = ("The 2021 paper already reported temporal preference changes. "
        "Prediction does not identify memory or heritable host lock-in. Only one population and no host-removal or adult survival response.")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receipt", type=Path, required=True)
    args = ap.parse_args()
    raw, source, errors = get_original()
    receipt = {"protocol":"docs/exploratory/EEDITHA_SEQUENTIAL_OVIPOSITION_PROTOCOL_V01.json",
               "source":source, "source_tries":errors, "source_md5_verified":bool(raw)}
    if raw is None:
        receipt.update({"status":"SOURCE_ACCESS_BLOCKED","effect_estimated":False})
    else:
        try:
            dat = parse_rows(raw)
            receipt.update({"status":"SOURCE_VERIFIED_ANALYSIS_COMPLETE", "md5":hashlib.md5(raw).hexdigest(),
                            "effect_estimated":True, "analysis":run(dat)})
        except ValueError as err:
            receipt.update({"status":"SOURCE_SCHEMA_OR_SEQUENCE_GATE_FAILED","effect_estimated":False,
                            "error":str(err)})
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2,allow_nan=False))


if __name__=="__main__":
    main()

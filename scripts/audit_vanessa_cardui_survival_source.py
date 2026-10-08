#!/usr/bin/env python3
"""Independent already-published Vanessa cardui larval-survival SOURCE gate.

Never interpret alive/censored records as completed adult emergence.
The host effects and larval source effects were reported by Saldivar & Wilson Rankin (2024).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

REQUIRED = (
    "time", "status", "species", "source", "pot", "trial.no",
    "plant.trial", "year"
)
PLANTS = {
    "ABPA": "Abutilon palmeri",
    "MACL": "Malacothamnus clementinus",
    "MAFA": "Malacothamnus fasciculatus",
    "NIGL": "Nicotiana glauca",
    "SPAM": "Sphaeralcea ambigua",
}
SOURCES = [
    "https://zenodo.org/records/10689711/files/larv.plants.csv?download=1",
    "https://zenodo.org/api/records/10689711/files/larv.plants.csv/content",
]
HORIZON = 10.0


def retrieve(expected_md5: str):
    errors = []
    for url in SOURCES:
        for attempt in range(2):
            try:
                req = urllib.request.Request(
                    url, headers={
                        "User-Agent": "chocho-butterfly-data-source-audit/0.1",
                        "Accept": "text/csv,text/plain,application/octet-stream,*/*",
                    }
                )
                with urllib.request.urlopen(req, timeout=40) as resp:
                    raw = resp.read()
                observed = hashlib.md5(raw).hexdigest()
                if observed != expected_md5:
                    raise RuntimeError(
                        "Original dataset digest mismatch: "
                        + observed + " != " + expected_md5
                    )
                return raw, {"url": url, "md5": observed, "bytes": len(raw)}
            except Exception as exc:
                errors.append(
                    f"{url}, attempt {attempt+1}: "
                    f"{type(exc).__name__}: {str(exc)[:120]}"
                )
                time.sleep(1)
    raise RuntimeError("The original source could not be verified: " + repr(errors))


def km_at(rows, horizon):
    """Event=death/status 2; status 1 right-censored, not adult eclosion."""
    if not rows:
        return {"km_survival_day10": None, "reason": "empty"}
    n = len(rows)
    events = Counter()
    censored = Counter()
    for row in rows:
        if row["event"]:
            events[row["time"]] += 1
        else:
            censored[row["time"]] += 1
    risk = n
    surv = 1.0
    for t in sorted(set(events) | set(censored)):
        if t > horizon:
            break
        d = events[t]
        if d > risk:
            raise RuntimeError("KM risk-set accounting failed")
        if d:
            surv *= 1.0 - d / risk
        risk -= d + censored[t]
    # If all subjects were censored before day 10, day-10 survival is unknown.
    latest = max(row["time"] for row in rows)
    if latest < horizon and surv > 0:
        return {
            "km_survival_day10": None,
            "reason": "follow-up_ends_before_day10_without_complete_failure",
        }
    return {
        "km_survival_day10": round(surv, 6),
        "at_risk_just_after_day10": risk,
        "reason": "observed_or_estimable",
    }


def analyze(raw, protocol):
    read = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if not read:
        raise RuntimeError("Original CSV empty")
    if not set(REQUIRED).issubset(read[0]):
        raise RuntimeError("Required columns absent: " + repr(sorted(set(REQUIRED) - set(read[0]))))
    natural = set(protocol["inclusion"]["plant_codes"])
    if natural != set(PLANTS):
        raise RuntimeError("The prespecified exact plant classes drifted")
    if set(protocol["inclusion"]["controls"]) != {"ARTI", "ARTH", "THAR"}:
        raise RuntimeError("Control classes drifted")

    rows = []
    dropped = Counter()
    raw_classes = Counter()
    for entry in read:
        host = str(entry["species"]).strip()
        raw_classes[host] += 1
        status = str(entry["status"]).strip()
        try:
            duration = float(entry["time"])
        except (ValueError, TypeError):
            dropped["bad_time"] += 1
            continue
        if not math.isfinite(duration) or duration < 0:
            dropped["bad_time"] += 1
            continue
        if status not in {"1", "2"}:
            dropped["unknown_event_status"] += 1
            continue
        rows.append({
            "host": host,
            "time": duration,
            "event": int(status == "2"),
            "source": str(entry["source"]).strip(),
            "trial": str(entry["trial.no"]).strip(),
            "plant_trial": str(entry["plant.trial"]).strip(),
            "pot": str(entry["pot"]).strip(),
            "year": str(entry["year"]).strip()
        })
    if dropped:
        raise RuntimeError("Fail closed on unrecognized original time/status: " + repr(dict(dropped)))
    if sum(raw_classes.values()) != len(rows):
        raise RuntimeError("Unexpected row attrition")
    unknown = set(raw_classes) - set(PLANTS) - {"ARTI","ARTH","THAR"}
    if unknown:
        raise RuntimeError("Unknown host class, do not silently relabel: " + repr(unknown))
    strata = defaultdict(list)
    source_counts = Counter()
    for row in rows:
        source_counts[row["source"]] += 1
        if row["host"] in natural:
            strata[(row["host"],row["source"])].append(row)
    plant_source_table = []
    for (host, source), rr in sorted(strata.items()):
        death = sum(r["event"] for r in rr)
        plant_source_table.append({
            "plant_code":host,
            "plant_species":PLANTS[host],
            "larva_source":source,
            "n":len(rr),
            "reported_death_events":death,
            "right_censored_or_alive":len(rr)-death,
            "time_min":min(r["time"] for r in rr),
            "time_max":max(r["time"] for r in rr),
            "distinct_trial_labels":len({r["trial"] for r in rr}),
            "distinct_plant_trial_labels":len({r["plant_trial"] for r in rr}),
            **km_at(rr, HORIZON)
        })
    if not plant_source_table:
        raise RuntimeError("No usable original natural-plant records")
    return {
        "status":"ORIGINAL_SOURCE_VERIFIED_EXPLORATORY_CALIBRATION_NOT_NEW_DISCOVERY",
        "butterfly":"Vanessa cardui",
        "study":"Saldivar and Wilson Rankin 2024 DOI 10.1002/ecs2.4810",
        "source_rows":len(read),
        "rows_retained":len(rows),
        "raw_diet_counts":dict(sorted(raw_classes.items())),
        "larva_source_counts":dict(sorted(source_counts.items())),
        "plant_source_strata":plant_source_table,
        "limits":[
            "status=1 may represent censoring/alive, never adult eclosion",
            "Source origin and host plant are nonrandomized and potentially confounded",
            "Multiple larvae in the same plant/trial need cluster-aware inference",
            "Experimental host-quality contrasts were already published",
            "No direct geographic colonization, regional resource prevalence or realized global fitness"
        ],
        "promotion":"Do NOT promote as an independent ecology paper or modify GEB main manuscript."
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--protocol",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    p=json.loads(args.protocol.read_text(encoding="utf-8"))
    if p["schema"]!="chocho_cardui_larval_survival_source_gate_v0.1":
        raise RuntimeError("Unexpected source protocol")
    raw, source=retrieve(p["source"]["md5"])
    result=analyze(raw,p)
    result["source_provenance"]=source
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "source_rows":result["source_rows"],
        "diet_counts":result["raw_diet_counts"],
        "source_counts":result["larva_source_counts"],
        "plant_source_strata":result["plant_source_strata"],
        "source_md5":source["md5"],
        "warning":result["limits"]
    },indent=2),flush=True)

if __name__=="__main__":
    main()

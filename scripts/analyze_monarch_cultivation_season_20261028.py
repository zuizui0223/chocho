#!/usr/bin/env python3
"""Within-native-host source-recorded monarch larva winter/cultivation contrasts.

This is a *photo host-composition* study, not a host-availability experiment.
All raw iNaturalist/author records are transient memory only; output contains
aggregate counts and explicitly conditional observational estimands.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import re
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from audit_monarch_annotated_host_source_20261028 import source_download, SHARE_URL
from audit_monarch_host_label_categories_20261028 import ORIGINAL_ZIP_SHA256, THINNED

OUTPUT_SCHEMA = "chocho_monarch_native_host_cultivation_winter_MH_20261028_v01"
WINTER = {11, 12, 1, 2}
SHOULDER = {3, 4, 9, 10}
SOURCE_ROWS = 12942
SEED = 20261028


def exact_binomial(name):
    """Exclude unresolved genus/complex/multispecies host labels."""
    return bool(re.fullmatch(r"[A-Z][a-z]+ [a-z][a-z-]*", name))


def read_source(archive):
    """Read only the author's spatial/week/host-thinned record set."""
    with archive.open(THINNED) as handle:
        reader = csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8-sig", newline=""))
        fields = set(reader.fieldnames or [])
        required = {"State", "year", "month", "host", "confidence", "cultivated",
                    "rearing", "native_status", "setting", "EW", "hex_id", "week"}
        if not required <= fields:
            raise ValueError(f"Source schema drift: {sorted(required-fields)}")
        items=[]
        quality=Counter()
        for row in reader:
            quality["source_rows"] += 1
            try:
                year = int(str(row["year"]).strip())
                month = int(str(row["month"]).strip())
            except (ValueError, TypeError):
                quality["invalid_calendar"] += 1
                continue
            if not 1 <= month <= 12 or not 1980 <= year <= 2025:
                quality["invalid_calendar"] += 1
                continue
            host = str(row["host"] or "").strip()
            state = str(row["State"] or "").strip()
            hexagon = str(row["hex_id"] or "").strip()
            confidence = str(row["confidence"] or "").strip().lower()
            cultivated = str(row["cultivated"] or "").strip().lower()
            rearing = str(row["rearing"] or "").strip().lower()
            origin = str(row["native_status"] or "").strip().lower()
            setting = str(row["setting"] or "").strip().lower()
            eastwest = str(row["EW"] or "").strip().lower()
            if not state or not hexagon or not exact_binomial(host):
                quality["missing_region_or_unresolved_host"] += 1
                continue
            if confidence not in {"certain", "high", "medium"}:
                quality["unqualified_host_confidence"] += 1
                continue
            if rearing != "no":
                quality["source_reared_or_rearing_ambiguous"] += 1
                continue
            if eastwest == "na" or setting == "obscured":
                quality["source_geography_ambiguous_or_obscured"] += 1
                continue
            if origin not in {"native", "non-native"}:
                quality["host_origin_unknown"] += 1
                continue
            if cultivated not in {"yes", "no", "likely", "unlikely", "maybe"}:
                quality["source_cultivation_missing"] += 1
                continue
            quality["common_source_eligible"] += 1
            items.append({
                "state":state, "year":year, "hex":hexagon, "host":host,
                "origin":origin, "cultivated":cultivated,
                "confidence":confidence, "winter":month in WINTER,
                "month":month,
            })
    if quality["source_rows"] != SOURCE_ROWS:
        raise ValueError("Author-thinned monarch source number changed")
    return items,dict(quality)


def table_strata(rows, *, origin="native", label_mode="strict",
                 keys=("host", "state", "year"), shoulder=False,
                 comparison="cultivation", high_confidence=False):
    groups = defaultdict(lambda: Counter())
    qc = Counter()
    for row in rows:
        if high_confidence and row["confidence"]=="medium":
            continue
        if shoulder and row["month"] not in WINTER | SHOULDER:
            continue
        if comparison=="cultivation":
            if origin is not None and row["origin"] != origin:
                continue
            if label_mode=="strict":
                target = {"yes":1,"no":0}.get(row["cultivated"])
            elif label_mode=="expanded":
                target = {"yes":1,"likely":1,"no":0,"unlikely":0}.get(row["cultivated"])
            else:
                raise ValueError("Unknown cultivation label mode")
        elif comparison=="origin":
            if row["cultivated"]!="no":
                continue
            target={"non-native":1,"native":0}.get(row["origin"])
        else:
            raise ValueError("Unknown comparison")
        if target is None:
            qc["ambiguous_treatment_labels"]+=1
            continue
        group=tuple(row[k] for k in keys)
        groups[group][(target, int(row["winter"]))]+=1
        qc["source_qualified_records"]+=1
    strata=[]
    for key, c in sorted(groups.items()):
        a,b,c_,d=(c[(1,1)],c[(1,0)],c[(0,1)],c[(0,0)])
        n=a+b+c_+d
        if min(a+b,c_+d)==0:
            continue
        if min(a+c_,b+d)==0:
            continue
        keyinfo=dict(zip(keys,key))
        strata.append({"a":a,"b":b,"c":c_,"d":d, **keyinfo, "n":n})
    return strata, dict(qc)


def mh_odds_ratio(strata):
    """Mantel-Haenszel no pseudo-count; a=treated winter, c=control winter."""
    top=bottom=0.
    for s in strata:
        a,b,c,d=(s[k] for k in "abcd")
        n=a+b+c+d
        if n<=0: raise ValueError("Empty stratum")
        top += a*d/n
        bottom += b*c/n
    if top<=0 or bottom<=0:
        return None
    return top/bottom


def group_bootstrap(strata, iterations=999):
    if not strata: return {"valid_draws":0,"interval":None}
    by_state=defaultdict(list)
    for s in strata:
        if "state" not in s: return {"valid_draws":0,"interval":None}
        by_state[s["state"]].append(s)
    states=sorted(by_state)
    rng=random.Random(SEED)
    estimates=[]
    for _ in range(iterations):
        sampled=[states[rng.randrange(len(states))] for _ in states]
        odds=mh_odds_ratio([s for state in sampled for s in by_state[state]])
        if odds is not None and math.isfinite(odds):
            estimates.append(math.log(odds))
    estimates.sort()
    if len(estimates)<math.ceil(.95*iterations):
        return {"valid_draws":len(estimates),"interval":None,
                "gate":"INCOMPLETE_BOOTSTRAP_FREQUENT_SEPARATION"}
    def percentile(q):
        idx=(len(estimates)-1)*q
        lo=int(idx); hi=min(len(estimates)-1,lo+1)
        return math.exp(estimates[lo]*(hi-idx)+estimates[hi]*(idx-lo))
    return {"valid_draws":len(estimates),
            "interval":[percentile(.025),percentile(.975)],
            "gate":"DESCRIPTIVE_STATE_CLUSTER_ONLY_NOT_CAUSAL"}


def summary(strata, *, keys, require_state_bootstrap=True, iterations=999):
    count=Counter()
    by_state=Counter()
    taxa=set()
    states=set()
    for s in strata:
        for k in "abcd": count[k]+=s[k]
        by_state[s["state"]]+=s["n"]
        states.add(s["state"])
        if "host" in s: taxa.add(s["host"])
    n=sum(count.values())
    strict_gate = (len(strata)>=20 and len(taxa)>=10 and len(states)>=10
                   and count["a"]>=15 and count["c"]>=15 and
                   (not by_state or max(by_state.values())/n<.5))
    odds=mh_odds_ratio(strata)
    return {
        "candidate_informative_strata":len(strata),
        "n_informative_source_photo_records":n,
        "source_cultivated_winter":count["a"],
        "source_cultivated_other_months":count["b"],
        "source_uncultivated_winter":count["c"],
        "source_uncultivated_other_months":count["d"],
        "native_host_taxa_in_informative_strata":len(taxa) if "host" in keys else None,
        "state_units":len(states),
        "largest_state_record_fraction":max(by_state.values())/n if n else None,
        "source_information_gate":("PASS_PREDECLARED_SOURCE_SUPPORT"
          if strict_gate else "FAIL_PREDECLARED_SOURCE_SUPPORT"),
        "mh_recording_composition_odds_ratio":odds,
        "bootstrap":group_bootstrap(strata,iterations=iterations) if require_state_bootstrap
                     else {"valid_draws":0,"interval":None},
        "not_plant_or_butterfly_abundance":True,
    }


def analyze(rows, bootstrap_iterations=999):
    # Main native-host species x state/province x year conditional comparison.
    primary, q1 = table_strata(rows)
    strict_hex, qhex=table_strata(rows, keys=("host","hex","state","year"))
    relaxed, q2=table_strata(rows,label_mode="expanded")
    highconf, q3=table_strata(rows,high_confidence=True)
    shoulder, q4=table_strata(rows,shoulder=True)
    alien, q5=table_strata(rows,origin="non-native")
    origin, q6=table_strata(rows, comparison="origin",keys=("state","year"))
    return {
        "schema":OUTPUT_SCHEMA,
        "raw_source_is_original_thinned_12942_photo_records":True,
        "primary_native_host_cultivation_strict":summary(
            primary,keys=("host","state","year"),iterations=bootstrap_iterations),
        "sensitivity_within_100km_host_year":summary(
            strict_hex,keys=("host","hex","state","year"),iterations=bootstrap_iterations),
        "sensitivity_likely_unlikely_included":summary(
            relaxed,keys=("host","state","year"),iterations=bootstrap_iterations),
        "sensitivity_high_confidence_only":summary(
            highconf,keys=("host","state","year"),iterations=bootstrap_iterations),
        "sensitivity_equal_length_shoulder_vs_winter":summary(
            shoulder,keys=("host","state","year"),iterations=bootstrap_iterations),
        "secondary_nonnative_host_cultivation":summary(
            alien,keys=("host","state","year"),iterations=bootstrap_iterations),
        "secondary_uncultivated_native_vs_nonnative":summary(
            origin,keys=("state","year"),iterations=bootstrap_iterations),
        "qc_filtered_source": {
            "main":q1,"hex":qhex,"expanded":q2,"high_confidence":q3,
            "shoulder":q4,"nonnative":q5,"origin":q6
        },
    }


def run(outfile, bootstrap_iterations=999):
    with TemporaryDirectory(prefix="chocho-public-monarch-") as d:
        file=Path(d)/"source.zip"
        source_download(SHARE_URL,file)
        sha=hashlib.sha256(file.read_bytes()).hexdigest()
        if sha!=ORIGINAL_ZIP_SHA256:
            raise ValueError("Public input source changed: stop rather than substitute")
        with ZipFile(file) as archive:
            rows, quality=read_source(archive)
    result=analyze(rows,bootstrap_iterations)
    result["source"]={
        "original_zip_sha256":sha,
        "public_citation":"Freedman 2026 doi:10.64898/2026.04.06.716508",
        "author_thinned_rows":SOURCE_ROWS,
        "field_eligible_rows_after_rearing_origin_cultivation_geography":len(rows),
        "source_quality":quality,
    }
    primary=result["primary_native_host_cultivation_strict"]
    result["decision"]=(
        "SOURCE_SUPPORTED_RECORD_COMPOSITION_NOT_ECOLOGICAL_CAUSAL"
        if primary["source_information_gate"]=="PASS_PREDECLARED_SOURCE_SUPPORT"
           and primary["mh_recording_composition_odds_ratio"] is not None
        else "INSUFFICIENT_MATCHED_NATIVE_HOST_CULTIVATION_RECORDING_EVIDENCE"
    )
    result["interpretation_red_lines"]=[
        "Freedman (2026) already reported nonnative milkweed winter host geography; not a chocho first discovery",
        "All data are monarch larvae already detected, not a controlled numerator over searched host stems",
        "Host cultivation is an observer/image annotation, NOT leaf phenology or experimentally assigned irrigation",
        "State×host×year matching does not remove seasonal observer effort or other unmeasured plant conditions",
        "Repeated photos, habitat, and geography are correlated; bootstrap is descriptive",
        "Do not infer larvae completing development, population increase, niche evolution or benefit of alien plants",
        "No individual coordinates, observer names, photos, or observation bodies are exported",
    ]
    outfile.parent.mkdir(parents=True,exist_ok=True)
    outfile.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"decision":result["decision"],
                      "primary":{k:primary[k] for k in (
                         "candidate_informative_strata","n_informative_source_photo_records",
                         "native_host_taxa_in_informative_strata","state_units",
                         "source_cultivated_winter","source_uncultivated_winter",
                         "source_information_gate","mh_recording_composition_odds_ratio")},
                      "other":{key:result[key]["mh_recording_composition_odds_ratio"]
                               for key in ("sensitivity_within_100km_host_year",
                                           "sensitivity_likely_unlikely_included",
                                           "secondary_nonnative_host_cultivation")},
                      },sort_keys=True))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--bootstrap-iterations",type=int,default=999)
    args=p.parse_args()
    if not 0<=args.bootstrap_iterations<=10000:
        raise ValueError("Invalid bootstrap bound")
    run(args.output,args.bootstrap_iterations)

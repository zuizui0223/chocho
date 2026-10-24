#!/usr/bin/env python3
"""Outcome-blind pilot for strict georeferenced larval feeding evidence in GloBI.

This intentionally reports only aggregated taxon/source-field counts. No raw
locations, dates, observation IDs or source-row copies are persisted.
See provenance/reviewer_defenses/protocols/GLOBI_LARVAL_HOST_SOURCE_FEASIBILITY_20261024.md.
"""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
EXPECTED_BLOB = "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"
BASE = "https://api.globalbioticinteractions.org/interaction.csv"
LIMIT = 1024
FIELDS = (
    "source_taxon_name", "interaction_type", "target_taxon_name",
    "target_taxon_path", "source_specimen_life_stage", "event_date",
    "latitude", "longitude", "study_source_id", "study_citation",
    "study_source_citation",
)
IMMATURE = re.compile(r"larva|larval|caterpillar|juvenile|instar", re.I)
SUSPECT = re.compile(r"(^|[/_.-])(nhm|hosts|naturalhistorymuseum)([/_.-]|$)", re.I)
IS_PLANT = re.compile(r"(?:^|\W)(plantae|viridiplantae)(?:\W|$)", re.I)
CONTROLS = ("Danaus plexippus", "Pieris brassicae")


def species():
    data = INPUT.read_bytes()
    git_blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()
    if git_blob != EXPECTED_BLOB:
        raise ValueError("frozen focal species source changed; abort")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    names = [r["species"] for r in reader]
    if len(names) != 239 or len(set(names)) != 239:
        raise ValueError("focal source must contain 239 unique names")
    def score(n):
        return hashlib.sha256(("20261024|" + n).encode()).digest()
    sample = sorted(names, key=score)[:20]
    out = sorted(set(sample).union(CONTROLS))
    return sample, out, git_blob


def date_known(text):
    match = re.match(r"^(\d{4})(?:-(\d{2})-(\d{2}))?", str(text or ""))
    if not match:
        return False, False
    year = int(match.group(1))
    if not 1800 <= year <= 2026:
        return False, False
    full = False
    if match.group(2):
        try:
            dt.date(year, int(match.group(2)), int(match.group(3)))
            full = True
        except (ValueError, TypeError):
            pass
    return True, full


def georeferenced(row):
    try:
        lat, lon = float(row.get("latitude") or "nan"), float(row.get("longitude") or "nan")
    except ValueError:
        return False
    return (math.isfinite(lat) and math.isfinite(lon) and
            -90 <= lat <= 90 and -180 <= lon <= 180 and (lat, lon) != (0, 0))


def trophic(name):
    s = str(name or "").lower().replace("_", "").replace(" ", "")
    # Do not silently treat pollinates/visits/hostOf as larval feeding.
    return s in {"eats", "feedson", "eatenby"} and s != "eatenby"


def source_independent(source_id):
    name = str(source_id or "").strip()
    return bool(name) and not SUSPECT.search(name)


def summarize_rows(rows, species_name):
    counts = collections.Counter()
    types = collections.Counter()
    sources = collections.Counter()
    qualifying_sources = set()
    stage_terms = collections.Counter()
    for row in rows:
        counts["rows"] += 1
        if row.get("source_taxon_name", "").strip() != species_name:
            counts["other_source_taxon"] += 1
            continue
        counts["exact_focal"] += 1
        typ = row.get("interaction_type", "").strip()
        types[typ or "[missing]"] += 1
        plant = bool(IS_PLANT.search(row.get("target_taxon_path", "")))
        stage = row.get("source_specimen_life_stage", "").strip()
        immature = bool(IMMATURE.search(stage))
        location = georeferenced(row)
        yr, full = date_known(row.get("event_date", ""))
        sid = row.get("study_source_id", "").strip()
        independent = source_independent(sid)
        source_key = sid or "[missing]"
        sources[source_key] += 1
        if stage:
            stage_terms[stage[:80]] += 1
        for k, ok in (
            ("plant_target", plant),
            ("explicit_trophic", trophic(typ)),
            ("immature_stage", immature),
            ("coordinates", location),
            ("dated_year", yr),
            ("dated_full", full),
            ("source_not_HOSTS_marked", independent),
        ):
            if ok:
                counts[k] += 1
        if plant and trophic(typ):
            counts["plant_trophic"] += 1
            if immature:
                counts["plant_trophic_immature"] += 1
                if location:
                    counts["plant_trophic_immature_georeferenced"] += 1
                    if independent:
                        counts["strict_candidate"] += 1
                        qualifying_sources.add(sid)
                        if yr:
                            counts["strict_candidate_dated"] += 1
    return {
        "counts": dict(sorted(counts.items())),
        "interaction_type_counts": dict(types.most_common(20)),
        "observed_stage_terms": dict(stage_terms.most_common(10)),
        "distinct_study_source_ids": len(sources),
        "study_source_id_top10": dict(sources.most_common(10)),
        "strict_qualifying_source_ids": sorted(qualifying_sources),
    }


def fetch(name):
    params = {
        "sourceTaxon": name,
        "includeObservations": "true",
        "limit": str(LIMIT),
        "type": "csv",
        "fields": ",".join(FIELDS),
    }
    url = BASE + "?" + urlencode(params)
    result = {"species": name, "url": url, "status": "NOT_RETRIEVED"}
    req = Request(
        url,
        headers={
            "Accept": "text/csv",
            "User-Agent": "chocho-research-source-feasibility/0.1 (+https://github.com/zuizui0223/chocho)",
        },
    )
    try:
        with urlopen(req, timeout=25) as response:
            status = response.status
            data = response.read(4_000_001)
            result["http_status"] = status
            result["content_type"] = response.headers.get("Content-Type", "")
        if status != 200:
            result["status"] = "HTTP_NOT_OK"
            return result
        if len(data) > 4_000_000:
            result["status"] = "RESPONSE_TOO_LARGE"
            return result
        result["source_response_sha256"] = hashlib.sha256(data).hexdigest()
        result["response_bytes"] = len(data)
        decoded = data.decode("utf-8-sig", errors="strict")
        if not decoded.strip():
            result["status"] = "EMPTY_BODY_UNVERIFIED"
            return result
        reader = csv.DictReader(io.StringIO(decoded))
        head = reader.fieldnames or []
        result["returned_header"] = head
        essentials = {"source_taxon_name", "target_taxon_name", "interaction_type", "study_source_id"}
        if not essentials.issubset(head):
            result["status"] = "WRONG_SOURCE_SCHEMA"
            return result
        records = list(reader)
        result["rows_returned"] = len(records)
        if len(records) >= LIMIT:
            result["status"] = "PAGINATION_REQUIRED"
            # keep bounded descriptive counts but do not call complete
        else:
            result["status"] = "COMPLETE"
        result.update(summarize_rows(records, name))
        return result
    except HTTPError as e:
        result["status"] = "HTTP_ERROR"
        result["http_status"] = e.code
    except (URLError, TimeoutError, OSError, UnicodeError, csv.Error) as e:
        result["status"] = "SOURCE_ERROR"
        result["error_type"] = type(e).__name__
    return result


def run(pilot_only=False):
    sample, names, sha = species()
    if pilot_only:
        # Transport check (not readiness gate), choose 2 hashed focal names
        names = sorted(set(sample[:2]).union(CONTROLS))
    out = []
    for i, name in enumerate(names):
        out.append(fetch(name))
        if i < len(names)-1:
            time.sleep(0.25)
    selected = [z for z in out if z["species"] in sample]
    complete = (len(selected) == 20 and all(z["status"] == "COMPLETE" for z in selected))
    taxa_with_strict = [z for z in selected if z.get("counts", {}).get("strict_candidate", 0) > 0]
    distinct_qualifying_sources = set()
    for s in taxa_with_strict:
        distinct_qualifying_sources.update(s.get("strict_qualifying_source_ids", []))
    if pilot_only:
        status = "TRANSPORT_SMOKE_ONLY_NOT_GATE"
    elif not complete:
        status = "SOURCE_GATE_INCOMPLETE_NOT_NEGATIVE_EVIDENCE"
    elif len(taxa_with_strict) >= 10 and len(distinct_qualifying_sources) >= 5:
        status = "SOURCE_GATE_FEASIBLE_NEEDS_VERSIONED_REPLICATION"
    else:
        status = "SOURCE_GATE_BELOW_MINIMUM_NO_REALIZED_MACRO_CLAIM"
    return {
        "schema": "chocho_globi_immature_georeferenced_source_pilot_v01",
        "status": status,
        "date_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source": "GloBI live API (unstable feasibility probe, NOT pinned archival release)",
        "source_versioned_archive": "https://doi.org/10.5281/zenodo.3950589",
        "focal_git_blob": sha,
        "focal_n": 239,
        "random_focal_selected": sample,
        "positive_controls": list(CONTROLS),
        "source_requests": len(out),
        "response_blind_selected_queries": len(selected),
        "complete_random_focal": complete,
        "random_focal_with_strict_candidate": len(taxa_with_strict),
        "strict_candidate_distinct_source_ids": len(distinct_qualifying_sources),
        "next_decision": status,
        "important_limits": [
            "Source IDs are not a fully verified independence check; may reuse HOSTS citations",
            "Larval georeferenced feeding is not successful adult recruitment",
            "Georeferenced interactions do not establish nonnative plant at encounter",
            "Live API data unstable; full study must pin a versioned GloBI archive",
            "Errors, missingness and truncation are not ecological absence",
        ],
        "taxa": out,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    r = run(pilot_only=args.smoke)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("SOURCE_STATUS", r["status"])
    print("GLOBI_PROBE_SUMMARY", json.dumps({
        "requested": r["source_requests"],
        "status_counts": dict(collections.Counter(x["status"] for x in r["taxa"])),
        "strict_positive": r["random_focal_with_strict_candidate"],
        "distinct_source_ids": r["strict_candidate_distinct_source_ids"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fixed source query pagination v0.2 after truncated first GloBI pilot.

Reuses v0.1 original panel, stage/location/provenance definitions, and gate.
No independent biological effect is estimated.
"""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import sys
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import audit_globi_larval_host_macro_source_20261024 as old

MAX_FOCAL_PAGES = 12
REQUEST_BYTES_LIMIT = 4_000_000
RESULT_NAME = "provenance/reviewer_defenses/results/globi_larval_host_feasibility_paginated_20261024.json"


def fetch_pages(name: str, max_pages: int) -> dict:
    result = {
        "species": name,
        "status": "NOT_RETRIEVED",
        "paging_cap": max_pages,
        "pagination": [],
    }
    rows = []
    seen_hashes = set()
    first_header = None
    for page in range(max_pages):
        skip = page * old.LIMIT
        params = {
            "sourceTaxon": name,
            "includeObservations": "true",
            "limit": str(old.LIMIT),
            "skip": str(skip),
            "type": "csv",
            "fields": ",".join(old.FIELDS),
        }
        url = old.BASE + "?" + urlencode(params)
        req = Request(url, headers={
            "Accept": "text/csv",
            "User-Agent": "chocho-research-source-feasibility/0.2 (+https://github.com/zuizui0223/chocho)",
        })
        try:
            with urlopen(req, timeout=25) as response:
                http_status = response.status
                data = response.read(REQUEST_BYTES_LIMIT + 1)
                content_type = response.headers.get("Content-Type", "")
            info = {"skip": skip, "http_status": http_status, "content_type": content_type,
                    "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
            result["pagination"].append(info)
            if http_status != 200:
                result["status"] = "HTTP_NOT_OK_INCOMPLETE"
                break
            if len(data) > REQUEST_BYTES_LIMIT:
                result["status"] = "TOO_LARGE_INCOMPLETE"
                break
            if info["sha256"] in seen_hashes:
                result["status"] = "PAGINATION_DUPLICATE_RESPONSE_INCOMPLETE"
                break
            seen_hashes.add(info["sha256"])
            text = data.decode("utf-8-sig")
            if not text.strip():
                result["status"] = "EMPTY_RESPONSE_INCOMPLETE"
                break
            reader = csv.DictReader(io.StringIO(text))
            header = reader.fieldnames or []
            essential = {"source_taxon_name", "interaction_type", "target_taxon_name", "study_source_id"}
            if not essential.issubset(header):
                result["status"] = "WRONG_SOURCE_SCHEMA_INCOMPLETE"
                break
            if first_header is not None and header != first_header:
                result["status"] = "HEADER_DRIFT_INCOMPLETE"
                break
            first_header = header
            page_rows = list(reader)
            info["rows"] = len(page_rows)
            rows.extend(page_rows)
            if len(page_rows) < old.LIMIT:
                result["status"] = "COMPLETE"
                break
        except HTTPError as e:
            result["status"] = "HTTP_ERROR_INCOMPLETE"
            result["http_status"] = e.code
            break
        except (URLError, TimeoutError, OSError, UnicodeError, csv.Error) as e:
            result["status"] = "SOURCE_ERROR_INCOMPLETE"
            result["error_type"] = type(e).__name__
            break
        time.sleep(0.25)
    else:
        result["status"] = "PAGINATION_CAP_INCOMPLETE"
    result["source_response_sha256"] = hashlib.sha256(
        "|".join(p["sha256"] for p in result["pagination"]).encode()).hexdigest()
    result["response_bytes"] = sum(p.get("bytes", 0) for p in result["pagination"])
    result["rows_returned"] = len(rows)
    if rows:
        result.update(old.summarize_rows(rows, name))
    return result


def run():
    focal, names, blob = old.species()
    selected = []
    out = []
    for i, name in enumerate(names):
        page_cap = MAX_FOCAL_PAGES if name in focal else 1
        z = fetch_pages(name, page_cap)
        out.append(z)
        if name in focal:
            selected.append(z)
        if i < len(names)-1:
            time.sleep(0.25)
    complete = (len(selected) == 20 and all(z["status"] == "COMPLETE" for z in selected))
    taxa_strict = [z for z in selected if z.get("counts", {}).get("strict_candidate", 0) > 0]
    source_ids = sorted(set(s for z in taxa_strict for s in z.get("strict_qualifying_source_ids", [])))
    if not complete:
        status = "SOURCE_GATE_INCOMPLETE_AFTER_FIXED_PAGINATION"
    elif len(taxa_strict) >= 10 and len(source_ids) >= 5:
        status = "SOURCE_GATE_FEASIBLE_NEEDS_VERSIONED_ARCHIVE_AND_SOURCE_INDEPENDENCE"
    else:
        status = "SOURCE_GATE_BELOW_MINIMUM_NO_REALIZED_MACRO_CLAIM"
    return {
        "schema": "chocho_globi_immature_georeferenced_source_pilot_paginated_v02",
        "status": status,
        "date_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source": "Live GloBI API (not a version-pinned archive)",
        "source_versioned_archive": "https://doi.org/10.5281/zenodo.3950589",
        "base_v01_receipt": "provenance/reviewer_defenses/results/globi_larval_host_feasibility_20261024.json",
        "protocol": "provenance/reviewer_defenses/protocols/GLOBI_LARVAL_HOST_SOURCE_FEASIBILITY_20261024.md",
        "focal_git_blob": blob,
        "focal_n": 239,
        "random_focal_selected": focal,
        "positive_controls": list(old.CONTROLS),
        "source_requests": len(out),
        "complete_random_focal": complete,
        "random_focal_with_strict_candidate": len(taxa_strict),
        "strict_candidate_distinct_source_ids": len(source_ids),
        "strict_candidate_source_ids": source_ids,
        "all_strict_candidate_rows": sum(z.get("counts", {}).get("strict_candidate", 0) for z in selected),
        "sampling": "20 response-blind SHA256-min focal taxa; positive controls excluded",
        "paging_limit": {"rows_per_page": old.LIMIT, "max_focal_pages": MAX_FOCAL_PAGES,
                         "controls_max_pages": 1},
        "important_limits": [
            "GloBI data may reuse source HOSTS records even with different study_source_id",
            "Live pages not a transactional version-stable snapshot",
            "Source rows are not independent individuals or confirmed larval feeding events",
            "Stage-labelled geographic trophic record not adult emergence, colonization or establishment",
            "No introduced or native status of observed individual plant known",
            "Any failed, capped, or missing source query cannot establish ecological absence",
        ],
        "taxa": out,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print("GLOBI_PAGED", json.dumps({
        "status": result["status"],
        "focal_complete": sum(z["status"] == "COMPLETE" for z in result["taxa"]
                              if z["species"] in result["random_focal_selected"]),
        "focal_strict": result["random_focal_with_strict_candidate"],
        "source_ids": result["strict_candidate_distinct_source_ids"],
        "rows": result["all_strict_candidate_rows"],
        "source_statuses": dict(collections.Counter(z["status"] for z in result["taxa"])),
    }))


if __name__ == "__main__":
    main()

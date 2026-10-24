#!/usr/bin/env python3
"""Acquire a NEW country-grain full-index GBIF facet dataset.

This does not repair, replace, or pool with the original bounded WGSRPD3
occurrence sample. Every request has limit=0: no occurrence records are
downloaded. Country totals describe indexed source records, never survey
effort, detection probability, biological absence, or host acceptance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import tempfile
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.parse import urlencode


API = "https://api.gbif.org/v1/occurrence/search"
SCHEMA = "chocho_GBIF_country_full_index_facets_2010_2025_20261020_v01"
ORIGINAL_SCHEMA = "chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01"
ORIGINAL_STATUSES = {"COMPLETE": 135, "PARTIAL_FETCH_NOT_ELIGIBLE": 96, "TAXON_UNRESOLVED": 8}
WINDOWS = (("2010_2017", "2010,2017"), ("2018_2025", "2018,2025"))
BASIS_SCOPES = {
    "all_accepted_bases": (
        "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
        "MACHINE_OBSERVATION", "OCCURRENCE",
    ),
    "core": ("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"),
}
FACET_LIMIT = 300
USER_AGENT = "chocho-country-facet-macroecology/0.1 (https://github.com/zuizui0223/chocho)"
TRANSIENT_HTTP = {408, 425, 429, 500, 502, 503, 504}
TRANSIENT_CURL = {5, 6, 7, 18, 28, 35, 52, 55, 56}
QUERY_IDENTITY = (
    "query_id", "species", "source_status", "taxon_key", "basis_scope",
    "window", "year_filter", "url", "query_parameters", "original_manifest_sha256",
)


class FacetValidationError(ValueError):
    """A response cannot establish complete country facet acquisition."""


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def write_atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_bytes(data)
    temp.replace(path)


def write_json(path, value):
    write_atomic(path, (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode())


def load_original_manifest(path):
    path = Path(path)
    raw = path.read_bytes()
    data = json.loads(raw)
    if data.get("schema") != ORIGINAL_SCHEMA:
        raise ValueError("Original source schema does not match the original 239-species manifest")
    rows = data.get("species_receipts")
    if not isinstance(rows, list) or len(rows) != 239 or data.get("species_requested") != 239:
        raise ValueError("Original manifest must contain exactly 239 source species")
    status_counts = dict(Counter(row.get("status") for row in rows))
    if status_counts != ORIGINAL_STATUSES or data.get("species_statuses") != ORIGINAL_STATUSES:
        raise ValueError("Original receipt status counts must be exactly 135 COMPLETE / 96 PARTIAL / 8 UNRESOLVED")
    names = set()
    for row in rows:
        name = row.get("species")
        if not isinstance(name, str) or not name.strip() or name in names:
            raise ValueError(f"Invalid or duplicate original species: {name!r}")
        names.add(name)
        if row["status"] != "TAXON_UNRESOLVED":
            key = (row.get("taxon") or {}).get("key")
            if type(key) is not int or key <= 0:
                raise ValueError(f"Invalid original taxon key for {name}")
    data["original_manifest_sha256"] = sha256(raw)
    data["original_manifest_path"] = str(path.resolve())
    return data


def build_queries(source):
    queries = []
    for row in source["species_receipts"]:
        if row["status"] == "TAXON_UNRESOLVED":
            continue
        for scope, bases in BASIS_SCOPES.items():
            for window, years in WINDOWS:
                params = [
                    ["taxonKey", str(row["taxon"]["key"])],
                    ["hasCoordinate", "true"], ["hasGeospatialIssue", "false"],
                    ["occurrenceStatus", "PRESENT"], ["year", years],
                    ["limit", "0"], ["facet", "country"], ["facetLimit", str(FACET_LIMIT)],
                ] + [["basisOfRecord", basis] for basis in bases]
                query = {
                    "species": row["species"], "source_status": row["status"],
                    "taxon_key": row["taxon"]["key"], "basis_scope": scope,
                    "window": window, "year_filter": years,
                    "query_parameters": params, "url": API + "?" + urlencode([tuple(pair) for pair in params]),
                    "original_manifest_sha256": source["original_manifest_sha256"],
                }
                query["query_id"] = sha256(canonical_bytes(query))[:32]
                queries.append(query)
    return queries


def nonnegative_int(value, label):
    if type(value) is not int or value < 0:
        raise FacetValidationError(f"{label} must be a nonnegative integer")
    return value


def parse_country_facet(raw):
    try:
        data = json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise FacetValidationError(f"Invalid or truncated JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise FacetValidationError("JSON response must be an object")
    total = nonnegative_int(data.get("count"), "GBIF count")
    if data.get("limit") != 0 or data.get("results") != []:
        raise FacetValidationError("Expected limit=0 and no occurrence results")
    facets = data.get("facets")
    if not isinstance(facets, list):
        raise FacetValidationError("Missing GBIF facets list")
    if not facets and total == 0:
        facets = [{"field": "COUNTRY", "counts": []}]
    if len(facets) != 1 or not isinstance(facets[0], dict) or facets[0].get("field") != "COUNTRY":
        raise FacetValidationError("Expected exactly one COUNTRY facet")
    entries = facets[0].get("counts")
    if not isinstance(entries, list):
        raise FacetValidationError("COUNTRY facet counts must be a list")
    if len(entries) >= FACET_LIMIT:
        raise FacetValidationError("COUNTRY facet reached facetLimit; possible truncation is not complete")
    counts = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise FacetValidationError("Malformed COUNTRY facet entry")
        country = entry.get("name")
        if not isinstance(country, str) or re.fullmatch(r"[A-Z]{2}", country) is None:
            raise FacetValidationError(f"Country must be a two-letter uppercase code: {country!r}")
        if country in counts:
            raise FacetValidationError(f"Duplicate COUNTRY facet code: {country}")
        counts[country] = nonnegative_int(entry.get("count"), f"Country {country} count")
    mapped = sum(counts.values())
    if mapped > total:
        raise FacetValidationError("COUNTRY facet count sum exceeds GBIF reported count")
    return {
        "reported_count": total, "country_count": len(counts),
        "country_counts": dict(sorted(counts.items())), "country_records": mapped,
        "missing_country_count": total - mapped,
    }


def save_complete_query(output_dir, query, raw, *, attempts):
    summary = parse_country_facet(raw)
    relative = f"raw_facets/{query['query_id']}.json"
    receipt = {
        **query, "schema": SCHEMA, "status": "COMPLETE", "complete": True,
        "completed_utc": utc_now(), "attempts": attempts,
        "body_sha256": sha256(raw), "body_bytes": len(raw), "raw_path": relative,
        **summary,
    }
    write_atomic(Path(output_dir) / relative, raw)
    write_json(Path(output_dir) / "checkpoints" / f"{query['query_id']}.json", receipt)
    return receipt


def read_complete_checkpoint(output_dir, query):
    root = Path(output_dir)
    checkpoint = root / "checkpoints" / f"{query['query_id']}.json"
    try:
        receipt = json.loads(checkpoint.read_bytes())
        if (not isinstance(receipt, dict) or receipt.get("schema") != SCHEMA or receipt.get("status") != "COMPLETE"
                or receipt.get("complete") is not True
                or any(receipt.get(key) != query[key] for key in QUERY_IDENTITY)):
            return None
        expected_raw = f"raw_facets/{query['query_id']}.json"
        if receipt.get("raw_path") != expected_raw:
            return None
        raw = (root / expected_raw).read_bytes()
        if receipt.get("body_sha256") != sha256(raw) or receipt.get("body_bytes") != len(raw):
            return None
        summary = parse_country_facet(raw)
        if any(receipt.get(key) != value for key, value in summary.items()):
            return None
        return receipt
    except (OSError, ValueError, TypeError, KeyError):
        return None


class GlobalRateLimiter:
    """Serialize starts across all workers, including retries."""

    def __init__(self, requests_per_second=4.0):
        if not 0 < requests_per_second <= 4:
            raise ValueError("Global request rate must be greater than zero and at most 4/s")
        self.interval = 1.0 / requests_per_second
        self.lock = threading.Lock()
        self.next_request = 0.0

    def wait(self):
        while True:
            with self.lock:
                now = time.monotonic()
                delay = self.next_request - now
                if delay <= 0:
                    self.next_request = now + self.interval
                    return
            # A response on another worker can extend the cooldown while we
            # sleep. Recheck it under lock before authorizing this request.
            time.sleep(delay)

    def cooldown(self, seconds):
        """Honor a server delay across every worker, not just its caller."""
        with self.lock:
            self.next_request = max(self.next_request, time.monotonic() + seconds)


def parse_retry_after(raw_headers, *, now=None):
    """Read final HTTP response Retry-After, in seconds or HTTP-date form."""
    blocks = re.split(rb"\r?\n\r?\n", raw_headers.strip())
    matches = re.findall(rb"(?im)^retry-after:\s*([^\r\n]+)", blocks[-1])
    if not matches:
        return None
    value = matches[-1].decode("ascii", errors="replace").strip()
    if value.isdigit():
        return float(value)
    try:
        date = parsedate_to_datetime(value)
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        return max(0.0, (date - (now or datetime.now(timezone.utc))).total_seconds())
    except (ValueError, TypeError, OverflowError):
        return None


def fetch_query(output_dir, query, limiter, *, timeout=30.0, tries=4):
    attempts = []
    checkpoint = Path(output_dir) / "checkpoints" / f"{query['query_id']}.json"
    try:
        prior = json.loads(checkpoint.read_bytes())
        if (isinstance(prior, dict) and prior.get("status") == "FAILED"
                and all(prior.get(key) == query[key] for key in QUERY_IDENTITY)):
            attempts = list(prior.get("attempts", []))
            archive_name = f"{query['query_id']}.{sha256(canonical_bytes(prior))[:16]}.json"
            write_json(Path(output_dir) / "failed_query_receipts" / archive_name, prior)
    except (OSError, ValueError, TypeError):
        pass
    final_error = "No response"
    for attempt_number in range(1, tries + 1):
        limiter.wait()
        started = time.monotonic()
        attempt = {"attempt": len(attempts) + 1, "attempt_in_current_run": attempt_number,
                   "requested_utc": utc_now(),
                   "global_max_requests_per_second": 1.0 / limiter.interval}
        cmd = [
            "curl", "--location", "--silent", "--show-error", "--compressed",
            "--connect-timeout", "8", "--max-time", str(timeout),
            "--header", f"User-Agent: {USER_AGENT}", "--header", "Accept: application/json",
            "--write-out", "\n%{http_code}", query["url"],
        ]
        raw = b""
        retryable = False
        try:
            with tempfile.TemporaryDirectory(prefix="chocho-gbif-headers-") as header_dir:
                header_path = Path(header_dir) / "response.headers"
                cmd[1:1] = ["--dump-header", str(header_path)]
                call = subprocess.run(cmd, capture_output=True, timeout=timeout + 3, check=False)
                headers = header_path.read_bytes() if header_path.exists() else b""
            retry_after = parse_retry_after(headers)
            body, separator, code = call.stdout.rpartition(b"\n")
            status = int(code) if separator and code.isdigit() else 0
            raw = body if separator else call.stdout
            attempt.update(http_status=status, curl_returncode=call.returncode,
                           retry_after_seconds=retry_after)
            if call.returncode == 0 and status == 200:
                try:
                    parse_country_facet(raw)
                except FacetValidationError as exc:
                    final_error = f"HTTP 200 facet validation failed: {exc}"
                else:
                    attempt.update(
                        completed_utc=utc_now(), duration_seconds=round(time.monotonic() - started, 4),
                        body_bytes=len(raw), body_sha256=sha256(raw),
                    )
                    attempts.append(attempt)
                    return save_complete_query(output_dir, query, raw, attempts=attempts)
            else:
                final_error = f"HTTP {status}; curl exit {call.returncode}: {call.stderr.decode(errors='replace')[:400]}"
                retryable = (status in TRANSIENT_HTTP or
                             (status in (0, 200) and call.returncode in TRANSIENT_CURL))
                if retryable and (status == 429 or retry_after is not None):
                    cooldown = retry_after if retry_after is not None else 5.0
                    limiter.cooldown(cooldown)
                    attempt["shared_cooldown_seconds"] = cooldown
        except subprocess.TimeoutExpired:
            final_error = "Whole-request subprocess timeout"
            attempt.update(http_status=0, curl_returncode=None)
            retryable = True
        except OSError as exc:
            final_error = f"Could not execute curl: {exc}"
            attempt.update(http_status=0, curl_returncode=None)
        attempt.update(
            completed_utc=utc_now(), duration_seconds=round(time.monotonic() - started, 4),
            body_bytes=len(raw), body_sha256=sha256(raw), error=final_error,
            transient_retry_eligible=retryable,
        )
        if raw:
            relative = f"failed_responses/{query['query_id']}.{time.time_ns()}.body"
            write_atomic(Path(output_dir) / relative, raw)
            attempt["raw_failure_path"] = relative
        attempts.append(attempt)
        if not retryable or attempt_number >= tries:
            break
        time.sleep(min(20.0, 2.0 ** (attempt_number - 1)))
    receipt = {
        **query, "schema": SCHEMA, "status": "FAILED", "complete": False,
        "completed_utc": utc_now(), "attempts": attempts, "error": final_error,
        "missingness_meaning": "SOURCE_QUERY_FAILED_NOT_ABSENCE",
    }
    write_json(Path(output_dir) / "checkpoints" / f"{query['query_id']}.json", receipt)
    archive_name = f"{query['query_id']}.{sha256(canonical_bytes(receipt))[:16]}.json"
    write_json(Path(output_dir) / "failed_query_receipts" / archive_name, receipt)
    return receipt


def accepted_key_collisions(source):
    names = defaultdict(list)
    for row in source["species_receipts"]:
        if row["status"] != "TAXON_UNRESOLVED":
            names[row["taxon"]["key"]].append(row["species"])
    return [{"taxon_key": key, "species": sorted(species)}
            for key, species in sorted(names.items()) if len(species) > 1]


def assemble_outputs(source, queries, supplied_receipts):
    received = {}
    expected = {q["query_id"]: q for q in queries}
    for receipt in supplied_receipts:
        query_id = receipt.get("query_id")
        if query_id not in expected or query_id in received:
            raise ValueError(f"Unexpected or duplicate query receipt: {query_id}")
        if any(receipt.get(key) != expected[query_id][key] for key in QUERY_IDENTITY):
            raise ValueError(f"Query identity mismatch in receipt: {query_id}")
        if receipt.get("status") not in ("COMPLETE", "FAILED", "MISSING"):
            raise ValueError(f"Invalid query status: {query_id}")
        if (receipt.get("status") == "COMPLETE") != (receipt.get("complete") is True):
            raise ValueError(f"Contradictory query completeness: {query_id}")
        received[query_id] = receipt
    query_receipts = [received.get(q["query_id"], {
        **q, "status": "MISSING", "complete": False,
        "missingness_meaning": "SOURCE_QUERY_NOT_FETCHED_NOT_ABSENCE",
    }) for q in queries]
    lookup = {(q["species"], q["basis_scope"], q["window"]): q for q in query_receipts}
    species_receipts, csv_rows, eligible = [], [], set()
    for original in source["species_receipts"]:
        item = {
            "species": original["species"], "source_status": original["status"],
            "taxon": original.get("taxon"), "basis_scopes": {},
        }
        if original["status"] == "TAXON_UNRESOLVED":
            item.update(status="TAXON_UNRESOLVED", taxon_key=None)
            item["basis_scopes"] = {scope: {"status": "TAXON_UNRESOLVED", "windows": {}}
                                    for scope in BASIS_SCOPES}
        else:
            item["taxon_key"] = original["taxon"]["key"]
            for scope in BASIS_SCOPES:
                windows = {window: lookup[(original["species"], scope, window)] for window, _ in WINDOWS}
                complete = all(q["status"] == "COMPLETE" and q["complete"] for q in windows.values())
                item["basis_scopes"][scope] = {
                    "status": "COMPLETE" if complete else "PARTIAL_FETCH_NOT_ELIGIBLE",
                    "windows": windows,
                }
                if complete:
                    eligible.add((original["species"], scope))
                    for window, query in windows.items():
                        for country, count in sorted(query["country_counts"].items()):
                            csv_rows.append({
                                "species": original["species"], "window": window, "country": country,
                                "basis_scope": scope, "record_count": count,
                            })
            item["status"] = (
                "COMPLETE" if all(v["status"] == "COMPLETE" for v in item["basis_scopes"].values())
                else "PARTIAL_FETCH_NOT_ELIGIBLE"
            )
        species_receipts.append(item)
    totals, discrepancies = {}, []
    for scope in BASIS_SCOPES:
        totals[scope] = {}
        for window, _ in WINDOWS:
            records = [q for q in query_receipts if q["basis_scope"] == scope and q["window"] == window
                       and (q["species"], scope) in eligible]
            by_key = defaultdict(list)
            for query in records:
                by_key[query["taxon_key"]].append(query)
            representatives = [sorted(group, key=lambda q: q["species"])[0] for group in by_key.values()]
            for key, group in by_key.items():
                if len({canonical_bytes([q["reported_count"], q["country_counts"]]) for q in group}) > 1:
                    discrepancies.append({
                        "taxon_key": key, "basis_scope": scope, "window": window,
                        "species_counts": {q["species"]: q["reported_count"] for q in group},
                        "reason": "Shared accepted-key responses differ across retrieval times",
                    })
            totals[scope][window] = {}
            for grain, subset in (("per_source_name", records), ("per_unique_accepted_key", representatives)):
                totals[scope][window][grain] = {
                    "source_names_count": len(subset),
                    "reported_count": sum(q["reported_count"] for q in subset),
                    "country_records": sum(q["country_records"] for q in subset),
                    "missing_country_count": sum(q["missing_country_count"] for q in subset),
                }
                if grain == "per_unique_accepted_key":
                    totals[scope][window][grain]["source_name_representatives"] = sorted(q["species"] for q in subset)
    result = {
        "schema": SCHEMA,
        "status": "COMPLETE_RESOLVED_SOURCE_NAMES" if all(q["status"] == "COMPLETE" for q in query_receipts)
                  else "PARTIAL_SOURCE_QUERY_MISSINGNESS",
        "source_grain": "GBIF indexed country code, full search-index facets",
        "separate_dataset_not_recovery_of_original_bounded_sample": True,
        "original_manifest_path": source.get("original_manifest_path"),
        "original_manifest_sha256": source["original_manifest_sha256"],
        "original_species_statuses": source["species_statuses"],
        "species_requested": len(species_receipts), "resolved_source_names": len(queries) // 4,
        "resolved_unique_accepted_keys": len({q["taxon_key"] for q in queries}),
        "species_statuses": dict(sorted(Counter(r["status"] for r in species_receipts).items())),
        "species_statuses_by_basis": {
            scope: dict(sorted(Counter(r["basis_scopes"][scope]["status"] for r in species_receipts).items()))
            for scope in BASIS_SCOPES
        },
        "query_statuses": dict(sorted(Counter(q["status"] for q in query_receipts).items())),
        "expected_query_count": len(queries), "query_receipts": query_receipts,
        "species_receipts": species_receipts, "record_totals": totals,
        "record_total_scope": "Only source-name/basis pairs with both windows complete; bases overlap and must not be summed",
        "unique_key_representative_rule": "Lexicographically first eligible original source name for each accepted taxon key",
        "accepted_key_collisions": accepted_key_collisions(source),
        "accepted_key_count_discrepancies": discrepancies,
        "basis_scopes": BASIS_SCOPES, "windows": dict(WINDOWS),
        "filters": {"hasCoordinate": True, "hasGeospatialIssue": False, "occurrenceStatus": "PRESENT"},
        "occurrence_records_downloaded": 0, "facet_limit": FACET_LIMIT,
        "source_response_set_sha256": sha256(canonical_bytes(sorted(
            [q["query_id"], q["body_sha256"]] for q in query_receipts if q.get("body_sha256")
        ))),
        "limitations": [
            "Index snapshot queries are live, not a DOI occurrence download or fixed global database snapshot.",
            "Country facets are a new coarser geographic grain; do not pool with bounded WGSRPD3 samples.",
            "Counts are occurrence records, not independent individuals, effort, detection probability, or survey non-detections.",
            "Missing or failed country/window/species queries never imply biological absence.",
            "Reported total minus country facet sum is retained as missing-country records and is not mapped to a country.",
            "API geospatial filters differ from the original coordinate-level custom filters; no 100-km uncertainty cutoff is added.",
            "Country can reflect GBIF interpreted country metadata and is not a point-in-polygon remapping.",
            "All accepted bases contains core; these scopes are sensitivity alternatives, not additive counts.",
            "Both original aliases of a shared accepted key remain in source receipts and require downstream collision handling.",
        ],
    }
    return result, csv_rows


def write_outputs(output_dir, source, queries, receipts, *, started_utc, resumed):
    result, rows = assemble_outputs(source, queries, receipts)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=("species", "window", "country", "basis_scope", "record_count"))
    writer.writeheader()
    writer.writerows(rows)
    raw_csv = buffer.getvalue().encode()
    root = Path(output_dir)
    csv_path = root / "gbif_country_counts.csv"
    receipt_path = root / "gbif_country_facet_receipts.json"
    write_atomic(csv_path, raw_csv)
    result.update(
        started_utc=started_utc, completed_utc=utc_now(),
        acquisition_script_sha256=sha256(Path(__file__).read_bytes()),
        checkpoint_queries_reused=resumed,
        outputs={"country_counts_csv": str(csv_path.resolve()), "country_counts_sha256": sha256(raw_csv),
                 "country_count_rows": len(rows)},
    )
    write_json(receipt_path, result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 5))
    parser.add_argument("--requests-per-second", type=float, default=4.0)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--tries", type=int, default=4)
    parser.add_argument("--assemble-only", action="store_true", help="Revalidate saved complete checkpoints without network")
    args = parser.parse_args(argv)
    if args.timeout <= 0 or not 1 <= args.tries <= 6:
        parser.error("timeout must be positive and tries must be between 1 and 6")
    limiter = GlobalRateLimiter(args.requests_per_second)
    source = load_original_manifest(args.original_manifest)
    queries = build_queries(source)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    receipts, pending = [], []
    for query in queries:
        saved = read_complete_checkpoint(args.output_dir, query)
        (receipts if saved is not None else pending).append(saved if saved is not None else query)
    resumed = len(receipts)
    print(json.dumps({"event": "start", "queries": len(queries), "complete_checkpoints": resumed,
                      "remaining": len(pending), "workers": args.workers,
                      "global_max_requests_per_second": args.requests_per_second}), flush=True)
    if not args.assemble_only:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_queries = {
                executor.submit(fetch_query, args.output_dir, q, limiter, timeout=args.timeout, tries=args.tries): q
                for q in pending
            }
            for future in as_completed(future_queries):
                receipt = future.result()
                receipts.append(receipt)
                if len(receipts) % 25 == 0 or receipt["status"] != "COMPLETE" or len(receipts) == len(queries):
                    print(json.dumps({"event": "progress", "finished": len(receipts), "expected": len(queries),
                                      "query_statuses": dict(Counter(r["status"] for r in receipts)),
                                      "last_species": receipt["species"]}), flush=True)
    else:
        for query in pending:
            path = args.output_dir / "checkpoints" / f"{query['query_id']}.json"
            try:
                failed = json.loads(path.read_bytes())
                if (isinstance(failed, dict) and failed.get("status") == "FAILED"
                        and all(failed.get(key) == query[key] for key in QUERY_IDENTITY)):
                    receipts.append(failed)
            except (OSError, ValueError, TypeError):
                pass
    result = write_outputs(args.output_dir, source, queries, receipts, started_utc=started, resumed=resumed)
    print(json.dumps({"event": "complete", "status": result["status"],
                      "query_statuses": result["query_statuses"], "species_statuses": result["species_statuses"],
                      "outputs": result["outputs"], "source_response_set_sha256": result["source_response_set_sha256"]}), flush=True)
    return 0 if result["status"] == "COMPLETE_RESOLVED_SOURCE_NAMES" else 2


if __name__ == "__main__":
    raise SystemExit(main())

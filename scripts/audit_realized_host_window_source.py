#!/usr/bin/env python3
"""Audit independent Dryad source-byte access; never infer butterfly effects."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

URLS = {
    "plant_phenology": "https://datadryad.org/downloads/file_stream/79537",
    "sampling_areas": "https://datadryad.org/downloads/file_stream/79538",
}
MAX_BYTES = 8_000_000


def fetch_source(name: str, url: str) -> dict:
    item = {"name": name, "url": url, "status": "SOURCE_ACCESS_BLOCKED"}
    try:
        req = Request(url, headers={"User-Agent": "chocho-source-provenance-audit/1.0"})
        with urlopen(req, timeout=20) as response:
            raw = response.read(MAX_BYTES + 1)
            item["http_status"] = response.status
            item["mime_type"] = response.headers.get("Content-Type")
        if len(raw) > MAX_BYTES or not raw or raw.lstrip().lower().startswith((b"<!doctype html", b"<html")):
            item["reason"] = "missing, oversized, or HTML response"
            return item
        item["bytes"] = len(raw)
        item["sha256"] = hashlib.sha256(raw).hexdigest()
        item["line_count"] = len(raw.splitlines())
        try:
            lines = raw.decode("utf-8-sig").splitlines()
        except UnicodeDecodeError:
            item["reason"] = "unknown text encoding; inspect manually"
            return item
        if len(lines) < 2 or len(lines[0]) > 5000:
            item["reason"] = "missing tabular header or records"
            return item
        try:
            dialect = csv.Sniffer().sniff("\n".join(lines[:20]), delimiters="\t,;")
            columns = next(csv.reader([lines[0]], dialect=dialect))
        except csv.Error:
            item["reason"] = "delimiter/schema cannot be identified"
            return item
        if len(columns) < 2:
            item["reason"] = "not a multicolumn table"
            return item
        item["columns"] = [c.strip() for c in columns]
        item["delimiter"] = "TAB" if dialect.delimiter == "\t" else dialect.delimiter
        item["status"] = "BYTES_DOWNLOADED_SCHEMA_REQUIRES_MANUAL_AUDIT"
        item["warning"] = "Source bytes alone do not validate stage codes, zero-egg denominators, adult flight dates or any ecological effect."
        return item
    except HTTPError as err:
        item["http_status"] = err.code
        item["reason"] = f"HTTP {err.code}; no data inspected"
    except (URLError, TimeoutError, OSError) as err:
        item["reason"] = f"source access failed: {type(err).__name__}"
    return item


def audit(output: Path) -> dict:
    sources = [fetch_source(k, u) for k, u in URLS.items()]
    receipt = {
        "schema": "chocho_phenological_source_access_receipt_v01",
        "source_doi": "10.5061/dryad.jp328r6",
        "sources": sources,
        "source_bytes_complete": all(s["status"] == "BYTES_DOWNLOADED_SCHEMA_REQUIRES_MANUAL_AUDIT" for s in sources),
        "row_level_ecological_effect_estimated": False,
        "main_manuscript_changed": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    r = audit(args.receipt)
    print(json.dumps({"source_bytes_complete": r["source_bytes_complete"], "statuses": [s["status"] for s in r["sources"]]}))


if __name__ == "__main__":
    main()

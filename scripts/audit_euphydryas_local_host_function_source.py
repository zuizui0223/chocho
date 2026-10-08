#!/usr/bin/env python3
"""Frozen-source, response-blind schema audit for original Haan et al. (2021) data.

This never tests host effects or derives butterfly survival or population fitness.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
import urllib.request
from pathlib import Path


def fetch(name: str, digest: str) -> tuple[bytes, str]:
    urls = [
        f"https://zenodo.org/records/4318182/files/{name}?download=1",
        f"https://zenodo.org/api/records/4318182/files/{name}/content",
    ]
    errors = []
    for url in urls:
        for attempt in range(2):
            try:
                req = urllib.request.Request(
                    url,
                    headers={
                        "User-Agent": "chocho-source-gate/0.1 (github.com/zuizui0223/chocho)",
                        "Accept": "text/csv,text/plain,application/octet-stream,*/*",
                    },
                )
                with urllib.request.urlopen(req, timeout=40) as response:
                    raw = response.read()
                actual = hashlib.md5(raw).hexdigest()
                if actual != digest:
                    raise RuntimeError("original input checksum mismatch: " + actual)
                return raw, url
            except Exception as exc:
                errors.append(f"{url}: {type(exc).__name__}: {str(exc)[:100]}")
                time.sleep(1 + attempt)
    raise RuntimeError(f"Could not authenticate original {name}: {errors}")


def audit_file(record: dict[str, str]) -> dict[str, object]:
    raw, url = fetch(record["name"], record["md5"])
    try:
        data = raw.decode("utf-8-sig")
    except UnicodeError as exc:
        raise RuntimeError(f"Undecodable original CSV {record['name']}") from exc
    reader = csv.reader(io.StringIO(data))
    try:
        header = next(reader)
    except StopIteration as exc:
        raise RuntimeError(f"Missing CSV header in {record['name']}") from exc
    if not header:
        raise RuntimeError(f"Missing CSV header in {record['name']}")
    repeated = sorted({field for field in header if header.count(field) > 1})
    count = 0
    for row in reader:
        if len(row) != len(header):
            raise RuntimeError(
                f"CSV source has non-rectangular records in {record['name']}"
            )
        count += 1
    if count == 0:
        raise RuntimeError(f"Empty original source {record['name']}")
    return {
        "file": record["name"],
        "file_role": record["role"],
        "md5": hashlib.md5(raw).hexdigest(),
        "bytes": len(raw),
        "rows": count,
        "columns": header,
        "duplicate_column_names": repeated,
        "fit_readiness": "HOLD_COLUMN_AMBIGUITY" if repeated else "SCHEMA_ONLY",
        "source_url": url,
        "response_values_examined": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if protocol["schema"] != "chocho_euphydryas_host_fallback_independent_source_gate_v0.1":
        raise RuntimeError("Incorrect source protocol")
    audited = [audit_file(f) for f in protocol["source"]["files"]]
    result = {
        "schema": "chocho_euphydryas_host_function_raw_source_schema_v0.1",
        "status": (
            "RAW_SOURCE_VERIFIED_BUT_DUPLICATE_COLUMNS_HOLD"
            if any(f["duplicate_column_names"] for f in audited)
            else "RAW_SOURCE_VERIFIED_SCHEMA_ONLY_NOT_A_NEW_MECHANISM_TEST"
        ),
        "sources": audited,
        "source_population": "Taylor's checkerspot in Haan et al. (2021), not Nevada population of Singer & Parmesan (2018)",
        "causal_boundary": [
            "A mapped species-level host is not proven locally usable for a particular insect population.",
            "The original experiments are already published, no biological response data are new.",
            "Adult recruitment is not identified by larval mass or short-term feeding preference.",
            "Do not claim a comparative host-fallback effect across dissimilar populations.",
        ],
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

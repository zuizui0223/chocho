#!/usr/bin/env python3
"""Verify a publicly shared annotated monarch larva-host database source.

Only file names, schemas, row counts and cryptographic hashes leave the run.
Individual iNaturalist records, observer identifiers, host labels, coordinates,
photographs and raw public source rows are not republished.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from zipfile import BadZipFile, ZipFile

SOURCE_PAGE = "https://www.freedman-lab.ca/resources"
SOURCE_ARTICLE = "https://doi.org/10.64898/2026.04.06.716508"
SHARE_URL = (
    "https://www.dropbox.com/scl/fo/"
    "bzr7nma5eetx34clz6ze4/AGNT2DKDUz_yk3Edg9o3gPo"
    "?rlkey=rgj83d80vh7pvhqm9dbn6012m&dl=1"
)
MAX_ARCHIVE_BYTES = 150_000_000
MAX_UNCOMPRESSED_BYTES = 400_000_000
SOURCE_SCHEMA = "chocho_monarch_annotated_host_public_metadata_20261028_v01"


def safe_csv_metadata(zf, member):
    with zf.open(member) as stream:
        reader = csv.reader(io.TextIOWrapper(stream, encoding="utf-8-sig", errors="replace", newline=""))
        try:
            header = next(reader)
        except StopIteration:
            return {"row_count": 0, "columns": []}
        rows = sum(1 for _ in reader)
    return {"row_count": rows, "columns": [x.strip()[:100] for x in header][:100]}


def inspect_zip(path):
    with ZipFile(path) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise ValueError("Zip source CRC failure")
        members = [x for x in archive.infolist() if not x.is_dir()]
        if not members:
            raise ValueError("Empty public source archive")
        if sum(x.file_size for x in members) > MAX_UNCOMPRESSED_BYTES:
            raise ValueError("Public source ZIP decompression size limit exceeded")
        out = []
        for entry in members:
            name = entry.filename
            p = Path(name)
            if p.is_absolute() or ".." in p.parts:
                raise ValueError("Unsafe ZIP member")
            entry_meta = {
                "source_filename": name[-220:],
                "bytes": entry.file_size,
                "suffix": p.suffix.lower(),
            }
            if p.suffix.lower() in {".csv", ".tsv"} and entry.file_size < 130_000_000:
                entry_meta.update(safe_csv_metadata(archive, entry))
            out.append(entry_meta)
    has_tabular = any(x.get("row_count", 0) > 0 for x in out)
    return {
        "status": "PUBLIC_SHARED_ARCHIVE_TABULAR_SOURCE_IDENTIFIED"
        if has_tabular else "PUBLIC_SHARED_ARCHIVE_NO_CSV_ROWS",
        "source_members": out,
        "total_tabular_rows_all_files": sum(x.get("row_count", 0) for x in out),
    }


def source_download(url, path, max_bytes=MAX_ARCHIVE_BYTES):
    request = Request(url, headers={"User-Agent": "chocho-monarch-public-source/1.0",
                                    "Accept": "application/zip,application/octet-stream,*/*"})
    actual = 0
    with urlopen(request, timeout=90) as response, path.open("wb") as output:
        host = urlparse(response.geturl()).hostname or ""
        if not (host == "dropbox.com" or host.endswith(".dropbox.com")
                or host == "dropboxusercontent.com" or host.endswith(".dropboxusercontent.com")):
            raise ValueError("Public shared source redirected outside Dropbox")
        status = response.status
        if status != 200:
            raise ValueError(f"Public shared source HTTP {status}")
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            actual += len(chunk)
            if actual > max_bytes:
                raise ValueError("Public shared source exceeds bounded file size")
            output.write(chunk)
        return {
            "final_host": host, "HTTP_status": status,
            "content_type": response.headers.get("Content-Type", ""),
            "response_bytes": actual,
        }


def run(output):
    from tempfile import TemporaryDirectory
    result = {
        "schema": SOURCE_SCHEMA, "source_page": SOURCE_PAGE,
        "source_article": SOURCE_ARTICLE,
        "source_description": "Micah Freedman public monarch caterpillar host annotations",
        "prior_published_results": {
            "annotated_record_count": 38856,
            "high_medium_certain_host_records": 30320,
            "non_native_only_winter_grid_cells": 209,
            "any_winter_grid_cells": 348,
            "reported_novelty_warning":
                "Freedman 2026 already reports nonnative host winter grid expansion; "
                "this cannot be claimed as chocho novelty",
        },
        "raw_observation_records_exported": False,
    }
    try:
        with TemporaryDirectory(prefix="chocho-plex-host-source-") as temp:
            local = Path(temp) / "shared_monarch.zip"
            result["transport"] = source_download(SHARE_URL, local)
            result["source_archive_sha256"] = hashlib.sha256(local.read_bytes()).hexdigest()
            result.update(inspect_zip(local))
    except (OSError, BadZipFile, ValueError, TimeoutError) as error:
        result["status"] = "SOURCE_BLOCKED_OR_NOT_A_VALID_TABULAR_ZIP"
        result["blocker_type"] = type(error).__name__
        result["blocker_detail"] = str(error)[:350]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, sort_keys=True, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({
        k: result.get(k) for k in
        ("status", "source_archive_sha256", "total_tabular_rows_all_files",
         "blocker_type", "blocker_detail")
    }, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    outcome = run(args.output)
    if outcome["status"] != "PUBLIC_SHARED_ARCHIVE_TABULAR_SOURCE_IDENTIFIED":
        raise SystemExit(2)

#!/usr/bin/env python3
"""Fetch the published EuDiS v5 CSV with source identity and checksum gates.

This is a trait-source acquisition step, not a species range or ocean-current test.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

ARTICLE = 22251028
VERSION = 5
API = f"https://api.figshare.com/v2/articles/{ARTICLE}/versions/{VERSION}"


def get(url: str) -> bytes:
    error = None
    for attempt in range(3):
        try:
            req = Request(url, headers={"User-Agent": "chocho-eudis-source-audit/2026.11", "Accept": "*/*"})
            with urlopen(req, timeout=90) as response:
                return response.read()
        except Exception as exc:
            error = exc
            if attempt < 2:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"External source unavailable at {url}: {error}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", required=True, type=Path)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    metadata = json.loads(get(API))
    if int(metadata.get("version", -1)) != VERSION or int(metadata.get("id", -1)) != ARTICLE:
        raise RuntimeError("Figshare returned a non-pinned EuDiS version or article")
    files = [f for f in metadata.get("files", []) if f["name"].lower().endswith(".csv")
             and "eudis" in f["name"].lower()]
    if len(files) != 1:
        raise RuntimeError(f"Expected exactly one EuDiS v5 CSV, found {[f['name'] for f in files]}")
    f = files[0]
    payload = get(f["download_url"])
    md5 = hashlib.md5(payload).hexdigest()
    if md5.lower() != str(f.get("computed_md5") or f.get("supplied_md5") or "").lower():
        raise RuntimeError("EuDiS CSV does not match Figshare's file MD5")
    dest = a.out_dir / "EuDiS.csv"
    dest.write_bytes(payload)
    with dest.open("r", encoding="utf-8-sig", newline="") as inp:
        cols = next(csv.reader(inp))
    if not {"Species_POWO", "Thalassochorous", "Further_clarification_needed"}.issubset(set(cols)):
        raise RuntimeError(f"Published EuDiS columns differ: {cols}")
    receipt = {
        "article_id": ARTICLE, "version": VERSION,
        "metadata_url": API, "doi": "10.6084/m9.figshare.22251028.v5",
        "file_id": f["id"], "original_filename": f["name"],
        "source_download_url": f["download_url"],
        "bytes": len(payload), "source_md5": md5,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "columns": cols, "status": "PINNED_SOURCE_TRAITS_NOT_OBSERVED_DISPERSAL"
    }
    (a.out_dir / "eudis_v5_source_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("EUDIS_PINNED_SOURCE", json.dumps({
        "version": VERSION, "file_id": f["id"], "bytes": len(payload), "sha256": receipt["sha256"],
        "columns": len(cols)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

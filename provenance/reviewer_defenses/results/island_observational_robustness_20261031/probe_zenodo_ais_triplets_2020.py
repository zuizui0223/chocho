#!/usr/bin/env python3
"""Validate provenance/schema of an open 2020 AIS ship-port triplet archive.

DOI 10.5281/zenodo.10783698, CC BY 4.0. Data from 2020 DO NOT establish
2001-2010 host introduction history. Inspect ZIP metadata and header only;
never publish ship identifiers, raw trajectories or individual voyages.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

ARCHIVE = Path("_private_source_stage/triplet_2020.zip")
EXPECTED_MD5 = "af27e1bef1f706e2e29ed53aee47ff10"
OUT = Path("provenance/reviewer_defenses/results/island_observational_robustness_20261031/"
           "zenodo_ais_2020_archive_schema.json")


def inspect(archive=ARCHIVE):
    if not archive.exists():
        raise ValueError("Expected public AIS ZIP not downloaded")
    h = hashlib.md5()
    with archive.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    if h.hexdigest() != EXPECTED_MD5:
        raise ValueError("Zenodo archive MD5 not matching source-record metadata")
    manifest = []
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ValueError("Invalid ZIP member checksum")
        for file in z.infolist():
            item = {
                "filename": file.filename, "size_uncompressed": file.file_size,
                "compressed_size": file.compress_size, "is_dir": file.is_dir(),
            }
            if not file.is_dir() and file.file_size:
                with z.open(file) as handle:
                    sample = handle.read(35000).decode("utf-8-sig", errors="replace")
                nonblank = [line for line in sample.splitlines() if line.strip()]
                first = nonblank[0] if nonblank else ""
                item["first_line_fields"] = len(first.split(",")) if first else 0
                item["first_line_separator_counts"] = {
                    "comma": first.count(","), "semicolon": first.count(";"),
                    "tab": first.count("\t"), "pipe": first.count("|"),
                }
                item["first_line_lengths"] = [len(x) for x in nonblank[:3]]
                # Export ONLY field names, NEVER the data rows. The source
                # proves semicolon-delimited CSV, not a comma-separated file.
                dialect = ";" if first.count(";") > first.count(",") else ","
                header = [x.strip().strip('"').strip() for x in first.split(dialect)]
                if header and not header[-1]:
                    header.pop()  # Source ends each record with a trailing semicolon.
                # The verified CSV has 12+ recognizable schema words and ZERO
                # numeric characters in the initial row, which is unambiguously
                # a header despite diacritics and unusual separators.
                keywords = ("date", "time", "ship", "port", "mmsi", "imo",
                            "latitude", "longitude", "lon", "lat", "country", "flag",
                            "vessel", "type", "arrival", "departure")
                hits = sum(any(k in t.lower() for k in keywords) for t in header)
                is_header = (len(header) >= 6 and len(set(header)) == len(header)
                             and hits >= 3 and not any(ch.isdigit() for t in header for ch in t)
                             and max(map(len, header)) < 90)
                item["has_verified_header"] = is_header
                item["column_names"] = header if is_header else []
                item["delimiter"] = dialect
                # The rare header with spaces/diacritics may fail strict regex.
                # Classify columns without ever exporting an individual row.
                item["first_line_token_profile"] = [
                    {"index": j + 1, "length": len(t),
                     "letters": sum(ch.isalpha() for ch in t),
                     "digits": sum(ch.isdigit() for ch in t),
                     "punctuation": sum(not ch.isalnum() and not ch.isspace() for ch in t),
                     "header_keyword": next(
                         (k for k in ("date","time","ship","port","mmsi","imo",
                                      "latitude","longitude","lon","lat","country","flag",
                                      "vessel","type","arrival","departure")
                          if k in t.lower()), None)}
                    for j, t in enumerate(header)
                ]
                item["header_keyword_hits"] = sum(
                    bool(x["header_keyword"]) for x in item["first_line_token_profile"])
                item["first_line_is_probably_observations"] = (
                    not is_header and item["header_keyword_hits"] < 2)

                # No sample observations, individual vessel IDs or MMSI.
                item["sample_raw_values_exported"] = False
            manifest.append(item)
    return {
        "source": "Duvallet et al. 2024",
        "doi": "10.5281/zenodo.10783698",
        "source_url": "https://zenodo.org/records/10783698/files/triplet_2020.zip?download=1",
        "archive_md5": h.hexdigest(),
        "license": "CC BY 4.0",
        "archive_size_bytes": archive.stat().st_size,
        "members": manifest,
        "status": "SOURCE_SCHEMA_ONLY_NO_VOYAGE_OR_BUTTERFLY_TEST",
        "limits": [
            "2020 AIS indices are retrospective relative to some 2011-2022 butterflies.",
            "A date-ship-port triplet is not itself a validated consecutive-port voyage.",
            "Port IDs must be mapped to WPI 2019 coordinates/landmasses independently.",
            "The global completeness and ship types of AIS sources must be verified.",
            "No valid port-to-port island similarity effect can be concluded from this ZIP audit.",
        ],
    }


def main():
    result = inspect()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf8")
    print(json.dumps({k: v for k, v in result.items() if k != "members"},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()

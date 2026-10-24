#!/usr/bin/env python3
"""Verify recovery keys against the archived original GBIF source before fetching.

No live taxonomy lookup, occurrence acquisition, or source-file modification is
performed. The source manifest stays outside the repository; the CLI prints a
small validation receipt with input hashes only after every pair matches.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import sys


SCHEMA = "chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01"
EXPECTED_STATUSES = {
    "COMPLETE": 135,
    "PARTIAL_FETCH_NOT_ELIGIBLE": 96,
    "TAXON_UNRESOLVED": 8,
}
FIELDS = ("species", "accepted_gbif_species_key")


def positive_taxon_key(value, *, source, species):
    text = str(value)
    if not text.isascii() or not text.isdecimal() or int(text) <= 0:
        raise ValueError(f"Invalid {source} taxon key for {species}: {value!r}")
    return int(text)


def validate(original_manifest, keys):
    """Return a receipt if the CSV is exactly the original 96 partial taxa."""
    manifest_bytes = Path(original_manifest).read_bytes()
    key_bytes = Path(keys).read_bytes()
    manifest = json.loads(manifest_bytes)
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ValueError("Original source schema is not the archived 239-taxon schema")
    if (manifest.get("species_requested") != 239 or
            manifest.get("species_statuses") != EXPECTED_STATUSES):
        raise ValueError("Original source must retain the recorded 239 taxa and 135/96/8 status counts")
    receipts = manifest.get("species_receipts")
    if not isinstance(receipts, list) or len(receipts) != 239:
        raise ValueError("Original source receipts must contain exactly 239 taxa")

    names = set()
    status_counts = Counter()
    expected = {}
    source_key_names = {}
    for row in receipts:
        if not isinstance(row, dict):
            raise ValueError("Invalid original species receipt")
        species = row.get("species")
        if not isinstance(species, str) or not species or species != species.strip():
            raise ValueError(f"Invalid original species name: {species!r}")
        if species in names:
            raise ValueError(f"Duplicate original species: {species}")
        names.add(species)
        status = row.get("status")
        if status not in EXPECTED_STATUSES:
            raise ValueError(f"Invalid original species status for {species}: {status!r}")
        status_counts[status] += 1
        if status == "PARTIAL_FETCH_NOT_ELIGIBLE":
            taxon = row.get("taxon")
            key = positive_taxon_key(
                taxon.get("key") if isinstance(taxon, dict) else None,
                source="original", species=species,
            )
            if key in source_key_names:
                raise ValueError(
                    f"Duplicate original partial taxon key {key}: "
                    f"{source_key_names[key]}, {species}"
                )
            source_key_names[key] = species
            expected[species] = key
    if dict(status_counts) != EXPECTED_STATUSES:
        raise ValueError(f"Original receipt status counts disagree with 135/96/8: {dict(status_counts)}")

    reader = csv.DictReader(io.StringIO(key_bytes.decode("utf-8"), newline=""))
    if tuple(reader.fieldnames or ()) != FIELDS:
        raise ValueError(f"Recovery CSV schema must be exactly {','.join(FIELDS)}")
    recovered = {}
    recovery_key_names = {}
    for row in reader:
        if None in row:
            raise ValueError(f"Recovery CSV row {reader.line_num} has extra fields")
        species = row["species"]
        if not isinstance(species, str) or not species or species != species.strip():
            raise ValueError(f"Invalid recovery species name at row {reader.line_num}: {species!r}")
        if species in recovered:
            raise ValueError(f"Duplicate recovery species: {species}")
        key = positive_taxon_key(
            row["accepted_gbif_species_key"], source="recovery", species=species,
        )
        if key in recovery_key_names:
            raise ValueError(
                f"Duplicate recovery taxon key {key}: {recovery_key_names[key]}, {species}"
            )
        recovery_key_names[key] = species
        recovered[species] = key

    missing = sorted(set(expected) - set(recovered))
    unexpected = sorted(set(recovered) - set(expected))
    problems = []
    if missing:
        problems.append(f"Missing recovery species: {', '.join(missing)}")
    if unexpected:
        problems.append(f"Unexpected recovery species: {', '.join(unexpected)}")
    if problems:
        raise ValueError("; ".join(problems))
    mismatches = [
        f"{species}: original {expected[species]}, recovery {recovered[species]}"
        for species in sorted(expected) if expected[species] != recovered[species]
    ]
    if mismatches:
        raise ValueError("Recovery taxon key mismatch: " + "; ".join(mismatches))
    return {
        "status": "PASS_ORIGINAL_96_TAXON_KEYS",
        "original_species_verified": len(names),
        "recovery_species_verified": len(recovered),
        "original_status_counts": dict(status_counts),
        "original_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "recovery_keys_sha256": hashlib.sha256(key_bytes).hexdigest(),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--original-manifest", required=True, type=Path,
        help="Archived 239-species receipt from original GBIF run 39233336804",
    )
    parser.add_argument(
        "--keys", required=True, type=Path,
        help="CSV of the 96 incomplete species and their original accepted GBIF keys",
    )
    args = parser.parse_args(argv)
    try:
        receipt = validate(args.original_manifest, args.keys)
    except (OSError, ValueError, csv.Error) as exc:
        print(f"FAIL_ORIGINAL_96_TAXON_KEYS: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import zipfile
from pathlib import Path


EXPECTED_OCCURRENCES_SHA256 = (
    "1dc71938fd5fbed48756cc8bfd0e4bb9580f7d20eebfb652547febcc2893dcad"
)
EXPECTED_RECORDS = 53_434


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def occurrence_bytes(path: Path) -> bytes:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            return zf.read("occurrences.csv")
    return path.read_bytes()


def build_request(data: bytes, email: str) -> dict:
    digest = sha256_bytes(data)
    if digest != EXPECTED_OCCURRENCES_SHA256:
        raise RuntimeError(
            f"unexpected occurrences.csv SHA256 {digest}; "
            f"expected {EXPECTED_OCCURRENCES_SHA256}"
        )

    rows = list(csv.DictReader(data.decode("utf-8").splitlines()))
    keys = [str(row["gbif_key"]).strip() for row in rows]
    if len(keys) != EXPECTED_RECORDS:
        raise RuntimeError(f"expected {EXPECTED_RECORDS} rows, found {len(keys)}")
    if len(set(keys)) != len(keys):
        raise RuntimeError("GBIF keys are not unique")

    return {
        "notificationAddresses": [email],
        "sendNotification": True,
        "format": "SIMPLE_CSV",
        "predicate": {
            "type": "in",
            "key": "GBIF_ID",
            "values": keys,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--occurrences",
        type=Path,
        required=True,
        help="occurrences.csv or the historical summary ZIP containing occurrences.csv",
    )
    ap.add_argument("--email", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    payload = build_request(occurrence_bytes(args.occurrences), args.email)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "gbif_ids": len(payload["predicate"]["values"]),
                "predicate": "GBIF_ID",
                "format": payload["format"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

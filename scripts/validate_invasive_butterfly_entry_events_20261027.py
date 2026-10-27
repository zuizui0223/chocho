#!/usr/bin/env python3
"""Validate literature-anchored butterfly recipient events against original 239.

First detection and population establishment are separate observables.
The registry itself is a historical-source curation, not proof of source-host use.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

COLS = (
    "species", "in_original_239", "native_region_statement",
    "recipient_region", "recipient_code", "first_year_lower",
    "first_year_upper", "establishment_by_year", "date_interpretation",
    "arrival_mechanism_status", "source_confidence",
    "candidate_host_context", "primary_evidence_url", "interpretive_limits",
)
SCHEMA = "chocho_invasive_butterfly_evidence_register_20261027_v01"

def read_register(path):
    with path.open(encoding="utf-8", newline="") as f:
        rows = csv.DictReader(f)
        if tuple(rows.fieldnames or []) != COLS:
            raise ValueError("Unexpected frozen invasive butterfly entry event columns")
        return list(rows)


def evaluate(register_path, source_species):
    raw = register_path.read_bytes()
    events = read_register(register_path)
    if len(source_species) != 239 or len(set(source_species)) != 239:
        raise ValueError("Not the original 239-butterfly source taxon panel")
    ids, seen, by_species = [], set(), defaultdict(list)
    for i, row in enumerate(events):
        sp, code = row["species"], row["recipient_code"]
        if not sp or not code or (sp, code) in seen:
            raise ValueError("Missing or repeated recipient event identity")
        seen.add((sp, code))
        truth = row["in_original_239"]
        if truth not in ("true", "false"):
            raise ValueError("Original 239 boolean must be explicit")
        if (sp in source_species) != (truth == "true"):
            raise ValueError("Original source 239 species eligibility drift " + sp)
        try:
            first = int(row["first_year_lower"])
            end = int(row["first_year_upper"])
        except ValueError as e:
            raise ValueError("Missing or invalid initial documented year") from e
        if not 1800 <= first <= end <= 2026:
            raise ValueError("Impossible first reporting window " + sp)
        est = row["establishment_by_year"]
        if est:
            est_year = int(est)
            if est_year < first or est_year > 2026:
                raise ValueError("Established year predates first known record")
        u = urlparse(row["primary_evidence_url"])
        if not (u.scheme == "https" and u.netloc and u.path):
            raise ValueError("Missing external historical source")
        if row["source_confidence"] not in {"high", "medium", "low"}:
            raise ValueError("Unknown source confidence")
        if not row["native_region_statement"] or not row["date_interpretation"]:
            raise ValueError("Missing invasion/native geography qualification")
        by_species[sp].append(row)
        ids.append({
            "species": sp, "recipient_code": code,
            "year": [first, end],
            "source_239": truth == "true",
            "year_type": row["date_interpretation"],
            "route_status": row["arrival_mechanism_status"],
        })
    if len(events) != 13 or len(by_species) != 8 or \
       len({s for s in by_species if s in source_species}) != 6:
        raise ValueError("The preselected 8-species / 13-event roster changed")
    if set(by_species) - source_species != {"Cacyreus marshalli", "Chilades pandava"}:
        raise ValueError("Missing named source-offpanel test species")
    if any(e["species"] == "Danaus plexippus" and
           e["route_status"] == "confirmed_human_transport" for e in ids):
        raise ValueError("Uncertain monarch dispersal must not be classified human-caused")
    flags = Counter("core_239" if s in source_species else "offpanel"
                    for s in by_species)
    outcome = {
        "schema": SCHEMA,
        "status": "HISTORICAL_ENTRY_EVENTS_SOURCE_CHECKED_NOT_ECOLOGICAL_CAUSAL_PROOF",
        "input_csv_sha256": hashlib.sha256(raw).hexdigest(),
        "source_original_239_species_verified": len(source_species),
        "distinct_candidate_butterfly_species": len(by_species),
        "candidate_in_original_239": flags["core_239"],
        "candidate_outside_original_239": flags["offpanel"],
        "distinct_species_recipient_event_count": len(events),
        "events_by_species": dict(sorted(
            Counter(x["species"] for x in ids).items())),
        "event_identities_no_raw_coordinates": ids,
        "key_cautions": [
            "First report year is not actual introduction date or necessarily first successful establishment.",
            "Danaus plexippus remote island origins may be naturally dispersed: not assumed anthropogenic introduction.",
            "Host plants listed in the candidate registry are plausible context, not source-confirmed invader larval adoption.",
            "Original 239 subset excludes Cacyreus marshalli and Chilades pandava: add them as separate verified extensions.",
            "Published Pieris rapae Australia dates 1929 first specimen versus c.1937 established are different outcomes.",
            "No evidence that original plant native range, consumer native range and new observed host pairs are the same mechanism.",
        ],
    }
    return outcome


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--source-original-species", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    with args.source_original_species.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    source_species = [r["species"] for r in rows]
    result = evaluate(args.csv, source_species)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in (
        "status", "candidate_in_original_239",
        "candidate_outside_original_239", "distinct_candidate_butterfly_species",
        "distinct_species_recipient_event_count")}, sort_keys=True))


if __name__ == "__main__":
    main()

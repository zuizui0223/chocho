#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from acquire_butterfly_resource_envelope_occurrences import (
    atomic_json,
    combine_records,
    fetch_species_pages,
)
from butterfly_specialization_ecology.checkpointed_gbif_occurrence import (
    occurrence_window_chunks,
    species_state_key,
)

EXPECTED_SOURCE_SHA256 = (
    "b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa"
)
EXPECTED_SPECIES_COUNT = 239
EXPECTED_SPECIES_SHA256 = (
    "c972e145e160f0c34a5fcee9f423e11c8e38072038e9765a40e79c0c5cd06157"
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_panel(path: Path) -> tuple[str, ...]:
    if sha256_path(path) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("host-contribution source SHA drift")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    names = tuple(sorted(str(row["species"]).strip() for row in rows))
    if len(names) != EXPECTED_SPECIES_COUNT or len(set(names)) != EXPECTED_SPECIES_COUNT:
        raise RuntimeError("expected exact 239-species resource panel")
    digest = hashlib.sha256(("\n".join(names) + "\n").encode("utf-8")).hexdigest()
    if digest != EXPECTED_SPECIES_SHA256:
        raise RuntimeError("239-species panel digest drift")
    return names


def write_occurrences(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "species",
        "gbif_key",
        "latitude",
        "longitude",
        "year",
        "basis_of_record",
        "dataset_key",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def repair_completed_chunk_windows(
    species: str,
    state_root: Path,
    *,
    chunk_size: int,
) -> dict[str, int]:
    """Rebuild a fixed ordinal window from completed transport chunks.

    GBIF occurrence-search paging is not guaranteed to be stable across
    separately transported chunks. If the same GBIF key appears in more than
    one completed chunk, retain its first fetched record and record the
    duplicate count. This changes only transport assembly; the six frozen
    ordinal parent windows and all scientific filters remain unchanged.
    """
    state_dir = state_root / species_state_key(species)
    meta_path = state_dir / "metadata.json"
    if not meta_path.exists():
        return {"pages_rebuilt": 0, "duplicate_keys_removed": 0}
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta.get("status") == "REJECTED_GBIF_TAXON_MATCH":
        return {"pages_rebuilt": 0, "duplicate_keys_removed": 0}

    total = int(meta.get("total_coordinate_records_2010_2026") or 0)
    page_size = int(meta.get("page_size") or 300)
    pages_dir = state_dir / "pages"
    pages_dir.mkdir(exist_ok=True)
    rebuilt = 0
    duplicate_total = 0

    for offset in map(int, meta.get("page_offsets") or []):
        page_path = pages_dir / f"offset_{offset:06d}.json"
        if page_path.exists():
            continue
        window_size = max(0, min(page_size, total - offset))
        chunks = occurrence_window_chunks(
            offset,
            window_size,
            chunk_size=int(chunk_size),
        )
        records = []
        complete = True
        for chunk_offset, chunk_limit in chunks:
            chunk_path = (
                pages_dir
                / f"offset_{offset:06d}_chunks"
                / f"offset_{chunk_offset:06d}_limit_{chunk_limit:03d}.json"
            )
            if not chunk_path.exists():
                complete = False
                break
            payload = json.loads(chunk_path.read_text(encoding="utf-8"))
            records.extend(payload.get("records", []))
        if not complete:
            continue

        unique = {}
        duplicates = 0
        for record in records:
            key = int(record["key"])
            if key in unique:
                duplicates += 1
                continue
            unique[key] = record
        atomic_json(
            page_path,
            {
                "species": species,
                "usage_key": int(meta["usage_key"]),
                "offset": int(offset),
                "window_size": int(window_size),
                "transport_chunk_size": int(chunk_size),
                "records": list(unique.values()),
                "realm_audit_transport_repair": {
                    "completed_chunks_reassembled": True,
                    "duplicate_gbif_keys_removed": int(duplicates),
                    "scientific_window_changed": False,
                },
            },
        )
        rebuilt += 1
        duplicate_total += duplicates

    return {
        "pages_rebuilt": int(rebuilt),
        "duplicate_keys_removed": int(duplicate_total),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-csv", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--shard-count", type=int, required=True)
    ap.add_argument("--state-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--species-seconds", type=float, default=900.0)
    ap.add_argument("--request-seconds", type=float, default=90.0)
    ap.add_argument("--maximum-pages", type=int, default=6)
    ap.add_argument("--transport-chunk-size", type=int, default=300)
    ap.add_argument("--passes", type=int, default=2)
    args = ap.parse_args()

    if args.shard_count < 1 or not (0 <= args.shard_index < args.shard_count):
        raise RuntimeError("invalid shard specification")
    if args.maximum_pages != 6:
        raise RuntimeError("GEB realm audit is frozen at six ordinal windows")

    panel = load_panel(args.source_csv)
    species = tuple(
        name for i, name in enumerate(panel)
        if i % args.shard_count == args.shard_index
    )
    if not species:
        raise RuntimeError("empty shard")

    args.state_dir.mkdir(parents=True, exist_ok=True)
    final = {}
    repair_totals = {
        name: {"pages_rebuilt": 0, "duplicate_keys_removed": 0}
        for name in species
    }
    for pass_index in range(args.passes):
        for name in species:
            prior = final.get(name)
            if prior is not None and prior.get("status") in {
                "COMPLETE",
                "REJECTED_GBIF_TAXON_MATCH",
            }:
                continue
            result = fetch_species_pages(
                name,
                args.state_dir,
                species_seconds=args.species_seconds,
                maximum_pages=args.maximum_pages,
                request_seconds=args.request_seconds,
                transport_chunk_size=args.transport_chunk_size,
            )
            # A completed set of smaller transport chunks can contain repeated
            # GBIF keys because occurrence-search order is not a formal stable
            # paging contract. Reassemble that frozen parent window by key and
            # continue; do not change offsets, filters, or species membership.
            for _ in range(args.maximum_pages + 1):
                if result.get("status") != "PARTIAL_TRANSPORT":
                    break
                repaired = repair_completed_chunk_windows(
                    name,
                    args.state_dir,
                    chunk_size=args.transport_chunk_size,
                )
                repair_totals[name]["pages_rebuilt"] += repaired["pages_rebuilt"]
                repair_totals[name]["duplicate_keys_removed"] += repaired[
                    "duplicate_keys_removed"
                ]
                if repaired["pages_rebuilt"] == 0:
                    break
                result = fetch_species_pages(
                    name,
                    args.state_dir,
                    species_seconds=args.species_seconds,
                    maximum_pages=args.maximum_pages,
                    request_seconds=args.request_seconds,
                    transport_chunk_size=args.transport_chunk_size,
                )
            final[name] = result

    rows = combine_records(args.state_dir, species)
    write_occurrences(args.output_csv, rows)

    counts = {}
    for row in rows:
        name = str(row["species"])
        counts[name] = counts.get(name, 0) + 1

    species_ledger = []
    for name in species:
        row = dict(final[name])
        row["realm_audit_transport_repair"] = dict(repair_totals[name])
        species_ledger.append(row)

    ledger = {
        "schema": "chocho_geb_occurrence_realm_shard_v0.1",
        "status": (
            "COMPLETE"
            if all(
                row["status"] in {"COMPLETE", "REJECTED_GBIF_TAXON_MATCH"}
                for row in species_ledger
            )
            else "RESUMABLE_PARTIAL_TRANSPORT"
        ),
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "panel_species_count": len(panel),
        "panel_species_sha256": EXPECTED_SPECIES_SHA256,
        "species": species_ledger,
        "occurrence_rows": len(rows),
        "occurrence_rows_by_species": dict(sorted(counts.items())),
        "transport": {
            "period": "2010-2026",
            "maximum_ordinal_windows_per_species": 6,
            "page_size": 300,
            "transport_chunk_size": args.transport_chunk_size,
            "passes": args.passes,
            "scientific_filters_identical_to_existing_chocho_climate_route": True,
        },
        "claim_boundary": (
            "Post-hoc realized-geography generality audit only; no frozen ecological "
            "result or causal claim is altered."
        ),
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(
        json.dumps(ledger, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": ledger["status"],
        "shard": args.shard_index,
        "species": len(species),
        "rows": len(rows),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

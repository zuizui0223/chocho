#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

from acquire_butterfly_resource_envelope_occurrences import get_json
from butterfly_specialization_ecology.checkpointed_gbif_occurrence import (
    occurrence_window_chunks,
)


FIELDS = [
    "species",
    "gbif_key",
    "latitude",
    "longitude",
    "year",
    "basis_of_record",
    "dataset_key",
]


def load_rows(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        out = []
        for row in csv.DictReader(handle):
            out.append(
                {
                    "species": str(row["species"]),
                    "gbif_key": int(row["gbif_key"]),
                    "latitude": float(row["latitude"]),
                    "longitude": float(row["longitude"]),
                    "year": row.get("year"),
                    "basis_of_record": row.get("basis_of_record"),
                    "dataset_key": row.get("dataset_key"),
                }
            )
    return out


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    unique = {
        (str(row["species"]), int(row["gbif_key"])): row
        for row in rows
    }
    ordered = [unique[key] for key in sorted(unique)]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(ordered)


def recover_window(
    species: str,
    usage_key: int,
    offset: int,
    *,
    deadline: float,
    request_seconds: float,
    chunk_size: int,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for chunk_offset, chunk_limit in occurrence_window_chunks(
        offset, 300, chunk_size=chunk_size
    ):
        payload = get_json(
            "occurrence/search",
            {
                "taxonKey": int(usage_key),
                "hasCoordinate": "true",
                "hasGeospatialIssue": "false",
                "occurrenceStatus": "PRESENT",
                "year": "2010,2026",
                "limit": int(chunk_limit),
                "offset": int(chunk_offset),
            },
            deadline=deadline,
            request_seconds=request_seconds,
        )
        for item in payload.get("results") or []:
            key = int(item.get("key") or 0)
            lat = item.get("decimalLatitude")
            lon = item.get("decimalLongitude")
            if key <= 0 or lat is None or lon is None:
                continue
            records.append(
                {
                    "species": species,
                    "gbif_key": key,
                    "latitude": float(lat),
                    "longitude": float(lon),
                    "year": item.get("year"),
                    "basis_of_record": item.get("basisOfRecord"),
                    "dataset_key": item.get("datasetKey"),
                }
            )
    return records


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", type=Path, required=True)
    ap.add_argument("--input-ledger", type=Path, required=True)
    ap.add_argument("--recovery-contract", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    contract = json.loads(args.recovery_contract.read_text(encoding="utf-8"))
    if contract.get("schema") != "chocho_geb_realm_transport_recovery_v0.1":
        raise RuntimeError("unexpected recovery contract")
    if contract.get("status") != (
        "FROZEN_TECHNICAL_RECOVERY_AFTER_PARTIAL_TRANSPORT_BEFORE_FULL_REALM_RESULT"
    ):
        raise RuntimeError("recovery contract is not frozen")

    ledger = json.loads(args.input_ledger.read_text(encoding="utf-8"))
    if ledger.get("schema") != "chocho_geb_occurrence_realm_shard_v0.1":
        raise RuntimeError("unexpected source shard ledger")

    rows = load_rows(args.input_csv)
    cfg = contract["transport_only_changes"]
    passes = int(cfg["passes"])
    request_seconds = float(cfg["request_timeout_cap_seconds"])
    species_seconds = float(cfg["species_total_deadline_seconds"])
    chunk_size = int(cfg["chunk_size"])

    out_species = []
    recovery_events = []
    for original in ledger["species"]:
        state = dict(original)
        name = str(state["species"])
        original_status = str(state["status"])

        if original_status != "PARTIAL_TRANSPORT":
            out_species.append(state)
            continue

        usage_key = int(state.get("usage_key") or 0)
        if usage_key <= 0:
            raise RuntimeError(f"{name}: partial transport without usage_key")
        planned = tuple(map(int, state.get("planned_offsets") or ()))
        missing = list(map(int, state.get("missing_offsets") or ()))
        completed = set(map(int, state.get("completed_offsets") or ()))
        if not missing:
            raise RuntimeError(f"{name}: PARTIAL_TRANSPORT but no missing offsets")
        if any(offset not in planned for offset in missing):
            raise RuntimeError(f"{name}: missing offset outside fixed plan")

        deadline = time.monotonic() + species_seconds
        errors: dict[int, str] = {}
        recovered: set[int] = set()
        for pass_index in range(1, passes + 1):
            for offset in list(missing):
                try:
                    page_rows = recover_window(
                        name,
                        usage_key,
                        offset,
                        deadline=deadline,
                        request_seconds=request_seconds,
                        chunk_size=chunk_size,
                    )
                except Exception as exc:
                    errors[offset] = f"{type(exc).__name__}: {exc}"
                    continue
                rows.extend(page_rows)
                completed.add(offset)
                recovered.add(offset)
                missing.remove(offset)
                errors.pop(offset, None)
            if not missing:
                break

        state["completed_offsets"] = sorted(completed)
        state["missing_offsets"] = sorted(missing)
        state["status"] = "COMPLETE" if not missing else "PARTIAL_TRANSPORT"
        state["error"] = None if not missing else "; ".join(
            f"{offset}:{errors.get(offset, 'unrecovered')}"
            for offset in sorted(missing)
        )
        state["technical_recovery"] = {
            "contract": str(args.recovery_contract),
            "original_status": original_status,
            "recovered_offsets": sorted(recovered),
            "remaining_missing_offsets": sorted(missing),
            "chunk_size": chunk_size,
            "request_timeout_cap_seconds": request_seconds,
            "species_total_deadline_seconds": species_seconds,
            "passes": passes,
        }
        recovery_events.append(
            {
                "species": name,
                "original_missing_offsets": list(map(int, original["missing_offsets"])),
                "recovered_offsets": sorted(recovered),
                "remaining_missing_offsets": sorted(missing),
                "final_status": state["status"],
            }
        )
        out_species.append(state)

    write_rows(args.output_csv, rows)
    status = (
        "COMPLETE"
        if all(
            row["status"] in {"COMPLETE", "REJECTED_GBIF_TAXON_MATCH"}
            for row in out_species
        )
        else "RESUMABLE_PARTIAL_TRANSPORT"
    )
    payload = dict(ledger)
    payload["status"] = status
    payload["species"] = out_species
    payload["occurrence_rows"] = sum(1 for _ in csv.DictReader(
        args.output_csv.open(newline="", encoding="utf-8")
    ))
    payload["technical_recovery"] = {
        "schema": contract["schema"],
        "source_workflow_run_id": contract["source"]["workflow_run_id"],
        "events": recovery_events,
        "scientific_query_changed": False,
        "species_replaced": False,
        "realm_result_used_for_recovery": False,
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status": status,
                "shard": ledger["shard_index"],
                "partial_species_before": sum(
                    row["status"] == "PARTIAL_TRANSPORT" for row in ledger["species"]
                ),
                "partial_species_after": sum(
                    row["status"] == "PARTIAL_TRANSPORT" for row in out_species
                ),
                "recovery_events": recovery_events,
            },
            sort_keys=True,
        )
    )
    return 0 if status == "COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())

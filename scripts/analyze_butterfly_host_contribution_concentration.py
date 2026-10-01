#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


EXPECTED_DESCRIPTOR_SHA256 = (
    "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_descriptors(path: Path) -> dict[str, dict[str, object]]:
    if sha256_path(path) != EXPECTED_DESCRIPTOR_SHA256:
        raise RuntimeError("descriptor SHA-256 drift")
    out: dict[str, dict[str, object]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"species", "host_family_count", "host_wgsrpd3_unit_count"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in reader:
            name = str(row["species"]).strip()
            out[name] = {
                "host_family_count": float(row["host_family_count"]),
                "native_resource_units": int(row["host_wgsrpd3_unit_count"]),
            }
    if len(out) != 339:
        raise RuntimeError(f"expected 339 descriptor species; got {len(out)}")
    return out


def load_pairs(path: Path):
    by_butterfly: dict[str, set[str]] = defaultdict(set)
    plant_meta: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "insect_species",
            "accepted_plant_name_id",
            "accepted_name",
            "family",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for row in reader:
            butterfly = str(row["insect_species"]).strip()
            plant_id = str(row["accepted_plant_name_id"]).strip()
            accepted_name = str(row["accepted_name"]).strip()
            family = str(row["family"]).strip()
            if not butterfly or not plant_id:
                continue
            by_butterfly[butterfly].add(plant_id)
            genus = accepted_name.split()[0] if accepted_name else ""
            meta = {
                "accepted_name": accepted_name,
                "genus": genus,
                "family": family,
            }
            previous = plant_meta.setdefault(plant_id, meta)
            if previous != meta:
                raise RuntimeError(f"plant metadata conflict for {plant_id}")
    return by_butterfly, plant_meta


def load_units(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("distribution sidecar schema drift")
        for row in reader:
            plant_id = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if plant_id and unit:
                out[plant_id].add(unit)
    return out


def gini(values: list[float]) -> float:
    values = sorted(float(v) for v in values if v > 0)
    if not values:
        return 0.0
    n = len(values)
    total = sum(values)
    return (
        2.0 * sum((i + 1) * value for i, value in enumerate(values))
        / (n * total)
        - (n + 1) / n
    )


def summarize_credit(credit: dict[str, float]) -> dict[str, object]:
    ranked = sorted(credit.items(), key=lambda kv: (-kv[1], kv[0]))
    total = sum(value for _, value in ranked)
    shares = [value / total for _, value in ranked] if total else []
    top = {}
    for k in (1, 5, 10, 20, 50, 100):
        if ranked:
            top[str(k)] = sum(value for _, value in ranked[:k]) / total
    hhi = sum(share * share for share in shares)
    return {
        "contributors": len(ranked),
        "total_fractional_credit": total,
        "top_k_share": top,
        "gini": gini([value for _, value in ranked]),
        "hhi": hhi,
        "effective_contributor_number": None if hhi == 0 else 1.0 / hhi,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-plant-csv", type=Path, required=True)
    ap.add_argument("--output-genus-csv", type=Path, required=True)
    args = ap.parse_args()

    descriptors = load_descriptors(args.descriptors_csv)
    pairs, plant_meta = load_pairs(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)

    plant_credit: dict[str, float] = defaultdict(float)
    genus_credit: dict[str, float] = defaultdict(float)
    expanded_butterflies = 0
    total_added_units = 0

    for butterfly, descriptor in descriptors.items():
        if (
            float(descriptor["host_family_count"]) <= 0
            or int(descriptor["native_resource_units"]) <= 0
        ):
            continue
        hosts = sorted(pairs.get(butterfly, set()))
        if not hosts:
            raise RuntimeError(
                f"resource-eligible butterfly lacks resolved hosts: {butterfly}"
            )

        native_union: set[str] = set()
        contemporary_union: set[str] = set()
        added_by_host: dict[str, set[str]] = {}
        for host in hosts:
            n = set(native.get(host, set()))
            c = set(contemporary.get(host, set()))
            if not n <= c:
                raise RuntimeError(f"native not subset contemporary for plant {host}")
            native_union.update(n)
            contemporary_union.update(c)
            added_by_host[host] = c - n

        if len(native_union) != int(descriptor["native_resource_units"]):
            raise RuntimeError(
                f"native union drift for {butterfly}: "
                f"{len(native_union)} vs {descriptor['native_resource_units']}"
            )

        added_union = contemporary_union - native_union
        if not added_union:
            continue
        expanded_butterflies += 1
        total_added_units += len(added_union)

        for unit in added_union:
            contributing_hosts = [
                host for host in hosts if unit in added_by_host[host]
            ]
            if not contributing_hosts:
                raise RuntimeError(f"no plant contributor for {butterfly} / {unit}")
            plant_share = 1.0 / len(contributing_hosts)
            for host in contributing_hosts:
                plant_credit[host] += plant_share

            contributing_genera = sorted(
                {
                    plant_meta[host]["genus"]
                    for host in contributing_hosts
                    if plant_meta[host]["genus"]
                }
            )
            if not contributing_genera:
                raise RuntimeError(f"no genus contributor for {butterfly} / {unit}")
            genus_share = 1.0 / len(contributing_genera)
            for genus in contributing_genera:
                genus_credit[genus] += genus_share

    if expanded_butterflies != 206:
        raise RuntimeError(
            f"expected 206 expanded butterflies; got {expanded_butterflies}"
        )
    if total_added_units != 14553:
        raise RuntimeError(
            f"expected 14,553 added butterfly×WGSRPD3 units; got {total_added_units}"
        )
    if abs(sum(plant_credit.values()) - total_added_units) > 1e-8:
        raise RuntimeError("plant fractional credits do not preserve total added units")
    if abs(sum(genus_credit.values()) - total_added_units) > 1e-8:
        raise RuntimeError("genus fractional credits do not preserve total added units")

    plant_rows = []
    for plant_id, credit in sorted(
        plant_credit.items(), key=lambda kv: (-kv[1], kv[0])
    ):
        meta = plant_meta[plant_id]
        plant_rows.append(
            {
                "accepted_plant_name_id": plant_id,
                "accepted_name": meta["accepted_name"],
                "genus": meta["genus"],
                "family": meta["family"],
                "fractional_credit": credit,
                "share_of_all_added_units": credit / total_added_units,
            }
        )

    genus_family: dict[str, set[str]] = defaultdict(set)
    for meta in plant_meta.values():
        if meta["genus"] and meta["family"]:
            genus_family[meta["genus"]].add(meta["family"])
    genus_rows = []
    for genus, credit in sorted(
        genus_credit.items(), key=lambda kv: (-kv[1], kv[0])
    ):
        families = sorted(genus_family.get(genus, set()))
        genus_rows.append(
            {
                "genus": genus,
                "family": families[0] if len(families) == 1 else ";".join(families),
                "fractional_credit": credit,
                "share_of_all_added_units": credit / total_added_units,
            }
        )

    payload = {
        "schema": "chocho_butterfly_host_contribution_concentration_v0.1",
        "status": "SUCCESS_POSTHOC_HOST_CONTRIBUTION_CONCENTRATION",
        "credit_definition": (
            "Within each butterfly×WGSRPD3 unit added by introduced host distributions, "
            "credit is divided equally among contributing host species. Genus sensitivity "
            "first collapses contributing host species to unique genera within that unit, "
            "then divides one unit equally among genera. Both decompositions therefore sum "
            "exactly to the 14,553 added butterfly×WGSRPD3 units."
        ),
        "expanded_butterflies": expanded_butterflies,
        "added_butterfly_x_wgsrpd3_units": total_added_units,
        "plant_species": summarize_credit(plant_credit),
        "genera": summarize_credit(genus_credit),
        "top_plant_species": plant_rows[:20],
        "top_genera": genus_rows[:20],
        "claim_boundary": (
            "This is a descriptive decomposition of reconstructed added resource opportunity. "
            "It does not identify causal plant traits, realized larval use, or historical "
            "introduction pathways."
        ),
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with args.output_plant_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(plant_rows[0].keys()))
        writer.writeheader()
        writer.writerows(plant_rows)
    with args.output_genus_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(genus_rows[0].keys()))
        writer.writeheader()
        writer.writerows(genus_rows)

    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

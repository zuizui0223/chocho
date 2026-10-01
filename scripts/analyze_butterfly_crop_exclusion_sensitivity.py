#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

EXPECTED_DESCRIPTOR_SHA256 = (
    "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def average_ranks(values) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def spearman(a, b) -> float | None:
    if len(a) < 3 or len(a) != len(b):
        return None
    x = average_ranks(a)
    y = average_ranks(b)
    if np.std(x) <= np.sqrt(np.finfo(float).eps):
        return None
    if np.std(y) <= np.sqrt(np.finfo(float).eps):
        return None
    return float(np.corrcoef(x, y)[0, 1])


def base_binomial(name: str) -> str:
    parts = str(name).strip().split()
    if len(parts) < 2:
        return ""
    return f"{parts[0]} {parts[1]}"


def genus(name: str) -> str:
    parts = str(name).strip().split()
    return parts[0] if parts else ""


def load_descriptors(path: Path) -> dict[str, dict[str, object]]:
    if sha256_path(path) != EXPECTED_DESCRIPTOR_SHA256:
        raise RuntimeError("descriptor SHA drift")
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = str(row["species"]).strip()
            out[name] = {
                "species": name,
                "host_family_count": float(row["host_family_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
                "native_resource_units_descriptor": int(row["host_wgsrpd3_unit_count"]),
            }
    if len(out) != 339:
        raise RuntimeError(f"expected 339 descriptors, got {len(out)}")
    return out


def load_interactions(path: Path) -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    by_butterfly: dict[str, set[str]] = defaultdict(set)
    host_meta: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            insect = str(row["insect_species"]).strip()
            hid = str(row["accepted_plant_name_id"]).strip()
            if not insect or not hid:
                continue
            by_butterfly[insect].add(hid)
            meta = {
                "accepted_name": str(row.get("accepted_name") or "").strip(),
                "input_host_name": str(row.get("input_host_name") or "").strip(),
                "family": str(row.get("family") or "").strip(),
            }
            if hid in host_meta and host_meta[hid] != meta:
                # Multiple input synonyms can resolve to the same accepted host.
                old = host_meta[hid]
                if old["accepted_name"] != meta["accepted_name"]:
                    raise RuntimeError(f"accepted host identity drift for {hid}")
                if not old["input_host_name"] and meta["input_host_name"]:
                    host_meta[hid] = meta
            else:
                host_meta[hid] = meta
    return by_butterfly, host_meta


def load_units(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            hid = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if hid and unit:
                out[hid].add(unit)
    return out


def fractional_contributions(
    host_added: dict[str, set[str]],
    added_union: set[str],
) -> dict[str, object]:
    if not added_union:
        return {
            "maximum_single_host_fractional_share": None,
            "effective_contributor_number": None,
            "fractional_credits": {},
        }
    credit: dict[str, float] = defaultdict(float)
    for unit in added_union:
        contributors = [hid for hid, units in host_added.items() if unit in units]
        if not contributors:
            raise RuntimeError(f"added unit {unit} lacks contributing host")
        share = 1.0 / len(contributors)
        for hid in contributors:
            credit[hid] += share
    total = float(len(added_union))
    if abs(sum(credit.values()) - total) > 1e-9:
        raise RuntimeError("fractional credit does not sum to added union")
    p = np.asarray([v / total for v in credit.values()], dtype=float)
    return {
        "maximum_single_host_fractional_share": float(np.max(p)),
        "effective_contributor_number": float(1.0 / np.sum(p * p)),
        "fractional_credits": dict(credit),
    }


def classify_hosts(
    host_meta: dict[str, dict[str, str]],
    crop: dict,
    *,
    conservative_genus: bool,
) -> set[str]:
    exact = set(map(str, crop["exact_binomials"]))
    genera = set(map(str, crop["genus_wildcards"]))
    excluded = set()
    for hid, meta in host_meta.items():
        names = [meta.get("accepted_name", ""), meta.get("input_host_name", "")]
        if any(base_binomial(name) in exact for name in names if name):
            excluded.add(hid)
            continue
        if conservative_genus and any(genus(name) in genera for name in names if name):
            excluded.add(hid)
    return excluded


def reconstruct(
    descriptors: dict[str, dict[str, object]],
    pairs: dict[str, set[str]],
    native_by_host: dict[str, set[str]],
    contemporary_by_host: dict[str, set[str]],
    excluded_hosts: set[str],
) -> list[dict[str, object]]:
    rows = []
    for name, d in descriptors.items():
        if float(d["host_family_count"]) <= 0 or int(d["native_resource_units_descriptor"]) <= 0:
            continue
        original_hosts = set(pairs.get(name, set()))
        if not original_hosts:
            raise RuntimeError(f"resource-eligible species lacks resolved hosts: {name}")
        hosts = sorted(original_hosts - excluded_hosts)

        native_union: set[str] = set()
        contemporary_union: set[str] = set()
        added_by_host: dict[str, set[str]] = {}
        for hid in hosts:
            n = set(native_by_host.get(hid, set()))
            c = set(contemporary_by_host.get(hid, set()))
            if not n.issubset(c):
                raise RuntimeError(f"native host units not subset contemporary: {hid}")
            native_union.update(n)
            contemporary_union.update(c)
            added_by_host[hid] = c - n

        added_union = contemporary_union - native_union
        contribution = fractional_contributions(added_by_host, added_union)
        n_native = len(native_union)
        n_contemporary = len(contemporary_union)
        ratio = (n_contemporary / n_native) if n_native else None
        rows.append({
            "species": name,
            "host_family_count": float(d["host_family_count"]),
            "resolved_host_species": int(d["resolved_host_species"]),
            "host_taxonomy_lower_bound_adequate": (
                int(d["resolved_host_species"]) >= float(d["host_family_count"])
            ),
            "original_host_species": len(original_hosts),
            "remaining_host_species": len(hosts),
            "crop_hosts_removed": len(original_hosts & excluded_hosts),
            "native_resource_units": n_native,
            "contemporary_resource_units": n_contemporary,
            "introduced_added_units": len(added_union),
            "log_resource_expansion": math.log(ratio) if ratio and ratio > 0 else None,
            "maximum_single_host_fractional_share": contribution[
                "maximum_single_host_fractional_share"
            ],
            "effective_contributor_number": contribution[
                "effective_contributor_number"
            ],
            "_fractional_credits": contribution["fractional_credits"],
        })
    if len(rows) != 239:
        raise RuntimeError(f"expected 239 original resource-eligible species, got {len(rows)}")
    return rows


def concentration_half_count(rows: list[dict[str, object]]) -> dict[str, object]:
    aggregate: dict[str, float] = defaultdict(float)
    total = 0.0
    for row in rows:
        for hid, value in row["_fractional_credits"].items():
            aggregate[str(hid)] += float(value)
            total += float(value)
    ordered = sorted(aggregate.items(), key=lambda kv: (-kv[1], kv[0]))
    if total <= 0:
        return {
            "total_fractional_added_credit": 0.0,
            "contributing_host_species": 0,
            "hosts_for_half_added_credit": 0,
        }
    cumulative = 0.0
    half_count = 0
    for _, value in ordered:
        cumulative += value
        half_count += 1
        if cumulative >= 0.5 * total:
            break
    return {
        "total_fractional_added_credit": float(total),
        "contributing_host_species": len(ordered),
        "hosts_for_half_added_credit": int(half_count),
        "top_host_fractional_credit": float(ordered[0][1]),
        "top_host_share_of_all_added_credit": float(ordered[0][1] / total),
    }


def summarize(rows: list[dict[str, object]], *, adequate_only: bool) -> dict[str, object]:
    population = [
        row for row in rows
        if (not adequate_only or bool(row["host_taxonomy_lower_bound_adequate"]))
    ]
    evaluable = [row for row in population if int(row["native_resource_units"]) > 0]
    expanded = [row for row in evaluable if int(row["introduced_added_units"]) > 0]
    total_native = sum(int(row["native_resource_units"]) for row in population)
    total_contemporary = sum(int(row["contemporary_resource_units"]) for row in population)
    total_added = total_contemporary - total_native

    x = [float(row["host_family_count"]) for row in evaluable]
    expansion = [float(row["log_resource_expansion"]) for row in evaluable]

    arch = [
        row for row in expanded
        if row["effective_contributor_number"] is not None
        and row["maximum_single_host_fractional_share"] is not None
    ]
    ax = [float(row["host_family_count"]) for row in arch]
    eff = [float(row["effective_contributor_number"]) for row in arch]
    dom = [float(row["maximum_single_host_fractional_share"]) for row in arch]

    return {
        "original_population_species": len(population),
        "species_with_non_crop_native_resource": len(evaluable),
        "species_with_no_non_crop_native_resource": len(population) - len(evaluable),
        "species_expanded": len(expanded),
        "fraction_evaluable_species_expanded": (
            len(expanded) / len(evaluable) if evaluable else None
        ),
        "total_native_species_units": total_native,
        "total_contemporary_species_units": total_contemporary,
        "total_added_species_units": total_added,
        "aggregate_proportional_increase": (
            total_added / total_native if total_native else None
        ),
        "spearman_host_family_vs_log_resource_expansion": spearman(x, expansion),
        "portfolio_expanded_species": len(arch),
        "median_maximum_single_host_fractional_share": (
            float(np.median(dom)) if dom else None
        ),
        "median_effective_contributor_number": (
            float(np.median(eff)) if eff else None
        ),
        "spearman_host_family_vs_effective_contributor_number": spearman(ax, eff),
        "spearman_host_family_vs_maximum_single_host_fractional_share": spearman(ax, dom),
        "aggregate_host_contribution_concentration": concentration_half_count(expanded),
    }


def clean_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        {k: v for k, v in row.items() if k != "_fractional_credits"}
        for row in rows
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--crop-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "chocho_butterfly_crop_exclusion_sensitivity_v0.1":
        raise RuntimeError("unexpected crop sensitivity protocol")
    if protocol.get("status") != "FROZEN_BEFORE_CROP_EXCLUSION_RESULT":
        raise RuntimeError("crop sensitivity protocol is not frozen")

    crop = json.loads(args.crop_json.read_text(encoding="utf-8"))
    if crop.get("schema") != "chocho_fao_crop_taxa_v0.1":
        raise RuntimeError("unexpected crop taxonomy artifact")

    descriptors = load_descriptors(args.descriptors_csv)
    pairs, host_meta = load_interactions(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)

    baseline = reconstruct(descriptors, pairs, native, contemporary, set())
    baseline_full = summarize(baseline, adequate_only=False)
    baseline_adequate = summarize(baseline, adequate_only=True)
    expected = {
        "species_expanded": 206,
        "total_native_species_units": 26530,
        "total_contemporary_species_units": 41083,
    }
    for key, value in expected.items():
        if baseline_full[key] != value:
            raise RuntimeError(
                f"baseline reconstruction drift for {key}: "
                f"{baseline_full[key]} != {value}"
            )

    variants = {}
    csv_rows = []
    for label, conservative in (
        ("fao_exact_binomial", False),
        ("fao_exact_plus_spp_genus", True),
    ):
        excluded = classify_hosts(
            host_meta,
            crop,
            conservative_genus=conservative,
        )
        rows = reconstruct(descriptors, pairs, native, contemporary, excluded)
        variants[label] = {
            "crop_host_ids_excluded": len(excluded),
            "butterflies_with_at_least_one_crop_host_removed": sum(
                int(row["crop_hosts_removed"]) > 0 for row in rows
            ),
            "full_resource_eligible_panel": summarize(rows, adequate_only=False),
            "host_taxonomy_lower_bound_adequate_subset": summarize(
                rows, adequate_only=True
            ),
        }
        for row in clean_rows(rows):
            csv_rows.append({"variant": label, **row})

    payload = {
        "schema": "chocho_butterfly_crop_exclusion_sensitivity_result_v0.1",
        "status": "CROP_EXCLUSION_SENSITIVITY_COMPLETE",
        "protocol_sha256": sha256_path(args.protocol_json),
        "crop_taxa_sha256": sha256_path(args.crop_json),
        "crop_source_html_sha256": crop["source_html_sha256"],
        "baseline_reconstruction": {
            "full_resource_eligible_panel": baseline_full,
            "host_taxonomy_lower_bound_adequate_subset": baseline_adequate,
        },
        "variants": variants,
        "interpretation_boundary": {
            "primary_analysis_changed": False,
            "specialization_predictor_changed": False,
            "crop_exclusion_is_sensitivity_only": True,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    fields = list(csv_rows[0].keys())
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(csv_rows)

    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

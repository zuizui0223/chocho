#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_DESCRIPTOR_SHA256 = (
    "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
)
EXPECTED_RESOURCE_SPECIES = 239
EXPECTED_ADEQUATE_SPECIES = 215


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


def clean_text(value: object) -> str:
    text = str(value or "")
    text = text.replace("\u00d7", " ").replace("\u00a0", " ")
    return re.sub(r"\s+", " ", text).strip()


def binomial(name: str) -> str | None:
    text = clean_text(name)
    m = re.search(r"\b([A-Z][A-Za-z-]+)\s+([a-z][A-Za-z-]+)\b", text)
    if not m:
        return None
    genus, epithet = m.groups()
    if epithet.lower() in {"sp", "spp", "species", "of", "and"}:
        return None
    return f"{genus} {epithet}"


def genus(name: str) -> str | None:
    text = clean_text(name)
    m = re.search(r"\b([A-Z][A-Za-z-]+)\b", text)
    return None if m is None else m.group(1)


def parse_fao_crop_taxa(path: Path) -> dict[str, object]:
    tables = pd.read_html(path)
    crop = None
    for table in tables:
        cols = [clean_text(x) for x in table.columns]
        if any("Botanical name" in x for x in cols) and any("Crop name" in x for x in cols):
            table = table.copy()
            table.columns = cols
            crop = table
            break
    if crop is None:
        raise RuntimeError("could not find FAO crop botanical-name table")

    botanical_col = next(x for x in crop.columns if "Botanical name" in x)
    exact_species: set[str] = set()
    genus_level: set[str] = set()
    source_rows = 0

    for raw in crop[botanical_col].dropna().tolist():
        text = clean_text(raw)
        if not text:
            continue
        source_rows += 1
        parts = [clean_text(x) for x in re.split(r";", text) if clean_text(x)]
        for part in parts:
            # Genus-level entries are deliberately reserved for the conservative tier.
            for gm in re.finditer(r"\b([A-Z][A-Za-z-]+)\s+spp?\.", part):
                genus_level.add(gm.group(1))
            # Extract explicit binomials even from entries that also contain qualifiers.
            for m in re.finditer(r"\b([A-Z][A-Za-z-]+)\s+([a-z][A-Za-z-]+)\b", part):
                g, e = m.groups()
                if e.lower() in {"sp", "spp", "species", "of", "and"}:
                    continue
                exact_species.add(f"{g} {e}")

    if "Zea mays" not in exact_species or "Medicago sativa" not in exact_species:
        raise RuntimeError("FAO crop parsing failed canonical sentinel taxa")
    if "Avena" not in genus_level:
        raise RuntimeError("FAO crop parsing failed Avena genus-level sentinel")

    return {
        "source_table_rows_with_botanical_name": int(source_rows),
        "strict_species": sorted(exact_species),
        "conservative_genera": sorted(genus_level),
    }


def load_descriptors(path: Path) -> dict[str, dict[str, object]]:
    if sha256_path(path) != EXPECTED_DESCRIPTOR_SHA256:
        raise RuntimeError("descriptor SHA drift")
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species",
            "host_family_count",
            "host_wgsrpd3_unit_count",
            "resolved_host_species",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in reader:
            name = clean_text(row["species"])
            out[name] = {
                "species": name,
                "host_family_count": float(row["host_family_count"]),
                "native_descriptor_units": int(row["host_wgsrpd3_unit_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
            }
    if len(out) != 339:
        raise RuntimeError(f"expected 339 descriptor species; got {len(out)}")
    return out


def load_interactions(path: Path) -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    by_butterfly: dict[str, set[str]] = defaultdict(set)
    host_meta: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "insect_species",
            "accepted_plant_name_id",
            "input_host_name",
            "accepted_name",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for row in reader:
            insect = clean_text(row["insect_species"])
            hid = clean_text(row["accepted_plant_name_id"])
            if not insect or not hid:
                continue
            by_butterfly[insect].add(hid)
            meta = host_meta.setdefault(
                hid,
                {
                    "accepted_plant_name_id": hid,
                    "accepted_name": clean_text(row.get("accepted_name", "")),
                    "input_host_name": clean_text(row.get("input_host_name", "")),
                },
            )
            # Multiple HOSTS input names can resolve to one accepted WCVP id.
            incoming = clean_text(row.get("input_host_name", ""))
            if incoming and incoming not in meta["input_host_name"].split(" | "):
                meta["input_host_name"] = " | ".join(
                    x for x in (meta["input_host_name"], incoming) if x
                )
    return by_butterfly, host_meta


def load_units(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("distribution sidecar schema drift")
        for row in reader:
            hid = clean_text(row["accepted_plant_name_id"])
            unit = clean_text(row["area_code_l3"])
            if hid and unit:
                out[hid].add(unit)
    return out


def crop_host_sets(
    host_meta: dict[str, dict[str, str]],
    crop_taxa: dict[str, object],
) -> tuple[set[str], set[str], list[dict[str, object]]]:
    strict_names = set(map(str, crop_taxa["strict_species"]))
    broad_genera = set(map(str, crop_taxa["conservative_genera"]))
    strict: set[str] = set()
    conservative: set[str] = set()
    rows = []

    for hid, meta in sorted(host_meta.items()):
        names = [meta.get("accepted_name", "")]
        names.extend(x.strip() for x in meta.get("input_host_name", "").split(" | "))
        bins = sorted({x for x in (binomial(name) for name in names) if x})
        gens = sorted({x for x in (genus(name) for name in names) if x})
        is_strict = any(name in strict_names for name in bins)
        is_broad = is_strict or any(g in broad_genera for g in gens)
        if is_strict:
            strict.add(hid)
        if is_broad:
            conservative.add(hid)
        if is_broad:
            rows.append(
                {
                    "accepted_plant_name_id": hid,
                    "accepted_name": meta.get("accepted_name", ""),
                    "input_host_name": meta.get("input_host_name", ""),
                    "strict_crop": bool(is_strict),
                    "conservative_crop": bool(is_broad),
                    "matched_binomials": bins,
                    "matched_genera": gens,
                }
            )
    return strict, conservative, rows


def fractional_host_contributions(
    host_added_units: dict[str, set[str]],
    butterfly_added_union: set[str],
) -> dict[str, object]:
    if not butterfly_added_union:
        return {
            "introduced_added_units": 0,
            "contributing_host_species": 0,
            "maximum_single_host_fractional_share": None,
            "effective_contributor_number": None,
            "fractional_credits": {},
        }

    credit: dict[str, float] = defaultdict(float)
    for unit in sorted(butterfly_added_union):
        contributors = [
            host for host, units in host_added_units.items() if unit in units
        ]
        if not contributors:
            raise RuntimeError(f"added unit lacks a host contributor: {unit}")
        share = 1.0 / len(contributors)
        for host in contributors:
            credit[host] += share

    total = float(len(butterfly_added_union))
    if abs(sum(credit.values()) - total) > 1e-9:
        raise RuntimeError("fractional credits do not sum to added union")
    props = sorted((value / total for value in credit.values()), reverse=True)
    return {
        "introduced_added_units": int(total),
        "contributing_host_species": len(props),
        "maximum_single_host_fractional_share": float(props[0]),
        "effective_contributor_number": float(
            1.0 / sum(value * value for value in props)
        ),
        "fractional_credits": dict(credit),
    }


def host_concentration(
    pooled_credits: dict[str, float],
    host_meta: dict[str, dict[str, str]],
) -> dict[str, object]:
    ordered = sorted(pooled_credits.items(), key=lambda kv: (-kv[1], kv[0]))
    total = float(sum(value for _, value in ordered))
    cumulative = 0.0
    half = 0
    for i, (_, value) in enumerate(ordered, start=1):
        cumulative += value
        if half == 0 and cumulative >= 0.5 * total:
            half = i
    top = []
    for hid, value in ordered[:20]:
        top.append(
            {
                "accepted_plant_name_id": hid,
                "accepted_name": host_meta.get(hid, {}).get("accepted_name", ""),
                "fractional_added_units": float(value),
                "share": None if total == 0 else float(value / total),
            }
        )
    return {
        "contributing_hosts": len(ordered),
        "pooled_fractional_added_units": total,
        "hosts_to_50_percent": int(half),
        "top_20": top,
    }


def analyze_mode(
    mode: str,
    excluded_hosts: set[str],
    descriptors: dict[str, dict[str, object]],
    interactions: dict[str, set[str]],
    native_by_host: dict[str, set[str]],
    contemporary_by_host: dict[str, set[str]],
    host_meta: dict[str, dict[str, str]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows = []
    pooled_credit: dict[str, float] = defaultdict(float)

    for species, descriptor in sorted(descriptors.items()):
        host_family_count = float(descriptor["host_family_count"])
        if host_family_count <= 0 or int(descriptor["native_descriptor_units"]) <= 0:
            continue
        original_hosts = sorted(interactions.get(species, set()))
        if not original_hosts:
            raise RuntimeError(f"resource-eligible butterfly lacks resolved hosts: {species}")
        hosts = [h for h in original_hosts if h not in excluded_hosts]

        native_union: set[str] = set()
        contemporary_union: set[str] = set()
        added_by_host: dict[str, set[str]] = {}
        for hid in hosts:
            native = set(native_by_host.get(hid, set()))
            contemporary = set(contemporary_by_host.get(hid, set()))
            if not native.issubset(contemporary):
                raise RuntimeError(f"native not subset contemporary for host {hid}")
            native_union.update(native)
            contemporary_union.update(contemporary)
            added_by_host[hid] = contemporary - native

        added_union = contemporary_union - native_union
        contribution = fractional_host_contributions(added_by_host, added_union)
        for hid, value in contribution["fractional_credits"].items():
            pooled_credit[hid] += float(value)

        n_native = len(native_union)
        n_contemporary = len(contemporary_union)
        resolved_hosts = int(descriptor["resolved_host_species"])
        rows.append(
            {
                "species": species,
                "mode": mode,
                "host_family_count": host_family_count,
                "resolved_host_species_original": resolved_hosts,
                "host_taxonomy_lower_bound_adequate_original": (
                    resolved_hosts >= host_family_count
                ),
                "original_resolved_hosts_in_sidecar": len(original_hosts),
                "remaining_hosts_after_crop_exclusion": len(hosts),
                "excluded_crop_hosts": len(original_hosts) - len(hosts),
                "native_resource_units": n_native,
                "contemporary_resource_units": n_contemporary,
                "introduced_added_units": len(added_union),
                "log_resource_expansion": (
                    None
                    if n_native == 0
                    else math.log1p(n_contemporary) - math.log1p(n_native)
                ),
                "contributing_host_species": contribution["contributing_host_species"],
                "maximum_single_host_fractional_share": contribution[
                    "maximum_single_host_fractional_share"
                ],
                "effective_contributor_number": contribution[
                    "effective_contributor_number"
                ],
            }
        )

    if len(rows) != EXPECTED_RESOURCE_SPECIES:
        raise RuntimeError(f"{mode}: expected 239 resource species; got {len(rows)}")

    retained = [row for row in rows if int(row["native_resource_units"]) > 0]
    x = [float(row["host_family_count"]) for row in retained]
    y = [float(row["log_resource_expansion"]) for row in retained]

    adequate = [
        row
        for row in rows
        if bool(row["host_taxonomy_lower_bound_adequate_original"])
    ]
    if len(adequate) != EXPECTED_ADEQUATE_SPECIES:
        raise RuntimeError(f"{mode}: expected 215 original adequate species; got {len(adequate)}")
    adequate_expanded = [
        row
        for row in adequate
        if int(row["introduced_added_units"]) > 0
        and row["effective_contributor_number"] is not None
    ]

    summary = {
        "mode": mode,
        "excluded_host_taxa": len(excluded_hosts),
        "resource_eligible_species_original": len(rows),
        "species_retaining_nonzero_native_resource": len(retained),
        "species_with_zero_non_crop_resource": len(rows) - len(retained),
        "species_expanded_of_original_239": int(
            sum(int(row["introduced_added_units"]) > 0 for row in rows)
        ),
        "fraction_species_expanded_of_original_239": float(
            sum(int(row["introduced_added_units"]) > 0 for row in rows)
            / len(rows)
        ),
        "aggregate": {
            "native_species_units": int(sum(row["native_resource_units"] for row in rows)),
            "contemporary_species_units": int(
                sum(row["contemporary_resource_units"] for row in rows)
            ),
            "added_species_units": int(sum(row["introduced_added_units"] for row in rows)),
        },
        "host_family_vs_log_resource_expansion": {
            "n": len(retained),
            "rho": spearman(x, y),
        },
        "original_adequate_architecture": {
            "adequate_species_original": len(adequate),
            "expanded_species_after_exclusion": len(adequate_expanded),
            "host_family_vs_effective_contributor_number": spearman(
                [float(row["host_family_count"]) for row in adequate_expanded],
                [float(row["effective_contributor_number"]) for row in adequate_expanded],
            ),
            "host_family_vs_maximum_single_host_share": spearman(
                [float(row["host_family_count"]) for row in adequate_expanded],
                [
                    float(row["maximum_single_host_fractional_share"])
                    for row in adequate_expanded
                ],
            ),
            "median_effective_contributor_number": (
                None
                if not adequate_expanded
                else float(
                    np.median(
                        [
                            float(row["effective_contributor_number"])
                            for row in adequate_expanded
                        ]
                    )
                )
            ),
            "median_maximum_single_host_fractional_share": (
                None
                if not adequate_expanded
                else float(
                    np.median(
                        [
                            float(row["maximum_single_host_fractional_share"])
                            for row in adequate_expanded
                        ]
                    )
                )
            ),
        },
        "pooled_host_contribution_concentration": host_concentration(
            pooled_credit, host_meta
        ),
    }
    native_total = summary["aggregate"]["native_species_units"]
    added_total = summary["aggregate"]["added_species_units"]
    summary["aggregate"]["percent_increase"] = (
        None if native_total == 0 else float(100.0 * added_total / native_total)
    )
    return summary, rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--fao-html", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-species-csv", type=Path, required=True)
    ap.add_argument("--output-crop-hosts-csv", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "chocho_crop_exclusion_sensitivity_protocol_v0.1":
        raise RuntimeError("unexpected crop sensitivity protocol")
    if protocol.get("status") != "FROZEN_BEFORE_CROP_EXCLUSION_RESULT":
        raise RuntimeError("crop sensitivity protocol not frozen")

    descriptors = load_descriptors(args.descriptors_csv)
    interactions, host_meta = load_interactions(args.insect_host_csv)
    native_by_host = load_units(args.native_distribution_csv)
    contemporary_by_host = load_units(args.contemporary_distribution_csv)
    crop_taxa = parse_fao_crop_taxa(args.fao_html)
    strict_hosts, conservative_hosts, crop_host_rows = crop_host_sets(
        host_meta, crop_taxa
    )

    baseline, baseline_rows = analyze_mode(
        "baseline",
        set(),
        descriptors,
        interactions,
        native_by_host,
        contemporary_by_host,
        host_meta,
    )
    strict, strict_rows = analyze_mode(
        "strict",
        strict_hosts,
        descriptors,
        interactions,
        native_by_host,
        contemporary_by_host,
        host_meta,
    )
    conservative, conservative_rows = analyze_mode(
        "conservative",
        conservative_hosts,
        descriptors,
        interactions,
        native_by_host,
        contemporary_by_host,
        host_meta,
    )

    # Bind this sensitivity to the already frozen manuscript values.
    if baseline["aggregate"]["native_species_units"] != 26530:
        raise RuntimeError("baseline aggregate native resource units drift")
    if baseline["aggregate"]["contemporary_species_units"] != 41083:
        raise RuntimeError("baseline aggregate contemporary resource units drift")
    if baseline["species_expanded_of_original_239"] != 206:
        raise RuntimeError("baseline expanded-species count drift")
    if abs(float(baseline["host_family_vs_log_resource_expansion"]["rho"]) - 0.008) > 0.002:
        raise RuntimeError("baseline host-breadth expansion correlation drift")
    if baseline["pooled_host_contribution_concentration"]["hosts_to_50_percent"] != 38:
        raise RuntimeError("baseline host-contribution concentration drift")

    payload = {
        "schema": "chocho_crop_exclusion_sensitivity_result_v0.1",
        "status": "CROP_EXCLUSION_SENSITIVITY_COMPLETE",
        "protocol_sha256": sha256_path(args.protocol_json),
        "fao_source": {
            "url": protocol["crop_reference"]["source_url"],
            "html_sha256": sha256_path(args.fao_html),
            "crop_table_rows": crop_taxa["source_table_rows_with_botanical_name"],
            "strict_species_names": len(crop_taxa["strict_species"]),
            "genus_level_crop_entries": len(crop_taxa["conservative_genera"]),
        },
        "classified_hosts": {
            "strict_crop_hosts": len(strict_hosts),
            "conservative_crop_hosts": len(conservative_hosts),
        },
        "baseline": baseline,
        "strict": strict,
        "conservative": conservative,
        "claim_boundary": {
            "sensitivity_only": True,
            "host_family_breadth_not_redefined": True,
            "butterfly_panel_not_reselected": True,
            "wcvp_introduced_wild_distribution_only": True,
        },
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    all_rows = baseline_rows + strict_rows + conservative_rows
    fields = list(all_rows[0].keys())
    with args.output_species_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)

    crop_fields = [
        "accepted_plant_name_id",
        "accepted_name",
        "input_host_name",
        "strict_crop",
        "conservative_crop",
        "matched_binomials",
        "matched_genera",
    ]
    with args.output_crop_hosts_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=crop_fields)
        writer.writeheader()
        for row in crop_host_rows:
            row = dict(row)
            row["matched_binomials"] = ";".join(row["matched_binomials"])
            row["matched_genera"] = ";".join(row["matched_genera"])
            writer.writerow(row)

    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

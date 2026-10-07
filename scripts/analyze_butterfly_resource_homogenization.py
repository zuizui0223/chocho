#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import median


def load_descriptors(path: Path) -> dict[str, dict[str, float | int]]:
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"species", "host_family_count", "host_wgsrpd3_unit_count"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("descriptor schema drift")
        for row in reader:
            species = str(row["species"]).strip()
            out[species] = {
                "host_family_count": float(row["host_family_count"]),
                "native_resource_units": int(row["host_wgsrpd3_unit_count"]),
            }
    return out


def load_taxonomy(path: Path) -> dict[str, str]:
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"Species", "Family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("LepTraits taxonomy columns missing")
        for row in reader:
            species = str(row["Species"]).strip()
            family = str(row["Family"]).strip()
            if species and family:
                out[species] = family
    return out


def load_pairs(path: Path):
    by_butterfly: dict[str, set[str]] = defaultdict(set)
    plant_meta: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id", "accepted_name", "family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for row in reader:
            butterfly = str(row["insect_species"]).strip()
            plant_id = str(row["accepted_plant_name_id"]).strip()
            if not butterfly or not plant_id:
                continue
            by_butterfly[butterfly].add(plant_id)
            meta = {
                "accepted_name": str(row["accepted_name"]).strip(),
                "family": str(row["family"]).strip(),
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


def mean_pairwise_jaccard(masks: list[int]) -> tuple[float, float, int]:
    values = []
    for i in range(len(masks)):
        a = masks[i]
        for j in range(i + 1, len(masks)):
            b = masks[j]
            union = (a | b).bit_count()
            if union == 0:
                continue
            values.append((a & b).bit_count() / union)
    if not values:
        raise RuntimeError("no region pairs with non-empty union")
    return sum(values) / len(values), median(values), len(values)


def mean_pairwise_set_jaccard(sets: list[set[int]]) -> tuple[float, float, int]:
    values = []
    for i in range(len(sets)):
        a = sets[i]
        for j in range(i + 1, len(sets)):
            b = sets[j]
            union = len(a | b)
            if union == 0:
                continue
            values.append(len(a & b) / union)
    if not values:
        raise RuntimeError("no set pairs with non-empty union")
    return sum(values) / len(values), median(values), len(values)


def masks_from_species_units(
    species_order: list[str],
    region_order: list[str],
    units_by_species: dict[str, set[str]],
) -> list[int]:
    region_index = {unit: j for j, unit in enumerate(region_order)}
    masks = [0 for _ in region_order]
    for i, species in enumerate(species_order):
        bit = 1 << i
        for unit in units_by_species.get(species, set()):
            j = region_index.get(unit)
            if j is not None:
                masks[j] |= bit
    return masks


def quantile(values: list[float], p: float) -> float:
    xs = sorted(values)
    if not xs:
        raise RuntimeError("empty quantile")
    pos = p * (len(xs) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    frac = pos - lo
    return xs[lo] * (1 - frac) + xs[hi] * frac


def swap_null(
    native_by_species: list[set[int]],
    added_by_species: list[set[int]],
    native_masks: list[int],
    permutations: int,
    seed: int,
) -> dict[str, object]:
    rng = random.Random(seed)
    rows = [set(x) for x in added_by_species]
    edges = [(i, col) for i, cols in enumerate(rows) for col in cols]
    if not edges:
        raise RuntimeError("no added edges for null")

    added_masks = [0 for _ in native_masks]
    for i, cols in enumerate(rows):
        bit = 1 << i
        for col in cols:
            added_masks[col] |= bit

    def attempt_swaps(attempts: int) -> int:
        accepted = 0
        for _ in range(attempts):
            k1 = rng.randrange(len(edges))
            k2 = rng.randrange(len(edges))
            if k1 == k2:
                continue
            i, a = edges[k1]
            j, b = edges[k2]
            if i == j or a == b:
                continue
            if b in rows[i] or a in rows[j]:
                continue
            if b in native_by_species[i] or a in native_by_species[j]:
                continue

            rows[i].remove(a)
            rows[j].remove(b)
            rows[i].add(b)
            rows[j].add(a)
            edges[k1] = (i, b)
            edges[k2] = (j, a)

            bi = 1 << i
            bj = 1 << j
            added_masks[a] ^= bi
            added_masks[a] ^= bj
            added_masks[b] ^= bj
            added_masks[b] ^= bi
            accepted += 1
        return accepted

    burn_attempts = max(10000, 10 * len(edges))
    burn_accepted = attempt_swaps(burn_attempts)
    spacing_attempts = max(2000, 2 * len(edges))

    null_means = []
    null_species_means = []
    accepted_total = burn_accepted
    for _ in range(permutations):
        accepted_total += attempt_swaps(spacing_attempts)
        masks = [native_masks[j] | added_masks[j] for j in range(len(native_masks))]
        mean_sim, _, _ = mean_pairwise_jaccard(masks)
        null_means.append(mean_sim)
        contemporary_species_sets = [
            native_by_species[i] | rows[i] for i in range(len(rows))
        ]
        species_mean, _, _ = mean_pairwise_set_jaccard(contemporary_species_sets)
        null_species_means.append(species_mean)

    return {
        "permutations": permutations,
        "seed": seed,
        "edges": len(edges),
        "burn_attempts": burn_attempts,
        "spacing_attempts": spacing_attempts,
        "accepted_swaps_total": accepted_total,
        "mean_similarity_null_median": median(null_means),
        "mean_similarity_null_ci95": [quantile(null_means, 0.025), quantile(null_means, 0.975)],
        "mean_butterfly_niche_overlap_null_median": median(null_species_means),
        "mean_butterfly_niche_overlap_null_ci95": [
            quantile(null_species_means, 0.025),
            quantile(null_species_means, 0.975),
        ],
        "_null_means": null_means,
        "_null_species_means": null_species_means,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=499)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-host-csv", type=Path, required=True)
    args = ap.parse_args()

    descriptors = load_descriptors(args.descriptors_csv)
    taxonomy = load_taxonomy(args.leptraits_csv)
    pairs, plant_meta = load_pairs(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)

    species_order = [
        sp for sp, d in descriptors.items()
        if float(d["host_family_count"]) > 0 and int(d["native_resource_units"]) > 0
    ]
    if len(species_order) != 239:
        raise RuntimeError(f"expected 239 resource-eligible butterflies; got {len(species_order)}")
    for sp in species_order:
        if sp not in taxonomy:
            raise RuntimeError(f"missing butterfly family for {sp}")

    native_by_species_name: dict[str, set[str]] = {}
    contemporary_by_species_name: dict[str, set[str]] = {}
    added_by_species_name: dict[str, set[str]] = {}
    host_credit: dict[str, float] = defaultdict(float)
    host_butterflies: dict[str, set[str]] = defaultdict(set)
    host_butterfly_families: dict[str, set[str]] = defaultdict(set)
    host_family_credit: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))

    total_added = 0
    for sp in species_order:
        hosts = sorted(pairs.get(sp, set()))
        if not hosts:
            raise RuntimeError(f"eligible butterfly lacks resolved hosts: {sp}")

        n_union: set[str] = set()
        c_union: set[str] = set()
        added_by_host: dict[str, set[str]] = {}
        for host in hosts:
            n = set(native.get(host, set()))
            c = set(contemporary.get(host, set()))
            if not n <= c:
                raise RuntimeError(f"native not subset contemporary for {host}")
            n_union.update(n)
            c_union.update(c)
            added_by_host[host] = c - n

        if len(n_union) != int(descriptors[sp]["native_resource_units"]):
            raise RuntimeError(
                f"native union drift for {sp}: {len(n_union)} vs "
                f"{descriptors[sp]['native_resource_units']}"
            )

        added_union = c_union - n_union
        native_by_species_name[sp] = n_union
        contemporary_by_species_name[sp] = c_union
        added_by_species_name[sp] = added_union
        total_added += len(added_union)

        family = taxonomy[sp]
        for unit in added_union:
            contributors = [host for host in hosts if unit in added_by_host[host]]
            if not contributors:
                raise RuntimeError(f"no contributor for {sp} / {unit}")
            share = 1.0 / len(contributors)
            for host in contributors:
                host_credit[host] += share
                host_butterflies[host].add(sp)
                host_butterfly_families[host].add(family)
                host_family_credit[host][family] += share

    if total_added != 14553:
        raise RuntimeError(f"expected 14,553 added units; got {total_added}")

    all_regions = sorted(set().union(*native_by_species_name.values(), *contemporary_by_species_name.values()))
    active_regions = sorted(
        unit for unit in all_regions
        if any(unit in native_by_species_name[sp] for sp in species_order)
    )
    contemporary_only_regions = sorted(set(all_regions) - set(active_regions))

    native_masks = masks_from_species_units(species_order, active_regions, native_by_species_name)
    contemporary_masks = masks_from_species_units(species_order, active_regions, contemporary_by_species_name)
    native_mean, native_median, pair_count = mean_pairwise_jaccard(native_masks)
    contemporary_mean, contemporary_median, pair_count2 = mean_pairwise_jaccard(contemporary_masks)
    if pair_count != pair_count2:
        raise RuntimeError("region pair count drift")

    region_index = {unit: j for j, unit in enumerate(active_regions)}
    native_by_species_idx = [
        {region_index[u] for u in native_by_species_name[sp] if u in region_index}
        for sp in species_order
    ]
    added_by_species_idx = [
        {region_index[u] for u in added_by_species_name[sp] if u in region_index}
        for sp in species_order
    ]
    contemporary_by_species_idx = [
        native_by_species_idx[i] | added_by_species_idx[i]
        for i in range(len(species_order))
    ]
    native_species_mean, native_species_median, species_pair_count = (
        mean_pairwise_set_jaccard(native_by_species_idx)
    )
    contemporary_species_mean, contemporary_species_median, species_pair_count2 = (
        mean_pairwise_set_jaccard(contemporary_by_species_idx)
    )
    if species_pair_count != species_pair_count2:
        raise RuntimeError("butterfly pair count drift")


    null = swap_null(
        native_by_species_idx,
        added_by_species_idx,
        native_masks,
        args.permutations,
        args.seed,
    )
    null_means = list(null.pop("_null_means"))
    null_species_means = list(null.pop("_null_species_means"))
    p_greater = (1 + sum(x >= contemporary_mean for x in null_means)) / (len(null_means) + 1)
    p_species_greater = (
        1 + sum(x >= contemporary_species_mean for x in null_species_means)
    ) / (len(null_species_means) + 1)

    ranked_hosts = sorted(host_credit, key=lambda h: (-host_credit[h], h))
    total_credit = sum(host_credit.values())
    cumulative = 0.0
    half_hosts = []
    for host in ranked_hosts:
        half_hosts.append(host)
        cumulative += host_credit[host]
        if cumulative >= 0.5 * total_credit:
            break
    if len(half_hosts) != 38:
        raise RuntimeError(f"expected 38 half-credit hosts; got {len(half_hosts)}")

    top_added_by_species: dict[str, set[str]] = {}
    remainder_added_by_species: dict[str, set[str]] = {}
    half_set = set(half_hosts)
    for sp in species_order:
        hosts = sorted(pairs.get(sp, set()))
        n_union = native_by_species_name[sp]
        top_added: set[str] = set()
        rem_added: set[str] = set()
        for host in hosts:
            h_added = set(contemporary.get(host, set())) - set(native.get(host, set()))
            if host in half_set:
                top_added.update(h_added)
            else:
                rem_added.update(h_added)
        top_added_by_species[sp] = top_added - n_union
        remainder_added_by_species[sp] = rem_added - n_union

    top_contemporary = {
        sp: native_by_species_name[sp] | top_added_by_species[sp]
        for sp in species_order
    }
    rem_contemporary = {
        sp: native_by_species_name[sp] | remainder_added_by_species[sp]
        for sp in species_order
    }
    top_masks = masks_from_species_units(species_order, active_regions, top_contemporary)
    rem_masks = masks_from_species_units(species_order, active_regions, rem_contemporary)
    top_mean, top_median, _ = mean_pairwise_jaccard(top_masks)
    rem_mean, rem_median, _ = mean_pairwise_jaccard(rem_masks)

    host_rows = []
    for rank, host in enumerate(ranked_hosts, start=1):
        fam_credits = host_family_credit[host]
        total_h = host_credit[host]
        hhi = sum((v / total_h) ** 2 for v in fam_credits.values()) if total_h else 0.0
        host_rows.append({
            "rank": rank,
            "accepted_plant_name_id": host,
            "accepted_name": plant_meta[host]["accepted_name"],
            "plant_family": plant_meta[host]["family"],
            "fractional_credit": total_h,
            "share_of_all_added_units": total_h / total_credit,
            "contributing_butterfly_species": len(host_butterflies[host]),
            "contributing_butterfly_families": len(host_butterfly_families[host]),
            "butterfly_families": ";".join(sorted(host_butterfly_families[host])),
            "effective_butterfly_family_number": None if hhi == 0 else 1.0 / hhi,
            "in_half_credit_set": host in half_set,
        })

    top38_rows = [r for r in host_rows if r["in_half_credit_set"]]
    top38_credit = sum(float(r["fractional_credit"]) for r in top38_rows)
    credit_multifamily = sum(
        float(r["fractional_credit"]) for r in top38_rows
        if int(r["contributing_butterfly_families"]) >= 2
    )

    payload = {
        "schema": "chocho_butterfly_resource_homogenization_v0.1",
        "status": "SUCCESS_POSTHOC_RESOURCE_ASSEMBLAGE_HOMOGENIZATION",
        "question": (
            "Does anthropogenic host redistribution make regional butterfly larval-resource "
            "assemblages more compositionally similar, beyond the homogenization expected "
            "from the same butterfly-specific and region-specific amount of added opportunity?"
        ),
        "panel": {
            "butterflies": len(species_order),
            "active_wgsrpd3_regions": len(active_regions),
            "contemporary_only_regions_excluded_from_pairwise_test": len(contemporary_only_regions),
            "pairwise_region_comparisons": pair_count,
            "added_butterfly_x_region_units": total_added,
        },
        "observed": {
            "mean_pairwise_jaccard_native": native_mean,
            "mean_pairwise_jaccard_contemporary": contemporary_mean,
            "delta_mean_jaccard": contemporary_mean - native_mean,
            "relative_change_mean_jaccard": (contemporary_mean / native_mean - 1.0) if native_mean else None,
            "median_pairwise_jaccard_native": native_median,
            "median_pairwise_jaccard_contemporary": contemporary_median,
            "mean_butterfly_pair_resource_geography_jaccard_native": native_species_mean,
            "mean_butterfly_pair_resource_geography_jaccard_contemporary": contemporary_species_mean,
            "delta_mean_butterfly_resource_niche_overlap": contemporary_species_mean - native_species_mean,
            "relative_change_mean_butterfly_resource_niche_overlap": (
                contemporary_species_mean / native_species_mean - 1.0
                if native_species_mean else None
            ),
            "median_butterfly_pair_resource_geography_jaccard_native": native_species_median,
            "median_butterfly_pair_resource_geography_jaccard_contemporary": contemporary_species_median,
            "butterfly_pair_comparisons": species_pair_count,
        },
        "fixed_margin_null": {
            **null,
            "preserves": [
                "each butterfly's number of added WGSRPD3 regions",
                "each WGSRPD3 region's number of added butterfly resource opportunities",
                "native resource incidences as structural exclusions",
            ],
            "observed_minus_null_median": contemporary_mean - float(null["mean_similarity_null_median"]),
            "p_greater_or_equal": p_greater,
            "observed_butterfly_niche_overlap_minus_null_median": (
                contemporary_species_mean
                - float(null["mean_butterfly_niche_overlap_null_median"])
            ),
            "p_butterfly_niche_overlap_greater_or_equal": p_species_greater,
        },
        "top38_half_credit_hosts": {
            "hosts": len(top38_rows),
            "fractional_credit_share": top38_credit / total_credit,
            "hosts_contributing_to_2plus_butterfly_families": sum(
                int(r["contributing_butterfly_families"]) >= 2 for r in top38_rows
            ),
            "hosts_contributing_to_3plus_butterfly_families": sum(
                int(r["contributing_butterfly_families"]) >= 3 for r in top38_rows
            ),
            "hosts_contributing_to_4plus_butterfly_families": sum(
                int(r["contributing_butterfly_families"]) >= 4 for r in top38_rows
            ),
            "fraction_top38_credit_from_2plus_family_hosts": credit_multifamily / top38_credit,
            "mean_pairwise_jaccard_using_only_top38_introduced_additions": top_mean,
            "delta_mean_jaccard_using_only_top38": top_mean - native_mean,
            "mean_pairwise_jaccard_using_only_remaining_introduced_additions": rem_mean,
            "delta_mean_jaccard_using_only_remaining_hosts": rem_mean - native_mean,
        },
        "top_hosts": host_rows[:20],
        "claim_boundary": (
            "This is a post-hoc analysis of reconstructed larval-resource opportunity, not "
            "realized butterfly community composition. Butterfly-pair resource-geography "
            "Jaccard is a potential trophic-geographic niche-overlap metric, not evidence of "
            "realized competition. Higher regional similarity indicates "
            "homogenization of potential resource opportunity. The fixed-margin null controls "
            "the amount of expansion per butterfly and per region but does not establish a "
            "historical causal effect on butterfly colonization or fitness."
        ),
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with args.output_host_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(host_rows[0].keys()))
        writer.writeheader()
        writer.writerows(host_rows)

    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

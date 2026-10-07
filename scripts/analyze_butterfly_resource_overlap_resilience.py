#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median


def load_descriptors(path: Path) -> dict[str, dict[str, float | int]]:
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sp = str(row["species"]).strip()
            out[sp] = {
                "host_family_count": float(row["host_family_count"]),
                "native_resource_units": int(row["host_wgsrpd3_unit_count"]),
            }
    return out


def load_taxonomy(path: Path) -> dict[str, str]:
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sp = str(row["Species"]).strip()
            fam = str(row["Family"]).strip()
            if sp and fam:
                out[sp] = fam
    return out


def load_pairs(path: Path):
    by_butterfly: dict[str, set[str]] = defaultdict(set)
    meta: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sp = str(row["insect_species"]).strip()
            host = str(row["accepted_plant_name_id"]).strip()
            if not sp or not host:
                continue
            by_butterfly[sp].add(host)
            meta.setdefault(host, {
                "accepted_name": str(row["accepted_name"]).strip(),
                "family": str(row["family"]).strip(),
            })
    return by_butterfly, meta


def load_units(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            host = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if host and unit:
                out[host].add(unit)
    return out


def quantile(xs: list[float], p: float) -> float:
    ys = sorted(xs)
    if not ys:
        raise RuntimeError("empty quantile")
    pos = p * (len(ys) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(ys) - 1)
    frac = pos - lo
    return ys[lo] * (1 - frac) + ys[hi] * frac


def summarize_counts(xs: list[int]) -> dict[str, float | int]:
    ys = sorted(xs)
    return {
        "n": len(ys),
        "mean": sum(ys) / len(ys) if ys else 0.0,
        "median": median(ys) if ys else 0.0,
        "q25": quantile([float(x) for x in ys], 0.25) if ys else 0.0,
        "q75": quantile([float(x) for x in ys], 0.75) if ys else 0.0,
        "zero_fraction": (sum(x == 0 for x in ys) / len(ys)) if ys else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--randomizations", type=int, default=499)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-host-csv", type=Path, required=True)
    args = ap.parse_args()

    descriptors = load_descriptors(args.descriptors_csv)
    taxonomy = load_taxonomy(args.leptraits_csv)
    pairs, plant_meta = load_pairs(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)

    species = [
        sp for sp, d in descriptors.items()
        if float(d["host_family_count"]) > 0 and int(d["native_resource_units"]) > 0
    ]
    if len(species) != 239:
        raise RuntimeError(f"expected 239 butterflies; got {len(species)}")
    focal = set(species)
    species_index = {sp: i for i, sp in enumerate(species)}

    host_consumers: dict[str, list[str]] = defaultdict(list)
    for sp in species:
        if sp not in taxonomy:
            raise RuntimeError(f"missing butterfly family for {sp}")
        for host in pairs.get(sp, set()):
            host_consumers[host].append(sp)

    native_union: dict[str, set[str]] = {}
    contemporary_union: dict[str, set[str]] = {}
    added_union: dict[str, set[str]] = {}
    unit_contributors: dict[str, dict[str, set[str]]] = {}
    host_credit: dict[str, float] = defaultdict(float)
    host_consumer_species: dict[str, set[str]] = defaultdict(set)
    host_consumer_families: dict[str, set[str]] = defaultdict(set)

    for sp in species:
        hosts = sorted(pairs.get(sp, set()))
        n_union: set[str] = set()
        c_union: set[str] = set()
        added_by_host: dict[str, set[str]] = {}
        for host in hosts:
            n = set(native.get(host, set()))
            c = set(contemporary.get(host, set()))
            n_union |= n
            c_union |= c
            added_by_host[host] = c - n
        if len(n_union) != int(descriptors[sp]["native_resource_units"]):
            raise RuntimeError(f"native union drift for {sp}")
        added = c_union - n_union
        native_union[sp] = n_union
        contemporary_union[sp] = c_union
        added_union[sp] = added

        contrib_by_unit: dict[str, set[str]] = {}
        for unit in added:
            contributors = {h for h in hosts if unit in added_by_host[h]}
            if not contributors:
                raise RuntimeError(f"no contributor for {sp}/{unit}")
            contrib_by_unit[unit] = contributors
            share = 1.0 / len(contributors)
            for host in contributors:
                host_credit[host] += share
                host_consumer_species[host].add(sp)
                host_consumer_families[host].add(taxonomy[sp])
        unit_contributors[sp] = contrib_by_unit

    total_added = sum(len(v) for v in added_union.values())
    if total_added != 14553:
        raise RuntimeError(f"expected 14553 added units; got {total_added}")

    # Potential interspecific resource overlap: butterfly pairs sharing at least
    # one known host in the same region. This is exposure/opportunity, not evidence
    # of realized competition.
    native_pair_regions: dict[tuple[int, int], set[str]] = defaultdict(set)
    contemporary_pair_regions: dict[tuple[int, int], set[str]] = defaultdict(set)
    for host, consumers0 in host_consumers.items():
        consumers = sorted(set(consumers0) & focal)
        if len(consumers) < 2:
            continue
        pairs_idx = [
            tuple(sorted((species_index[a], species_index[b])))
            for a, b in combinations(consumers, 2)
        ]
        for pair in pairs_idx:
            native_pair_regions[pair].update(native.get(host, set()))
            contemporary_pair_regions[pair].update(contemporary.get(host, set()))

    native_pair_region_units = sum(len(v) for v in native_pair_regions.values())
    contemporary_pair_region_units = sum(len(v) for v in contemporary_pair_regions.values())
    novel_pair_regions: dict[tuple[int, int], set[str]] = {}
    for pair, regs in contemporary_pair_regions.items():
        novel = regs - native_pair_regions.get(pair, set())
        if novel:
            novel_pair_regions[pair] = novel
    novel_pair_region_units = sum(len(v) for v in novel_pair_regions.values())

    within_family_novel = 0
    between_family_novel = 0
    for (i, j), regs in novel_pair_regions.items():
        if taxonomy[species[i]] == taxonomy[species[j]]:
            within_family_novel += len(regs)
        else:
            between_family_novel += len(regs)

    # How crowded are native vs newly added resource opportunities, considering
    # only focal butterflies that share at least one exact host species?
    native_comp: dict[tuple[str, str], set[str]] = defaultdict(set)
    added_comp: dict[tuple[str, str], set[str]] = defaultdict(set)
    for host, consumers0 in host_consumers.items():
        consumers = sorted(set(consumers0) & focal)
        if len(consumers) < 2:
            continue
        for unit in native.get(host, set()):
            for sp in consumers:
                if unit in native_union[sp]:
                    native_comp[(sp, unit)].update(x for x in consumers if x != sp)
        for unit in contemporary.get(host, set()):
            for sp in consumers:
                if unit in added_union[sp]:
                    added_comp[(sp, unit)].update(x for x in consumers if x != sp)

    native_crowding = [
        len(native_comp.get((sp, unit), set()))
        for sp in species for unit in native_union[sp]
    ]
    added_crowding = [
        len(added_comp.get((sp, unit), set()))
        for sp in species for unit in added_union[sp]
    ]

    # Rank contributing hosts by fractional credit.
    ranked_hosts = sorted(host_credit, key=lambda h: (-host_credit[h], h))
    if len(ranked_hosts) != 670:
        raise RuntimeError(f"expected 670 contributing hosts; got {len(ranked_hosts)}")

    total_credit = sum(host_credit.values())
    cum = 0.0
    half_k = None
    for k, host in enumerate(ranked_hosts, start=1):
        cum += host_credit[host]
        if cum >= 0.5 * total_credit:
            half_k = k
            break
    if half_k != 38:
        raise RuntimeError(f"expected half-credit k=38; got {half_k}")

    host_bit = {h: 1 << i for i, h in enumerate(ranked_hosts)}
    unit_masks: list[tuple[int, int]] = []
    baseline_added_by_species = [0 for _ in species]
    for sp in species:
        si = species_index[sp]
        baseline_added_by_species[si] = len(added_union[sp])
        for contributors in unit_contributors[sp].values():
            mask = 0
            for h in contributors:
                mask |= host_bit[h]
            unit_masks.append((si, mask))

    def removal_effect(removed: set[str]) -> dict[str, object]:
        rmask = 0
        for h in removed:
            rmask |= host_bit[h]
        lost_by_species = [0 for _ in species]
        for si, mask in unit_masks:
            if mask & ~rmask == 0:
                lost_by_species[si] += 1
        lost = sum(lost_by_species)
        affected = sum(x > 0 for x in lost_by_species)
        frac25 = frac50 = frac100 = 0
        for si, lost_i in enumerate(lost_by_species):
            base = baseline_added_by_species[si]
            if base <= 0:
                continue
            frac = lost_i / base
            frac25 += frac >= 0.25
            frac50 += frac >= 0.50
            frac100 += frac >= 1.0
        return {
            "lost_added_butterfly_x_region_units": lost,
            "fraction_of_all_added_units_lost": lost / total_added,
            "butterflies_losing_any_added_opportunity": affected,
            "butterflies_losing_at_least_25pct_added_opportunity": frac25,
            "butterflies_losing_at_least_50pct_added_opportunity": frac50,
            "butterflies_losing_all_added_opportunity": frac100,
        }

    rng = random.Random(args.seed)
    ks = [1, 5, 10, 20, 38, 50, 100]
    removal = {}
    for k in ks:
        top_set = set(ranked_hosts[:k])
        observed = removal_effect(top_set)
        null_lost = []
        null_affected = []
        null_50 = []
        for _ in range(args.randomizations):
            chosen = set(rng.sample(ranked_hosts, k))
            x = removal_effect(chosen)
            null_lost.append(int(x["lost_added_butterfly_x_region_units"]))
            null_affected.append(int(x["butterflies_losing_any_added_opportunity"]))
            null_50.append(int(x["butterflies_losing_at_least_50pct_added_opportunity"]))
        observed["random_host_removal_null"] = {
            "randomizations": args.randomizations,
            "lost_units_median": median(null_lost),
            "lost_units_ci95": [quantile([float(x) for x in null_lost], 0.025), quantile([float(x) for x in null_lost], 0.975)],
            "p_random_lost_ge_observed": (1 + sum(x >= int(observed["lost_added_butterfly_x_region_units"]) for x in null_lost)) / (args.randomizations + 1),
            "affected_species_median": median(null_affected),
            "species_losing_50pct_median": median(null_50),
        }
        removal[str(k)] = observed

    host_rows = []
    for rank, host in enumerate(ranked_hosts, start=1):
        host_rows.append({
            "rank": rank,
            "accepted_plant_name_id": host,
            "accepted_name": plant_meta.get(host, {}).get("accepted_name", ""),
            "plant_family": plant_meta.get(host, {}).get("family", ""),
            "fractional_credit": host_credit[host],
            "share_of_all_added_units": host_credit[host] / total_credit,
            "focal_butterfly_species": len(host_consumer_species[host]),
            "focal_butterfly_families": len(host_consumer_families[host]),
            "butterfly_families": ";".join(sorted(host_consumer_families[host])),
        })

    payload = {
        "schema": "chocho_butterfly_resource_overlap_resilience_v0.1",
        "status": "SUCCESS_POSTHOC_RESOURCE_OVERLAP_AND_REMOVAL_SCENARIOS",
        "panel": {
            "butterflies": len(species),
            "contributing_hosts": len(ranked_hosts),
            "added_butterfly_x_region_units": total_added,
        },
        "potential_interspecific_resource_overlap": {
            "definition": (
                "A butterfly-pair×WGSRPD3 unit is counted when both focal butterflies are "
                "documented to use at least one exact same host species that is present in "
                "that region. This is potential resource-sharing exposure, not observed competition."
            ),
            "native_pair_x_region_units": native_pair_region_units,
            "contemporary_pair_x_region_units": contemporary_pair_region_units,
            "novel_pair_x_region_units_created_by_introduced_host_ranges": novel_pair_region_units,
            "relative_increase_pair_x_region_units": (
                contemporary_pair_region_units / native_pair_region_units - 1.0
                if native_pair_region_units else None
            ),
            "butterfly_pairs_with_native_shared_resource_geography": len(native_pair_regions),
            "butterfly_pairs_with_contemporary_shared_resource_geography": len(contemporary_pair_regions),
            "butterfly_pairs_gaining_at_least_one_novel_shared_resource_region": len(novel_pair_regions),
            "novel_pair_x_region_units_within_butterfly_family": within_family_novel,
            "novel_pair_x_region_units_between_butterfly_families": between_family_novel,
        },
        "resource_space_crowding": {
            "native_resource_units": summarize_counts(native_crowding),
            "introduced_added_resource_units": summarize_counts(added_crowding),
            "interpretation": (
                "Counts are numbers of other focal butterflies sharing at least one exact "
                "host species in the same regional resource opportunity. Zero does not mean "
                "absence of all ecological competitors outside the focal panel."
            ),
        },
        "introduced_host_removal_scenarios": {
            "definition": (
                "Remove only introduced-range contributions of ranked host plants while "
                "retaining their native ranges and all other known hosts. Loss occurs only "
                "where no unremoved contributing host can support the same butterfly×region unit."
            ),
            "top_host_ranking": "unit-preserving fractional contribution to the 14,553 added butterfly×region units",
            "scenarios": removal,
            "management_boundary": (
                "These coarse regional scenarios identify potential resource gaps and network "
                "redundancy; they do not justify retaining invasive plants or predict local "
                "population responses. Site-specific management requires local host use, host "
                "quality, abundance and demographic data."
            ),
        },
        "top_hosts": host_rows[:50],
        "claim_boundary": (
            "All competition language is potential resource overlap only. The data contain "
            "no butterfly abundance or density response and cannot estimate realized intra- "
            "or interspecific competition. Removal scenarios quantify resource-network "
            "sensitivity, not demographic effects."
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

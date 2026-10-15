#!/usr/bin/env python3
"""Plant-range fixed-degree null for butterfly shared-host and redundancy geography.

Post-hoc *structural* sensitivity, deliberately separate from the submitted
GEB paper. No larval success/competition inference is licensed here.

Swaps among introduced plant × region cells are constrained within plant
family × WGSRPD Level1. This preserves (a) recorded butterfly–host links and
their network degrees; (b) each plant's introduced-cell count *per continent*;
(c) each region's introduced-host count per plant family; (d) all native host
distributions as exclusions. The null does NOT preserve each butterfly's
added-region count or ecological suitability.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import random
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_host_geography(
    rows: list[dict[str, str]], regions: dict[str, int]
) -> dict[str, set[int]]:
    result: dict[str, set[int]] = defaultdict(set)
    for r in rows:
        h = r["accepted_plant_name_id"].strip()
        u = r["area_code_l3"].strip()
        if h and u:
            result[h].add(regions[u])
    return result


def percentiles(values: list[float]) -> tuple[float, float, float]:
    x = sorted(values)
    if not x:
        raise ValueError("empty null")
    def q(p: float) -> float:
        t = p * (len(x) - 1)
        lo = math.floor(t)
        hi = math.ceil(t)
        return x[lo] * (hi - t) + x[hi] * (t - lo) if lo != hi else x[lo]
    return q(0.025), float(median(x)), q(0.975)


def site_contribution_stats(
    species_hosts: list[list[str]],
    native_by_host: dict[str, int],
    introduced_by_host: dict[str, int],
    shared_pairs: list[tuple[int, int, list[str], int]],
    native_butterfly: list[int],
) -> dict[str, float | int]:
    """Bitmask union counting; unique host×region pairs count only once."""
    added_count = 0
    single_count = 0
    added_masks: list[int] = []
    for sp, hostlist in enumerate(species_hosts):
        seen = 0
        duplicate = 0
        baseline = native_butterfly[sp]
        for h in hostlist:
            here = introduced_by_host.get(h, 0) & ~baseline
            duplicate |= seen & here
            seen |= here
        added_masks.append(seen)
        added_count += seen.bit_count()
        single_count += (seen & ~duplicate).bit_count()

    # Same-denominator measure: for every newly gained butterfly × region
    # cell, count DISTINCT other focal butterflies able to share an exact
    # host there. A (butterfly_i, butterfly_j, region) pair contributes at
    # most one exposure even if several plant species are shared.
    shared_new = 0
    added_peer_exposures = 0
    added_with_peers = [0] * len(species_hosts)
    for i, j, common_hosts, baseline in shared_pairs:
        now = 0
        for h in common_hosts:
            now |= native_by_host.get(h, 0) | introduced_by_host.get(h, 0)
        shared_new += (now & ~baseline).bit_count()
        exposure_i = now & added_masks[i]
        exposure_j = now & added_masks[j]
        added_peer_exposures += exposure_i.bit_count() + exposure_j.bit_count()
        added_with_peers[i] |= exposure_i
        added_with_peers[j] |= exposure_j
    any_peer = sum(x.bit_count() for x in added_with_peers)

    return {
        "added_units": added_count,
        "single_units": single_count,
        "single_fraction": single_count / added_count if added_count else float("nan"),
        "novel_exact_shared_pair_regions": shared_new,
        "shared_pair_units_per_added_unit": (
            shared_new / added_count if added_count else float("nan")
        ),
        "added_cell_total_other_shared_butterflies": added_peer_exposures,
        "added_cell_mean_other_shared_butterflies": (
            added_peer_exposures / added_count if added_count else float("nan")
        ),
        "added_cell_fraction_with_any_shared_butterfly": (
            any_peer / added_count if added_count else float("nan")
        ),
        "added_cell_fraction_without_shared_butterfly": (
            (added_count - any_peer) / added_count
            if added_count else float("nan")
        ),
    }


def make_groups(
    native_sets: dict[str, set[int]],
    introduced_sets: dict[str, set[int]],
    family: dict[str, str],
    level1: dict[int, str],
) -> tuple[list[list[tuple[str, int]]], list[int], int]:
    groups: dict[tuple[str, str], list[tuple[str, int]]] = defaultdict(list)
    for h, cols in introduced_sets.items():
        if not h in family:
            raise ValueError("missing host-plant family")
        for u in cols:
            if u in native_sets.get(h, set()):
                raise ValueError("introduced/native overlap")
            groups[(family[h], level1[u])].append((h, u))
    viable = [
        list(sorted(g))
        for g in groups.values()
        if len(g) >= 2 and len(set(h for h, _ in g)) >= 2
        and len(set(u for _, u in g)) >= 2
    ]
    viable.sort(key=lambda x: (x[0][0], x[0][1]))
    cumulative = []
    total = 0
    for g in viable:
        total += len(g)
        cumulative.append(total)
    return viable, cumulative, sum(len(g) for g in groups.values())


def attempt_range_swaps(
    groups: list[list[tuple[str, int]]],
    cumulative: list[int],
    native_sets: dict[str, set[int]],
    introduced_sets: dict[str, set[int]],
    attempts: int,
    rng: random.Random,
) -> int:
    """One reversible double-edge swap on a single family×continent stratum."""
    if not cumulative:
        return 0
    accepted = 0
    for _ in range(attempts):
        k = rng.randrange(cumulative[-1])
        gi = bisect.bisect_right(cumulative, k)
        lower = cumulative[gi - 1] if gi else 0
        group = groups[gi]
        i = k - lower
        j = rng.randrange(len(group))
        if i == j:
            continue
        a, x = group[i]
        b, y = group[j]
        if a == b or x == y:
            continue
        if y in native_sets.get(a, set()) or x in native_sets.get(b, set()):
            continue
        if y in introduced_sets[a] or x in introduced_sets[b]:
            continue
        introduced_sets[a].remove(x)
        introduced_sets[b].remove(y)
        introduced_sets[a].add(y)
        introduced_sets[b].add(x)
        group[i] = (a, y)
        group[j] = (b, x)
        accepted += 1
    return accepted


def build_model(
    descriptors: list[dict[str, str]],
    butterfly_host_rows: list[dict[str, str]],
    native_rows: list[dict[str, str]],
    contemporary_rows: list[dict[str, str]],
    level1_map: dict[str, str],
) -> dict[str, object]:
    species = [
        r["species"] for r in descriptors
        if float(r["host_family_count"]) > 0
        and int(r["host_wgsrpd3_unit_count"]) > 0
    ]
    if len(species) != 239:
        raise ValueError(f"unexpected focal species count: {len(species)}")
    sp_index = {s: i for i, s in enumerate(species)}
    species_hosts = [set() for _ in species]
    host_families = {}
    host_consumers: dict[str, set[int]] = defaultdict(set)
    for row in butterfly_host_rows:
        s = row["insect_species"].strip()
        h = row["accepted_plant_name_id"].strip()
        fam = row["family"].strip()
        if not h or s not in sp_index:
            continue
        if h in host_families and host_families[h] != fam:
            raise ValueError("inconsistent plant-family assignment")
        host_families[h] = fam
        i = sp_index[s]
        species_hosts[i].add(h)
        host_consumers[h].add(i)
    species_hosts_ordered = [sorted(h) for h in species_hosts]
    if any(not h for h in species_hosts_ordered):
        raise ValueError("focal butterfly without resolved hosts")
    hosts = set(host_consumers)
    relevant_native = [r for r in native_rows if r["accepted_plant_name_id"] in hosts]
    relevant_cont = [
        r for r in contemporary_rows if r["accepted_plant_name_id"] in hosts
    ]
    region_codes = sorted(
        {r["area_code_l3"] for r in relevant_native + relevant_cont}
    )
    if not set(region_codes) <= set(level1_map):
        raise ValueError("WGSRPD Level1 map misses resource units")
    region_index = {u: j for j, u in enumerate(region_codes)}
    lev = {region_index[u]: str(level1_map[u]) for u in region_codes}
    native_sets = read_host_geography(relevant_native, region_index)
    contemporary_sets = read_host_geography(relevant_cont, region_index)
    introduced_sets = {
        h: contemporary_sets.get(h, set()) - native_sets.get(h, set())
        for h in hosts
    }
    if any(
        not native_sets.get(h, set()) <= contemporary_sets.get(h, set())
        for h in hosts
    ):
        raise ValueError("native subset contemporary invariant violated")

    def mask(vals: set[int]) -> int:
        result = 0
        for u in vals:
            result |= 1 << u
        return result

    native_bits = {h: mask(native_sets.get(h, set())) for h in hosts}
    introduced_bits = {h: mask(introduced_sets[h]) for h in hosts}

    native_species = []
    for i, hostlist in enumerate(species_hosts_ordered):
        nm = 0
        for h in hostlist:
            nm |= native_bits[h]
        if nm.bit_count() != int(
            next(x["host_wgsrpd3_unit_count"] for x in descriptors
                 if x["species"] == species[i])
        ):
            raise ValueError("native species-resource union drift")
        native_species.append(nm)

    shared_hosts: dict[tuple[int, int], list[str]] = defaultdict(list)
    for h, users in host_consumers.items():
        if len(users) >= 2:
            for i, j in combinations(sorted(users), 2):
                shared_hosts[(i, j)].append(h)
    shared_pairs = []
    for (i, j), common in sorted(shared_hosts.items()):
        baseline = 0
        for h in common:
            baseline |= native_bits[h]
        shared_pairs.append((i, j, sorted(common), baseline))

    groups, cumulative, total_edges = make_groups(
        native_sets, introduced_sets, host_families, lev
    )
    return {
        "species": species,
        "species_hosts": species_hosts_ordered,
        "native_bits": native_bits,
        "introduced_sets": introduced_sets,
        "native_sets": native_sets,
        "shared_pairs": shared_pairs,
        "native_species": native_species,
        "groups": groups,
        "cumulative": cumulative,
        "total_introduced_edges": total_edges,
        "host_families": host_families,
        "level1": lev,
    }


def run_null(model: dict[str, object], reps: int, seed: int,
             burn_factor: int, step_factor: int) -> dict[str, object]:
    if reps < 19 or burn_factor < 1 or step_factor < 1:
        raise ValueError("insufficient declared null runs")
    hosts = model["species_hosts"]
    native_bits = model["native_bits"]
    introduced_sets = model["introduced_sets"]
    shared = model["shared_pairs"]
    native_species = model["native_species"]
    groups = model["groups"]
    cumulative = model["cumulative"]
    native_sets = model["native_sets"]
    E = model["total_introduced_edges"]
    rng = random.Random(seed)

    def introduced_bits() -> dict[str, int]:
        return {
            h: sum(1 << u for u in cols)
            for h, cols in introduced_sets.items()
        }
    def evaluate() -> dict[str, float | int]:
        return site_contribution_stats(
            hosts, native_bits, introduced_bits(), shared, native_species
        )
    observed = evaluate()
    # Frozen receipts from the exact primary reconstruction, NOT null outcomes.
    if observed["added_units"] != 14553 or observed["single_units"] != 8574 \
        or observed["novel_exact_shared_pair_regions"] != 79412:
        raise ValueError(f"empirical reconstruction drift: {observed}")
    # Frozen earlier calculation: number of OTHER focal butterflies sharing
    # exact host plants at each newly added butterfly × region opportunity.
    if abs(observed["added_cell_mean_other_shared_butterflies"]
           - 5.712567855424998) > 1e-10:
        raise ValueError("observed new-resource co-user mean drift")
    if abs(observed["added_cell_fraction_without_shared_butterfly"]
           - 0.09894867037724181) > 1e-10:
        raise ValueError("observed new-resource co-user zero fraction drift")
    baseline_edges = {h: frozenset(cols) for h, cols in introduced_sets.items()}
    burn_attempts = burn_factor * E
    accepted_burn = attempt_range_swaps(
        groups, cumulative, native_sets, introduced_sets, burn_attempts, rng
    )
    changed = sum(
        len(baseline_edges[h] ^ introduced_sets[h])
        for h in introduced_sets
    )
    # Each moved edge is removed from one set and added to another, so this
    # fraction represents disturbed incidence cells, not number of pairs.
    moved_fraction = changed / (2 * E) if E else 0.
    if moved_fraction < 0.03 or accepted_burn < max(100, E // 20):
        return {
            "status": "NOT_EVALUABLE_INADEQUATE_MIXING",
            "observed": observed,
            "total_introduced_edges": E,
            "viable_edges": cumulative[-1] if cumulative else 0,
            "burn_accepted": accepted_burn,
            "changed_fraction_after_burn": moved_fraction,
            "rule": "burn must accept >=max(100,E/20) swaps and move >=3% of all introduced incidences",
        }
    null: dict[str, list[float]] = {
        name: [] for name in observed
    }
    accepted_step = []
    for _ in range(reps):
        accepted_step.append(attempt_range_swaps(
            groups, cumulative, native_sets, introduced_sets, step_factor * E, rng
        ))
        ans = evaluate()
        for name in null:
            null[name].append(float(ans[name]))

    results = {}
    for k, values in null.items():
        qlo, center, qhi = percentiles(values)
        obs = float(observed[k])
        results[k] = {
            "observed": obs,
            "null_q025": qlo,
            "null_median": center,
            "null_q975": qhi,
            "observed_minus_null_median": obs - center,
            "p_ge": (1 + sum(x >= obs for x in values)) / (reps + 1),
            "p_le": (1 + sum(x <= obs for x in values)) / (reps + 1),
        }
    return {
        "status": "EXPLORATORY_PLANT_RANGE_FIXED_DEGREE_NULL",
        "observed": observed,
        "total_introduced_edges": E,
        "viable_edges": cumulative[-1] if cumulative else 0,
        "seed": seed,
        "replicates": reps,
        "burn_proposals": burn_attempts,
        "burn_accepted": accepted_burn,
        "changed_fraction_after_burn": moved_fraction,
        "step_proposals_per_rep": step_factor * E,
        "step_accepted_min": min(accepted_step),
        "step_accepted_median": median(accepted_step),
        "null_summary": results,
        "preserves": [
            "butterfly-host interaction identities and degrees",
            "introduced geographic cells per plant per WGSRPD Level1",
            "introduced plant richness per region per plant family",
            "all native plant occurrences as forbidden introduced cells",
        ],
        "does_not_preserve": [
            "individual butterfly added-resource region count",
            "realistic plant climate niche or time of introduction",
            "observed butterfly occupancy, larval success or host phenology",
        ],
        "interpretation": (
            "Post-hoc plant-family×continent range-shuffle diagnostic. "
            "The new-cell co-user density and zero-co-user fraction share "
            "an identical butterfly×region denominator and are preferable "
            "to the cross-unit shared pair×region / added-cell ratio. "
            "All results are structural opportunity, not local co-occurrence "
            "or competitive intensity."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    for key in (
        "descriptors_csv", "insect_host_csv", "native_distribution_csv",
        "contemporary_distribution_csv", "level1_json", "output_json"
    ):
        p.add_argument("--" + key.replace("_", "-"), type=Path, required=True)
    p.add_argument("--replicates", type=int, default=199)
    p.add_argument("--seed", type=int, default=20261015)
    p.add_argument("--burn-factor", type=int, default=5)
    p.add_argument("--step-factor", type=int, default=1)
    args = p.parse_args()
    lookup = json.loads(args.level1_json.read_text(encoding="utf-8"))["units"]
    model = build_model(
        read_csv(args.descriptors_csv), read_csv(args.insect_host_csv),
        read_csv(args.native_distribution_csv),
        read_csv(args.contemporary_distribution_csv), lookup
    )
    out = run_null(
        model, args.replicates, args.seed, args.burn_factor, args.step_factor
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Post-hoc plant-range fixed-margin null for butterfly geographic host redundancy.

Randomizes the *locations* of recorded introductions, NOT butterfly host identity.
Preserves, within each WCVP plant family and WGSRPD level-1 continent:
- each host plant's number of introduced WGSRPD3 records,
- each WGSRPD3 region's number of introduced host plants in that family,
- all existing native plant occurrences as structural exclusions.
Records containing the same plant in a region are prohibited.

Requires frozen WCVP/HOSTS sidecars from build_wcvp_hosts_null_sidecars.R.
Run within the isolated analysis branch; NOT a submitted-manuscript primary test.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


def csv_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_units(path: Path):
    values = defaultdict(set)
    for r in csv_rows(path):
        host = r["accepted_plant_name_id"].strip()
        region = r["area_code_l3"].strip()
        if host and region:
            values[host].add(region)
    return dict(values)


def read_hosts(path: Path):
    by_species = defaultdict(set)
    host_family = {}
    for r in csv_rows(path):
        sp, h, family = (r[k].strip() for k in
                         ("insect_species", "accepted_plant_name_id", "family"))
        if not sp or not h or not family:
            continue
        if h in host_family and host_family[h] != family:
            raise RuntimeError("Conflicting WCVP family ID for host " + h)
        by_species[sp].add(h)
        host_family[h] = family
    return dict(by_species), host_family


def read_descriptors(path: Path):
    return {r["species"]: r for r in csv_rows(path)}


def quantile(xs, p):
    ys = sorted(xs)
    if not ys:
        raise ValueError("Cannot take quantile of empty data")
    pos = p * (len(ys) - 1)
    low = int(pos)
    high = min(low + 1, len(ys) - 1)
    return ys[low] + (ys[high] - ys[low]) * (pos - low)


def init_edges(native, contemporary, host_family, continent):
    """Use the full resolved HOSTS-WCVP pool, not only focal butterfly plants."""
    introduced = {}
    for h in host_family:
        n = native.get(h, set())
        c = contemporary.get(h, set())
        if not n <= c:
            raise RuntimeError("Native host region not contemporary for " + h)
        introduced[h] = set(c - n)

    blocks = defaultdict(list)
    for host, regs in sorted(introduced.items()):
        for region in sorted(regs):
            if region not in continent:
                raise RuntimeError("Missing fixed WGSRPD1 lookup for " + region)
            blocks[(host_family[host], continent[region])].append((host, region))
    return introduced, dict(blocks)


def geographic_resources(species, species_hosts, native, introduced, selected):
    """Recount exact introduced host contributors on regions lacking any native host."""
    metrics = {}
    for sp in species:
        n_union = set()
        for h in species_hosts.get(sp, set()):
            n_union.update(native.get(h, set()))
        counts = Counter()
        for h in species_hosts.get(sp, set()):
            for r in introduced.get(h, set()):
                if r not in n_union:
                    counts[r] += 1
        units = len(counts)
        singles = sum(v == 1 for v in counts.values())
        metrics[sp] = (len(n_union), units, singles)

    summaries = {}
    for name, members in selected.items():
        vals = [metrics[sp] for sp in members]
        expanded = [v for v in vals if v[1] > 0]
        total = sum(v[1] for v in vals)
        single = sum(v[2] for v in vals)
        summaries[name] = {
            "species": len(members),
            "expanded_species": len(expanded),
            "native_units": sum(v[0] for v in vals),
            "added_units": total,
            "single_host_added_units": single,
            "single_host_share": single / total if total else None,
            "mean_species_single_share": (
                sum(v[2] / v[1] for v in expanded) / len(expanded)
                if expanded else None
            ),
        }
    return summaries, metrics


def build_swappable_blocks(blocks):
    out = {}
    for key in sorted(blocks):
        edges = blocks[key]
        if len(edges) > 1 and len({h for h, _ in edges}) > 1 and len({r for _, r in edges}) > 1:
            out[key] = edges
    return out


def swaps(blocks, native, introduced, rng, sweeps):
    """Double-edge switch; preserves host and region degrees per family/continent."""
    attempts = accepted = 0
    for edges in blocks.values():
        n = len(edges)
        for _ in range(n * sweeps):
            attempts += 1
            i = rng.randrange(n)
            j = rng.randrange(n)
            if i == j:
                continue
            h1, r1 = edges[i]
            h2, r2 = edges[j]
            if h1 == h2 or r1 == r2:
                continue
            if (r2 in introduced[h1] or r1 in introduced[h2] or
                    r2 in native.get(h1, ()) or r1 in native.get(h2, ())):
                continue
            introduced[h1].remove(r1)
            introduced[h2].remove(r2)
            introduced[h1].add(r2)
            introduced[h2].add(r1)
            edges[i] = (h1, r2)
            edges[j] = (h2, r1)
            accepted += 1
    return attempts, accepted


def degree_signature(blocks):
    return {
        key: (
            tuple(sorted(Counter(h for h, _ in edges).items())),
            tuple(sorted(Counter(r for _, r in edges).items())),
        )
        for key, edges in blocks.items()
    }


def summarize_metric(observed, replicates, key):
    actual = observed[key]
    vals = [v[key] for v in replicates]
    if actual is None or any(v is None for v in vals):
        raise RuntimeError("Missing null metric: " + key)
    return {
        "observed": actual,
        "null_median": median(vals),
        "null_q025": quantile(vals, 0.025),
        "null_q975": quantile(vals, 0.975),
        "difference_from_null_median": actual - median(vals),
        "p_observed_ge_null": (1 + sum(actual <= v for v in vals)) / (len(vals) + 1),
        "p_observed_le_null": (1 + sum(actual >= v for v in vals)) / (len(vals) + 1),
    }


def analyze(
    species_hosts, family, native, contemporary, continent, frozen,
    randomizations=199, seed=20261011, burn_sweeps=15, gap_sweeps=3
):
    panel = sorted(sp for sp, row in frozen.items()
                   if float(row["host_family_count"]) > 0 and
                   int(row["native_resource_units"]) > 0)
    if len(panel) != 239:
        raise RuntimeError("Frozen focal panel drift (expected 239)")
    adequate = sorted(sp for sp in panel
                      if int(frozen[sp]["introduced_added_units"]) > 0 and
                      frozen[sp]["host_taxonomy_lower_bound_adequate"] == "True")
    rich = sorted(sp for sp in adequate if int(frozen[sp]["resolved_host_species"]) >= 10)
    if len(adequate) != 191 or len(rich) != 80:
        raise RuntimeError("Frozen adequate/rich panel drift")
    for sp in panel:
        if sp not in species_hosts:
            raise RuntimeError("Focal species missing from fixed HOSTS sidecar: " + sp)

    selected = {"all_239": panel, "adequate_191": adequate, "rich_hosts_80": rich}
    introduced, blocks = init_edges(native, contemporary, family, continent)
    observed, per_sp = geographic_resources(panel, species_hosts, native, introduced, selected)
    expected = {
        "all_239": (26530, 14553, 8574, 206),
        "adequate_191": (None, 13529, 7733, 191),
        "rich_hosts_80": (None, 6224, 2823, 80),
    }
    for tag, (n_native, n_added, n_single, n_expanded) in expected.items():
        row = observed[tag]
        if (n_native is not None and row["native_units"] != n_native or
            row["added_units"] != n_added or
            row["single_host_added_units"] != n_single or
            row["expanded_species"] != n_expanded):
            raise RuntimeError("Frozen native/introduced reconstruction drift for " + tag +
                               ": " + json.dumps(row))

    # Check all individual native envelopes too.
    for sp in panel:
        if per_sp[sp][0] != int(frozen[sp]["native_resource_units"]):
            raise RuntimeError("Native host geography drift for " + sp)

    rng = random.Random(seed)
    movable = build_swappable_blocks(blocks)
    initial_signature = degree_signature(blocks)
    original_edges = {edge for v in blocks.values() for edge in v}
    total_edges = len(original_edges)
    movable_edges = sum(len(v) for v in movable.values())

    attempted, accepted = swaps(movable, native, introduced, rng, burn_sweeps)
    vals = defaultdict(list)
    interval_accepted = []
    altered = []
    for i in range(randomizations):
        at, ac = swaps(movable, native, introduced, rng, gap_sweeps)
        attempted += at
        accepted += ac
        interval_accepted.append(ac)
        summary, _ = geographic_resources(panel, species_hosts, native, introduced, selected)
        for tag, m in summary.items():
            vals[tag].append(m)
        if i in (0, randomizations // 2, randomizations - 1):
            altered.append({
                "replicate": i + 1,
                "number_original_introduction_edges_changed":
                    total_edges - sum(edge in original_edges
                                      for edges in blocks.values() for edge in edges),
            })

    if degree_signature(blocks) != initial_signature:
        raise RuntimeError("Fixed margins were not preserved")
    for h in introduced:
        if introduced[h] & native.get(h, set()):
            raise RuntimeError("Introduced edge collides with native range: " + h)

    metrics = ("added_units", "single_host_share", "mean_species_single_share")
    comparisons = {
        tag: {key: summarize_metric(observed[tag], vals[tag], key) for key in metrics}
        for tag in selected
    }

    return {
        "schema": "chocho_butterfly_geographic_introduction_redundancy_null_v0.1",
        "status": "POSTHOC_SPATIAL_CONFIGURATION_FIXED_MARGIN_NULL",
        "randomizations": randomizations,
        "seed": seed,
        "burn_sweeps_per_edge": burn_sweeps,
        "between_replica_sweeps_per_edge": gap_sweeps,
        "null": (
            "Swap introduced host-plant x WGSRPD3 incidences ONLY within the same "
            "WCVP plant family and WGSRPD level1 region. Conserve the number of "
            "introduced WGSRPD3 units per plant within level1 and number of "
            "introduced same-family plant taxa per WGSRPD3 region. Native ranges "
            "remain fixed and cannot receive introduced edges. Known butterfly-host "
            "identities, host species counts, and recorded host prominence are fixed."
        ),
        "baseline": observed,
        "comparison": comparisons,
        "mixing": {
            "introduced_edges_in_fixed_HOSTS_WCVP_pool": total_edges,
            "edges_in_eligible_blocks": movable_edges,
            "eligible_block_count": len(movable),
            "total_family_continent_blocks": len(blocks),
            "attempted_swaps": attempted,
            "accepted_swaps": accepted,
            "accepted_fraction": accepted / attempted if attempted else None,
            "median_accepted_swaps_between_replicates": median(interval_accepted),
            "example_changed_original_edge_counts": altered,
            "invariants_confirmed": True,
        },
        "claim_boundary": (
            "These are reconstructed plant ranges at coarse WGSRPD3 scale, not "
            "local butterfly host use, survival, or persistence. The null preserves "
            "marginal introduction counts within plant family and continent, not "
            "fine-scale climate matching, geographic contiguity, or sampling "
            "intensity. Test is exploratory post-hoc and excludes nulls that "
            "reassign butterfly host identities."
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--expansion-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--level1-lookup-json", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=199)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--burn-sweeps", type=int, default=15)
    ap.add_argument("--gap-sweeps", type=int, default=3)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()
    if args.permutations <= 0 or args.burn_sweeps <= 0 or args.gap_sweeps <= 0:
        raise ValueError("Permutation and sweep counts must be positive")
    species_hosts, fam = read_hosts(args.insect_host_csv)
    native = read_units(args.native_distribution_csv)
    contemporary = read_units(args.contemporary_distribution_csv)
    descriptions = read_descriptors(args.descriptors_csv)
    exp = read_descriptors(args.expansion_csv)
    if set(exp) - set(descriptions):
        raise RuntimeError("Frozen expansion species outside S1 descriptor panel")
    frozen = {}
    for sp, e in exp.items():
        frozen[sp] = dict(e)
        frozen[sp]["native_resource_units"] = descriptions[sp]["host_wgsrpd3_unit_count"]
    lookup = json.loads(args.level1_lookup_json.read_text(encoding="utf-8"))["units"]
    report = analyze(
        species_hosts, fam, native, contemporary, lookup, frozen,
        randomizations=args.permutations, seed=args.seed,
        burn_sweeps=args.burn_sweeps, gap_sweeps=args.gap_sweeps,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

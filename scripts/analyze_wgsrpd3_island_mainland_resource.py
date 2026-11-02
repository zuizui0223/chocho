#!/usr/bin/env python3
"""Q1: Butterfly host-resource opportunity on strict island vs mainland WGSRPD3 regions.

Inputs are pinned HOSTS + WCVP v13 sidecars and the frozen 239-species GEB panel.
This measures potential resource presence, never butterfly arrival or occupancy.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

# Positive list from the official pinned WGSRPD3 second-edition Level-3 area names.
# Codes are geographic islands/archipelagos, including large islands (e.g., Japan),
# and are not inferred from a region's butterfly records or outcome.
STRICT_ISLAND_CODES = frozenset("""
ALD ALU AND ARU ASC ASP ATP AZO BAH BAL BER BIS BOR BOU CAY CGS CHH
CKI CNY COM COO COR CPI CRL CRZ CTM CUB CVI CYP DOM DSV EAI EAS FAL
FIJ FOR GAL GGI GIL GNL GRB HAI HAW HBI HMD ICE IRE JAM JAP JAW
JNF KEG KER KRI KUR KZN LDV LEE LIN LSI MAQ MAU MCI MCS MDG
MDR MDV MOL MPE MRN MRQ MRS MXI NCB NFK NFL NLA NNS NRU NUE
NWC NWG NZN NZS OGA PEI PHI PHX PIT PUE REU ROD SAK SAM SAR
SCI SCS SCZ SEL SEY SGE SIC SOC SOL SRL SSA STH SUL SUM SVA
TAI TAS TCI TDC TOK TON TRT TUA TUB TUV VAN VNA WAK WAL WIN XMS
""".split())

# Entire Level-3 regions contain both mainland and islands (or mixed Gulf states).
# Exclude instead of making a post-hoc choice based on measured resource gain.
MIXED_CODES = frozenset("DEN EQG GRC GST NUN TAN".split())


def rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def host_distribution(path: Path) -> dict[str, set[str]]:
    out = defaultdict(set)
    for row in rows(path):
        host = row["accepted_plant_name_id"].strip()
        code = row["area_code_l3"].strip()
        if host and code:
            out[host].add(code)
    return out


def host_by_butterfly(path: Path) -> dict[str, set[str]]:
    out = defaultdict(set)
    for row in rows(path):
        butterfly = row["insect_species"].strip()
        host = row["accepted_plant_name_id"].strip()
        if butterfly and host:
            out[butterfly].add(host)
    return out


def wgsrpd_names(path: Path) -> dict[str, str]:
    out = {}
    with path.open(encoding="latin-1") as handle:
        for line in handle:
            parts = line.strip().split("*")
            if len(parts) > 2 and len(parts[0]) == 3 and parts[0] != "L3 ":
                out[parts[0]] = parts[1]
    return out


def avg(data):
    return None if not data else float(statistics.mean(data))


def med(data):
    return None if not data else float(statistics.median(data))


def metrics_for(regions):
    active = [r for r in regions if r["contemporary"] > 0]
    baselined = [r for r in regions if r["native"] > 0]
    n = sum(r["native"] for r in regions)
    c = sum(r["contemporary"] for r in regions)
    return {
        "all_regions": len(regions),
        "regions_with_contemporary_opportunity": len(active),
        "regions_with_native_baseline": len(baselined),
        "regions_with_new_opportunity": sum(r["added"] > 0 for r in regions),
        "regions_new_from_zero_baseline": sum(r["native"] == 0 and r["contemporary"] > 0 for r in regions),
        "native_butterfly_region_units": n,
        "contemporary_butterfly_region_units": c,
        "added_butterfly_region_units": c - n,
        "pooled_increase_pct_if_native_positive": None if n == 0 else 100 * (c - n) / n,
        "pooled_added_share_pct_if_contemporary_positive": None if c == 0 else 100 * (c - n) / c,
        "mean_added_share_per_active_region": avg([r["added_share"] for r in active]),
        "median_added_share_per_active_region": med([r["added_share"] for r in active]),
        "median_increase_pct_nonzero_native_regions": med(
            [100 * r["added"] / r["native"] for r in baselined]
        ),
    }


def adjusted_mean_diff(regions, labels=None):
    groups = defaultdict(lambda: {"island": [], "mainland": []})
    for i, row in enumerate(regions):
        label = row["class"] if labels is None else labels[i]
        if label in {"island", "mainland"}:
            groups[row["level1"]][label].append(row["added_share"])
    pairs = []
    for name, d in sorted(groups.items()):
        if d["island"] and d["mainland"]:
            weight = min(len(d["island"]), len(d["mainland"]))
            diff = avg(d["island"]) - avg(d["mainland"])
            pairs.append({"level1": name, "island_n": len(d["island"]),
                          "mainland_n": len(d["mainland"]), "mean_difference": diff,
                          "weight": weight})
    denom = sum(x["weight"] for x in pairs)
    value = (sum(x["mean_difference"] * x["weight"] for x in pairs) / denom
             if denom else None)
    return value, pairs


def main():
    p = argparse.ArgumentParser()
    for flag in ("panel-csv", "host-csv", "native-csv", "contemporary-csv",
                 "level1-json", "wgsrpd3-table", "out-dir"):
        p.add_argument("--" + flag, type=Path, required=True)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    butterfly_hosts = host_by_butterfly(a.host_csv)
    native_hosts = host_distribution(a.native_csv)
    contemporary_hosts = host_distribution(a.contemporary_csv)
    panel = list(rows(a.panel_csv))
    level1 = json.loads(a.level1_json.read_text(encoding="utf-8"))["units"]
    names = wgsrpd_names(a.wgsrpd3_table)

    missing_islands = sorted(STRICT_ISLAND_CODES - names.keys())
    # Guard against silently classifying explicitly named island regions
    # as mainland (the earlier v0.1 inadvertently did this for CNY Canary Is.).
    explicit_named_islands = {
        code for code, name in names.items()
        if " Is." in name or " Islands" in name or "Archipelago" in name
    }
    misclassified_named_islands = sorted(explicit_named_islands - STRICT_ISLAND_CODES)
    if misclassified_named_islands:
        raise RuntimeError(f"Named island regions classed mainland: {misclassified_named_islands}")
    if len(STRICT_ISLAND_CODES) != 126 or "CNY" not in STRICT_ISLAND_CODES or "RHO" in STRICT_ISLAND_CODES:
        raise RuntimeError("Incorrect 126-region geographical island scope")
    if missing_islands or len(panel) != 239:
        raise RuntimeError(f"Region codes not in official WGSRPD3={missing_islands}, panel_n={len(panel)}")
    known_codes = set(names) & set(level1)
    native_by_region, contemporary_by_region = Counter(), Counter()
    added_pairs = []
    per_species = []
    mismatches = []
    for row in panel:
        sp = row["species"]
        hosts = butterfly_hosts.get(sp, set())
        n = set().union(*(native_hosts.get(h, set()) for h in hosts))
        c = set().union(*(contemporary_hosts.get(h, set()) for h in hosts))
        expected_n, expected_c = int(row["native_resource_units"]), int(row["contemporary_resource_units"])
        if not n <= c or len(n) != expected_n or len(c) != expected_c:
            mismatches.append({"species": sp, "expected_native": expected_n,
                               "reconstructed_native": len(n),
                               "expected_contemporary": expected_c,
                               "reconstructed_contemporary": len(c),
                               "native_subset": n <= c})
        for code in n:
            native_by_region[code] += 1
        for code in c:
            contemporary_by_region[code] += 1
        for code in sorted(c - n):
            added_pairs.append({"species": sp, "area_code_l3": code,
                                "level1": level1.get(code, "UNKNOWN"),
                                "geographic_group": "island" if code in STRICT_ISLAND_CODES else
                                    ("mixed_excluded" if code in MIXED_CODES else "mainland")})
        per_species.append({"species": sp, "native": len(n), "contemporary": len(c),
                            "introduced_additions": len(c - n)})
    available = set(native_by_region) | set(contemporary_by_region)
    unknown_codes = sorted(available - known_codes)
    diagnostic = {
        "panel_species": len(panel),
        "native_total": sum(native_by_region.values()),
        "contemporary_total": sum(contemporary_by_region.values()),
        "new_total": len(added_pairs),
        "panel_mismatches": mismatches[:25],
        "n_panel_mismatches": len(mismatches),
        "unmapped_resource_codes": unknown_codes,
        "expected_totals": {"native": 26530, "contemporary": 41083, "added": 14553},
    }
    (a.out_dir / "reconstruction_audit.json").write_text(
        json.dumps(diagnostic, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if mismatches or unknown_codes or diagnostic["native_total"] != 26530 or diagnostic["contemporary_total"] != 41083:
        raise RuntimeError("Source-panel mismatch. Inspect reconstruction_audit.json before interpreting Q1.")

    region_rows = []
    for code in sorted(known_codes):
        n, c = native_by_region[code], contemporary_by_region[code]
        g = "island" if code in STRICT_ISLAND_CODES else ("mixed" if code in MIXED_CODES else "mainland")
        region_rows.append({
            "area_code_l3": code, "region_name": names[code], "level1": level1[code],
            "class": g, "native": n, "contemporary": c, "added": c - n,
            "added_share": None if c == 0 else (c - n) / c,
            "increase_pct": None if n == 0 else 100 * (c - n) / n,
        })

    def write_csv(path, values, fields):
        with path.open("w", encoding="utf-8", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields)
            w.writeheader()
            w.writerows(values)

    write_csv(a.out_dir / "region_metrics.csv", region_rows, list(region_rows[0].keys()))
    write_csv(a.out_dir / "added_species_region_pairs.csv", added_pairs, list(added_pairs[0].keys()))
    write_csv(a.out_dir / "species_reconstruction.csv", per_species, list(per_species[0].keys()))
    groups = {label: metrics_for([r for r in region_rows if r["class"] == label])
              for label in ("island", "mainland", "mixed")}
    active = [r for r in region_rows if r["contemporary"] > 0 and r["class"] != "mixed"]
    observed, stratum = adjusted_mean_diff(active)
    label_source = [r["class"] for r in active]
    by_l1 = defaultdict(list)
    for i, r in enumerate(active):
        by_l1[r["level1"]].append(i)
    rng = random.Random(20261102)
    draws = []
    for _ in range(4999):
        labels = label_source.copy()
        for inds in by_l1.values():
            shuffled = [labels[i] for i in inds]
            rng.shuffle(shuffled)
            for i, new_label in zip(inds, shuffled):
                labels[i] = new_label
        value, _ = adjusted_mean_diff(active, labels)
        if value is not None:
            draws.append(value)
    p_two = None if observed is None else (1 + sum(abs(x) >= abs(observed)
                   for x in draws)) / (1 + len(draws))
    outcome = {
        "schema": "chocho_wgsrpd3_strict_island_mainland_resource_q1_v0.2",
        "status": "RECONSTRUCTED_FROM_PINNED_WCVP_HOSTS",
        "design": "Descriptive binary butterfly x WGSRPD3 potential host-resource membership, not realized occupancy",
        "source_commits": {"wcvp": "65bed76bae9d644ccb6ad200c05f9f5071d89e05",
                           "hosts": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
                           "wgsrpd": "52da7828aba9d461dd133c27b3bd7a4407161f54"},
        "classification": {"corrected_named_island": "CNY Canary Is.; RHO Rhode I. is a US mainland state",
                           "official_name_island_test_count": len(explicit_named_islands),
                           "strict_island_codes": sorted(STRICT_ISLAND_CODES),
                           "excluded_mixed_codes": sorted(MIXED_CODES),
                           "mainland": "remaining Level-3 regions, excluding mixed"},
        "panel_total": diagnostic,
        "groups": groups,
        "stratified_contrast": {
            "metric": "island minus mainland mean (introduced_added / contemporary) among active WGSRPD3 regions",
            "strata": stratum,
            "weight": "min(island_regions, mainland_regions) per Level-1 area",
            "weighted_difference": observed,
            "descriptive_within_level1_label_permutation_p_two_sided": p_two,
            "permutations": len(draws),
            "warning": "No causal inference; region areas, species composition and macro-geographic spatial dependence remain uncontrolled."
        },
        "limitations": [
            "WGSRPD3 regions have highly unequal geographic areas.",
            "Islands without a reconstructed native resource baseline have undefined percentage growth; they remain in the added-share analysis.",
            "Positive island resource change does not demonstrate transport, port effects, colonization, or extant butterfly populations.",
            "Historical colonial exposure, ports and native origins of butterflies are not observed in this Q1 analysis."
        ],
    }
    (a.out_dir / "island_mainland_summary.json").write_text(
        json.dumps(outcome, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": outcome["status"], "groups": groups,
                      "contrast": outcome["stratified_contrast"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

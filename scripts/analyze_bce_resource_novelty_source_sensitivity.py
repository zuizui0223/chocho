#!/usr/bin/env python3
"""BCE/Clarke host-inventory sensitivity of geographic resource *novelty*.

Inputs are the exact evidence-filtered European candidate host pairs already
reported in the GEB supplementary Table S9. This is a post-hoc, source-
conditional *descriptive* sensitivity, NOT a new geographic null model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import median


SOURCE_SCHEMA = "chocho_bce_clarke92_wcvp_accepted_id_gap_v0.1"
SOURCE_SHA256 = "f46d18703fe5cc8d97699af6f084f652a704e46b0f9e8c68e32dd70decc32c53"
MATCH_STATUS = "ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY"
EXPECTED_BASE = {
    "focal_butterflies": 239,
    "species_with_expansion": 206,
    "native_cells": 26530,
    "contemporary_cells": 41083,
    "distinct_added_cells": 14553,
    "candidate_incidences": 75346,
    "native_absorbed_incidences": 46387,
    "redundant_new_cell_incidences": 14406,
}
EXPECTED_AUGMENTED = {
    "focal_butterflies": 239,
    "species_with_expansion": 213,
    "native_cells": 29852,
    "contemporary_cells": 44819,
    "distinct_added_cells": 14967,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def load_candidate_host_links(path: Path) -> tuple[dict[str, set[str]], dict[str, object]]:
    sha = digest(path)
    if sha != SOURCE_SHA256:
        raise ValueError("BCE source candidate IDs differ from verified archived bytes")
    obj = json.loads(path.read_text(encoding="utf-8"))
    if obj.get("schema") != SOURCE_SCHEMA or obj.get("exact_focal_butterflies") != 92:
        raise ValueError("BCE accepted-ID crosswalk schema/cohort changed")
    if obj.get("frozen_HOSTS") != "808e0b869f9ec1adf8efff87cf6a395adda103e0":
        raise ValueError("frozen HOSTS identity changed")
    if obj.get("frozen_WCVP") != "65bed76bae9d644ccb6ad200c05f9f5071d89e05":
        raise ValueError("frozen WCVP identity changed")

    additions: dict[str, set[str]] = defaultdict(set)
    matched = []
    nrows = 0
    raw_links = 0
    for species_record in obj["per_butterfly"]:
        species = str(species_record["species"]).strip()
        if not species or species in matched:
            raise ValueError("invalid or duplicate BCE species")
        matched.append(species)
        for entry in species_record["details"]:
            if entry["match_status"] != MATCH_STATUS:
                continue
            plant = str(entry["accepted_plant_name_id"]).strip()
            if not plant:
                raise ValueError("BCE match without accepted WCVP ID")
            nrows += 1
            additions[species].add(plant)
        raw_links += len(species_record["details"])
    nlinks = sum(len(x) for x in additions.values())
    unique_plants = len(set().union(*additions.values()))
    if (
        len(matched) != 92 or len(additions) != 81
        or nlinks != 1027 or nrows != 1027 or unique_plants != 679
    ):
        raise ValueError(
            f"archived source scope drift 92/81/1027/679: "
            f"{len(matched)}/{len(additions)}/{nlinks}/{unique_plants}"
        )
    meta = {
        "bce_accepted_edges_sha256": sha,
        "source_schema": SOURCE_SCHEMA,
        "source_exact_butterfly_matches": len(matched),
        "butterflies_with_candidate_hosts": len(additions),
        "candidate_additional_host_links": nlinks,
        "unique_additional_accepted_hosts": unique_plants,
        "evidence": "BCE-curated Clarke (2024), evidence levels 1-3; not independent field feeding observations",
        "source_workflow_run": 37795847333,
        "source_workflow_artifact_id": 11557274784,
    }
    return additions, {**meta, "BCE_exact_butterflies": sorted(matched)}


def host_geography(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for row in rows(path):
        host = row["accepted_plant_name_id"].strip()
        reg = row["area_code_l3"].strip()
        if host and reg:
            out[host].add(reg)
    return dict(out)


def combine_host_ranges(
    base: dict[str, set[str]],
    add: dict[str, set[str]],
    extra_hosts: set[str],
) -> dict[str, set[str]]:
    result = {k: set(v) for k, v in base.items()}
    for host in extra_hosts:
        r = add.get(host, set())
        if host in result:
            if result[host] != r:
                raise ValueError(f"WCVP range mismatch for overlapping accepted ID {host}")
        else:
            result[host] = set(r)
    return result


def count_species_opportunity(
    ids: set[str], native: dict[str, set[str]], contemporary: dict[str, set[str]]
) -> dict[str, object]:
    if not ids:
        raise ValueError("focal butterfly has no accepted host IDs")
    native_union: set[str] = set()
    contemporary_union: set[str] = set()
    introduced_per_host: dict[str, set[str]] = {}
    for host in ids:
        n = native.get(host, set())
        c = contemporary.get(host, set())
        if not n <= c:
            raise ValueError(f"native WCVP geography outside contemporary for {host}")
        native_union.update(n)
        contemporary_union.update(c)
        introduced_per_host[host] = c - n

    H = sum(len(u) for u in introduced_per_host.values())
    A = sum(len(u & native_union) for u in introduced_per_host.values())
    U = len(contemporary_union - native_union)
    R = H - A - U
    if R < 0 or A > H or not native_union <= contemporary_union:
        raise ValueError("resource-incidence identity failed")
    return {
        "native": len(native_union),
        "contemporary": len(contemporary_union),
        "added": U,
        "H": H,
        "A": A,
        "R": R,
        "native_regions": native_union,
        "contemporary_regions": contemporary_union,
        "added_regions": contemporary_union - native_union,
    }


def summation(metrics: list[dict[str, object]]) -> dict[str, int | float | None]:
    n = len(metrics)
    totals = {
        "focal_butterflies": n,
        "species_with_expansion": sum(int(x["added"]) > 0 for x in metrics),
        "native_cells": sum(int(x["native"]) for x in metrics),
        "contemporary_cells": sum(int(x["contemporary"]) for x in metrics),
        "distinct_added_cells": sum(int(x["added"]) for x in metrics),
        "candidate_incidences": sum(int(x["H"]) for x in metrics),
        "native_absorbed_incidences": sum(int(x["A"]) for x in metrics),
        "redundant_new_cell_incidences": sum(int(x["R"]) for x in metrics),
    }
    H = totals["candidate_incidences"]
    if H != (
        totals["native_absorbed_incidences"]
        + totals["redundant_new_cell_incidences"]
        + totals["distinct_added_cells"]
    ):
        raise ValueError("global H=A+R+U identity failed")
    totals["native_absorption_fraction"] = totals["native_absorbed_incidences"]/H if H else None
    totals["distinct_new_region_per_candidate_fraction"] = totals["distinct_added_cells"]/H if H else None
    totals["mean_introduced_support_per_new_cell"] = (
        (H-totals["native_absorbed_incidences"]) / totals["distinct_added_cells"]
        if totals["distinct_added_cells"] else None
    )
    return totals


def analyze(
    candidate: dict[str, set[str]],
    source: dict[str, object],
    metrics: list[dict[str, str]],
    baseline_links: list[dict[str, str]],
    baseline_native: dict[str, set[str]],
    baseline_contemporary: dict[str, set[str]],
    additional_native: dict[str, set[str]],
    additional_contemporary: dict[str, set[str]],
) -> dict[str, object]:
    metric_by_sp = {x["species"]: x for x in metrics}
    if len(metric_by_sp) != len(metrics) or len(metric_by_sp) != 239:
        raise ValueError("focal 239 butterfly metrics missing or duplicated")
    hosts_by_sp: dict[str, set[str]] = defaultdict(set)
    for row in baseline_links:
        sp = row["insect_species"].strip()
        host = row["accepted_plant_name_id"].strip()
        if sp in metric_by_sp and host:
            hosts_by_sp[sp].add(host)
    if any(not hosts_by_sp[sp] for sp in metric_by_sp):
        raise ValueError("missing original accepted hosts for focal butterfly")
    if not set(candidate).issubset(metric_by_sp):
        raise ValueError("candidate butterfly not in original focal 239")
    if any(candidate[sp] & hosts_by_sp[sp] for sp in candidate):
        raise ValueError("candidate crosswalk no longer absent from original HOSTS species")
    all_extra = set().union(*candidate.values())
    native = combine_host_ranges(baseline_native, additional_native, all_extra)
    cont = combine_host_ranges(baseline_contemporary, additional_contemporary, all_extra)

    base_species = {}
    augmented_species = {}
    by_sp = []
    source_species = set(source["BCE_exact_butterflies"])
    for sp in sorted(metric_by_sp):
        ids = hosts_by_sp[sp]
        b = count_species_opportunity(ids, baseline_native, baseline_contemporary)
        a = count_species_opportunity(ids | candidate.get(sp, set()), native, cont)
        expected = metric_by_sp[sp]
        if (b["native"], b["contemporary"], b["added"]) != (
            int(expected["native_resource_units"]),
            int(expected["contemporary_resource_units"]),
            int(expected["introduced_added_units"]),
        ):
            raise ValueError(f"original 239-species footprint drift: {sp}")
        base_species[sp] = b
        augmented_species[sp] = a
        old_gain = b["added_regions"]
        new_gain = a["added_regions"]
        gained_native = a["native_regions"]-b["native_regions"]
        # A cell can cease being introduced-added because it is now natively
        # supported, even though contemporary footprint only increases.
        lost_new = old_gain-new_gain
        if not lost_new <= gained_native:
            raise ValueError(f"unexpected lost new region not due to native reclassification: {sp}")
        added_new = new_gain-old_gain
        rec = {
            "species": sp,
            "BCE_candidate_host_links_added": len(candidate.get(sp, set())),
            "BCE_matched": sp in source_species,
            "base_native": b["native"],
            "base_added": b["added"],
            "augmented_native": a["native"],
            "augmented_added": a["added"],
            "new_native_regions": len(gained_native),
            "old_added_regions_reclassified_as_native": len(lost_new),
            "previously_unadded_regions_newly_added": len(added_new),
            "base_H": b["H"], "base_A": b["A"], "base_R": b["R"],
            "augmented_H": a["H"], "augmented_A": a["A"], "augmented_R": a["R"],
        }
        by_sp.append(rec)
    base_total = summation(list(base_species.values()))
    aug_total = summation(list(augmented_species.values()))
    for key, count in EXPECTED_BASE.items():
        if base_total[key] != count:
            raise ValueError(f"original baseline differs from frozen result {key}: {base_total[key]} != {count}")
    for key, count in EXPECTED_AUGMENTED.items():
        if aug_total[key] != count:
            raise ValueError(f"archived BCE union geography disagrees for {key}: {aug_total[key]} != {count}")

    def subset_aggregates(subset: set[str]) -> dict[str, object]:
        b = summation([base_species[sp] for sp in sorted(subset)])
        a = summation([augmented_species[sp] for sp in sorted(subset)])
        return {
            "species": len(subset),
            "base": b,
            "BCE_augmented": a,
            "native_absorption_share_change_pp": 100*(
                a["native_absorption_fraction"] - b["native_absorption_fraction"]
            ) if a["native_absorption_fraction"] is not None
              and b["native_absorption_fraction"] is not None else None,
        }
    changed_species = set(candidate)
    unaffected = set(metric_by_sp)-changed_species
    matched = set(source_species)
    output = {
        "schema":"chocho_butterfly_bce_host_inventory_resource_novelty_v0.1",
        "status":"SUCCESS_POSTHOC_EUROPE_BIASED_INVENTORY_SENSITIVITY_NOT_SPATIAL_NULL",
        "source": {
            **{k:v for k,v in source.items() if k!="BCE_exact_butterflies"},
            "primary_resource_baseline": "frozen HOSTS/WCVP 239 butterflies; 365 recorded resource regions",
            "WCVP": "65bed76bae9d644ccb6ad200c05f9f5071d89e05",
            "HOSTS": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
        },
        "base": base_total,
        "BCE_augmented_partial": aug_total,
        "change": {
            "delta_native_absorption_fraction_percentage_points":
                100*(aug_total["native_absorption_fraction"]-base_total["native_absorption_fraction"]),
            "delta_added_fraction_of_H_percentage_points":
                100*(aug_total["distinct_new_region_per_candidate_fraction"]-base_total["distinct_new_region_per_candidate_fraction"]),
            "delta_A": aug_total["native_absorbed_incidences"]-base_total["native_absorbed_incidences"],
            "delta_H": aug_total["candidate_incidences"]-base_total["candidate_incidences"],
            "delta_R": aug_total["redundant_new_cell_incidences"]-base_total["redundant_new_cell_incidences"],
            "delta_U": aug_total["distinct_added_cells"]-base_total["distinct_added_cells"],
            "old_added_cells_reclassified_native":sum(x["old_added_regions_reclassified_as_native"] for x in by_sp),
            "previously_unadded_regions_newly_added":sum(x["previously_unadded_regions_newly_added"] for x in by_sp),
        },
        "subgroups": {
            "81_BCE_candidate_augmented_butterflies":subset_aggregates(changed_species),
            "92_BCE_exact_binomial_matches":subset_aggregates(matched),
            "158_unmodified_host_inventory_butterflies":subset_aggregates(unaffected),
        },
        "affected_butterflies_sorted": [
            x for x in by_sp if x["BCE_candidate_host_links_added"] > 0
        ],
        "scientific_boundaries":[
            "BCE/Clarke source is Europe-biased: not a second globally representative HOSTS dataset.",
            "All candidate host associations have been projected across accepted WCVP ranges; they are not verified local larval feeding.",
            "Even augmented A/H is a descriptive fraction, not a statistical expectation under the new host inventory.",
            "The +271–276 original numerically qualified 1500-km spatial-null deviation CANNOT be transferred to BCE augmented hosts; that would require a scenario-specific new randomization and full mixing diagnostics.",
            "Candidate host additions can reclassify previously introduced-supported regions as natively supported; interpretation requires denominator labels.",
            "No causal conclusions about introduction pathways, host preference, recruitment, or butterfly conservation follow from these maps.",
        ],
    }
    return output


def main() -> None:
    p = argparse.ArgumentParser()
    for key in (
        "bce_accepted_edges_json", "frozen_species_metrics_csv",
        "original_insect_host_csv", "original_native_csv",
        "original_contemporary_csv", "bce_additional_native_csv",
        "bce_additional_contemporary_csv", "output_json",
    ):
        p.add_argument("--"+key.replace("_","-"), type=Path, required=True)
    a = p.parse_args()
    added, source=load_candidate_host_links(a.bce_accepted_edges_json)
    res=analyze(
        added, source, rows(a.frozen_species_metrics_csv),
        rows(a.original_insect_host_csv),
        host_geography(a.original_native_csv),
        host_geography(a.original_contemporary_csv),
        host_geography(a.bce_additional_native_csv),
        host_geography(a.bce_additional_contemporary_csv),
    )
    a.output_json.parent.mkdir(parents=True, exist_ok=True)
    a.output_json.write_text(
        json.dumps(res, indent=2, ensure_ascii=False)+"\n",encoding="utf-8"
    )
    print(json.dumps({
        "schema":res["schema"],
        "status":res["status"],
        "source":res["source"],
        "base":res["base"],
        "BCE_augmented_partial":res["BCE_augmented_partial"],
        "change":res["change"],
        "subgroups":res["subgroups"],
        "output":str(a.output_json),
    },indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()

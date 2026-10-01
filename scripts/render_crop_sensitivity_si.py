#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def fmt(value, digits=3):
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.{digits}f}"


def percent(value, digits=1):
    if value is None:
        return "NA"
    return f"{100.0 * float(value):.{digits}f}%"


def expansion_row(label: str, block: dict) -> str:
    return " | ".join([
        label,
        str(block["species_with_non_crop_native_resource"]),
        f'{block["species_expanded"]} ({percent(block["fraction_evaluable_species_expanded"])})',
        f'{block["total_native_species_units"]:,}',
        f'{block["total_contemporary_species_units"]:,}',
        f'{block["total_added_species_units"]:,}',
        percent(block["aggregate_proportional_increase"]),
        fmt(block["spearman_host_family_vs_log_resource_expansion"]),
    ])


def architecture_row(label: str, block: dict) -> str:
    conc = block["aggregate_host_contribution_concentration"]
    return " | ".join([
        label,
        str(block["portfolio_expanded_species"]),
        str(conc["hosts_for_half_added_credit"]),
        fmt(block["median_maximum_single_host_fractional_share"]),
        fmt(block["median_effective_contributor_number"]),
        fmt(block["spearman_host_family_vs_effective_contributor_number"]),
        fmt(block["spearman_host_family_vs_maximum_single_host_fractional_share"]),
    ])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-json", type=Path, required=True)
    ap.add_argument("--crop-json", type=Path, required=True)
    ap.add_argument("--output-md", type=Path, required=True)
    args = ap.parse_args()

    result = json.loads(args.result_json.read_text(encoding="utf-8"))
    crop = json.loads(args.crop_json.read_text(encoding="utf-8"))
    if result.get("schema") != "chocho_butterfly_crop_exclusion_sensitivity_result_v0.1":
        raise RuntimeError("unexpected crop-sensitivity result")
    if result.get("status") != "CROP_EXCLUSION_SENSITIVITY_COMPLETE":
        raise RuntimeError("crop-sensitivity result incomplete")

    baseline_full = result["baseline_reconstruction"]["full_resource_eligible_panel"]
    baseline_adequate = result["baseline_reconstruction"][
        "host_taxonomy_lower_bound_adequate_subset"
    ]
    exact = result["variants"]["fao_exact_binomial"]
    conservative = result["variants"]["fao_exact_plus_spp_genus"]

    exact_full = exact["full_resource_eligible_panel"]
    conservative_full = conservative["full_resource_eligible_panel"]
    exact_adequate = exact["host_taxonomy_lower_bound_adequate_subset"]
    conservative_adequate = conservative[
        "host_taxonomy_lower_bound_adequate_subset"
    ]

    exact_same_direction = (
        exact_full["spearman_host_family_vs_log_resource_expansion"] is not None
        and exact_adequate[
            "spearman_host_family_vs_effective_contributor_number"
        ] is not None
        and exact_adequate[
            "spearman_host_family_vs_maximum_single_host_fractional_share"
        ] is not None
        and exact_adequate[
            "spearman_host_family_vs_effective_contributor_number"
        ] > 0
        and exact_adequate[
            "spearman_host_family_vs_maximum_single_host_fractional_share"
        ] < 0
    )

    if exact_same_direction:
        interpretation = (
            "Excluding FAO-defined crop hosts did not change the qualitative "
            "conclusions: resource expansion remained widespread, host-family "
            "breadth remained only weakly associated with proportional expansion, "
            "and broader diets retained a more distributed host-contribution "
            "architecture."
        )
    else:
        interpretation = (
            "Excluding FAO-defined crop hosts materially altered at least one "
            "central pattern; the corresponding primary claim should therefore be "
            "qualified rather than described as crop-robust."
        )

    lines = [
        "# Supporting Information — Crop-host sensitivity",
        "",
        "WCVP introduced distributions represent wild non-native distributions rather "
        "than cultivated acreage. As a sensitivity analysis, we additionally removed "
        "host taxa classified as crops by the FAO Indicative Crop Classification "
        "before reconstructing butterfly-level native and contemporary resource "
        "envelopes. The primary crop definition used exact botanical binomial matches; "
        "an intentionally over-conservative variant also removed all congeners when "
        "the FAO list specified a crop at genus level (spp.). Butterfly host-family "
        "breadth was kept fixed as the specialization predictor.",
        "",
        f"The FAO parser identified {len(crop['exact_binomials'])} exact crop binomials "
        f"and {len(crop['genus_wildcards'])} genus-level crop entries. In the butterfly "
        f"host reconstruction, the exact rule excluded {exact['crop_host_ids_excluded']} "
        f"resolved host taxa affecting {exact['butterflies_with_at_least_one_crop_host_removed']} "
        f"butterfly species; the conservative rule excluded "
        f"{conservative['crop_host_ids_excluded']} host taxa affecting "
        f"{conservative['butterflies_with_at_least_one_crop_host_removed']} species.",
        "",
        "## Table S1. Sensitivity of resource expansion and host-contribution architecture to crop-host exclusion",
        "",
        "**A. Resource expansion in the original 239-species resource-eligible panel**",
        "",
        "Analysis | Evaluable species | Expanded species | Native species-units | Contemporary species-units | Added species-units | Aggregate increase | Spearman rho (host-family breadth vs log expansion)",
        "--- | ---: | ---: | ---: | ---: | ---: | ---: | ---:",
        expansion_row("Primary reconstruction", baseline_full),
        expansion_row("FAO exact-binomial crop exclusion", exact_full),
        expansion_row("FAO exact + genus-level conservative exclusion", conservative_full),
        "",
        "**B. Host-contribution architecture in the host-taxonomy-adequate subset**",
        "",
        "Analysis | Expanded species | Host taxa for 50% of added credit | Median maximum single-host share | Median effective contributors | Spearman rho (breadth vs effective contributors) | Spearman rho (breadth vs maximum share)",
        "--- | ---: | ---: | ---: | ---: | ---: | ---:",
        architecture_row("Primary reconstruction", baseline_adequate),
        architecture_row("FAO exact-binomial crop exclusion", exact_adequate),
        architecture_row(
            "FAO exact + genus-level conservative exclusion",
            conservative_adequate,
        ),
        "",
        interpretation,
        "",
        "The conservative genus-expanded analysis deliberately removes non-crop "
        "congeners in genera represented as spp. in the FAO source and is therefore "
        "treated as a stress test rather than an alternative primary definition.",
        "",
    ]

    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({
        "output": str(args.output_md),
        "exact_same_direction": exact_same_direction,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

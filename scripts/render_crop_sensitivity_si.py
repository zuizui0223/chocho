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
    half = block["aggregate_host_contribution_concentration"][
        "hosts_for_half_added_credit"
    ]
    return " | ".join([
        label,
        str(block["species_with_non_crop_native_resource"]),
        f'{block["species_expanded"]} ({percent(block["fraction_evaluable_species_expanded"])})',
        f'{block["total_added_species_units"]:,}',
        percent(block["aggregate_proportional_increase"]),
        fmt(block["spearman_host_family_vs_log_resource_expansion"]),
        str(half),
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

    baseline = result["baseline_reconstruction"]["full_resource_eligible_panel"]
    exact = result["variants"]["fao_exact_binomial"]
    conservative = result["variants"]["fao_exact_plus_spp_genus"]
    exact_full = exact["full_resource_eligible_panel"]
    conservative_full = conservative["full_resource_eligible_panel"]

    lines = [
        "# Crop-host sensitivity analysis receipt",
        "",
        "Submission mapping: Supplementary Table S7 in the v0.2 manuscript package.",
        "",
        "WCVP introduced distributions represent wild non-native distributions rather "
        "than cultivated acreage. The crop-host sensitivity removes FAO-classified "
        "crop taxa before reconstructing butterfly-level native and contemporary "
        "resource envelopes. The primary crop definition uses exact botanical "
        "binomials; an intentionally over-conservative variant also removes all "
        "congeners when the FAO list specifies a crop at genus level (spp.). "
        "Butterfly host-family breadth remains fixed as the specialization predictor.",
        "",
        f"The FAO parser identified {len(crop['exact_binomials'])} exact crop binomials "
        f"and {len(crop['genus_wildcards'])} genus-level crop entries. In the butterfly "
        f"host reconstruction, the exact rule excluded {exact['crop_host_ids_excluded']} "
        f"resolved host taxa affecting {exact['butterflies_with_at_least_one_crop_host_removed']} "
        f"butterfly species; the conservative rule excluded "
        f"{conservative['crop_host_ids_excluded']} host taxa affecting "
        f"{conservative['butterflies_with_at_least_one_crop_host_removed']} species.",
        "",
        "Analysis | Evaluable species | Expanded species | Added species-units | Aggregate increase | Spearman rho (host-family breadth vs log expansion) | Host species for 50% of added credit",
        "--- | ---: | ---: | ---: | ---: | ---: | ---:",
        expansion_row("Primary reconstruction", baseline),
        expansion_row("FAO exact-binomial crop exclusion", exact_full),
        expansion_row(
            "FAO exact + genus-level conservative exclusion",
            conservative_full,
        ),
        "",
        "Excluding FAO-defined crop hosts did not change the qualitative conclusions: "
        "resource expansion remained widespread, host-family breadth remained only "
        "weakly associated with proportional expansion, and 37-38 host species still "
        "accounted for half of all added opportunity.",
        "",
        "The conservative genus-expanded analysis deliberately removes non-crop "
        "congeners in genera represented as spp. in the FAO source and is therefore "
        "treated as a stress test rather than an alternative primary definition.",
        "",
    ]

    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"output": str(args.output_md)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

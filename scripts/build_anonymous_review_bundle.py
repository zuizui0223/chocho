#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

COPY_PATHS = (
    "pyproject.toml",
    "src",
    "scripts",
    "tests",
    "data/external",
    "data/frozen",
    "docs/exploratory",
    "docs/REPRODUCIBILITY.md",
    "benchmarks/exploratory",
    "provenance/reviewer_defenses",
    "provenance/climate",
    "provenance/crop",
    "provenance/archive/climate/operations",
    "provenance/archive/climate/pilot",
    "manuscript/generated/butterfly_specialization_ecology_blinded_v0.2.md",
    "manuscript/butterfly_specialization_claim_map_v0.2.json",
    "manuscript/butterfly_specialization_supplement_v0.2.md",
)

# Independent analyses from draft PR #43 and its continuations are excluded from the
# EXISTING GEB v0.2 double-anonymous reviewer snapshot. They are delivered via
# their separate analysis branches and review packages, never the submitted pack.
# Review the list whenever an independent analysis gains further files.
EXCLUDE_RELATIVE = {
    Path("docs/exploratory/BUTTERFLY_GLOBAL_RESOURCE_NOVELTY_PAPER_STORY_20261016.md"),
    Path("docs/exploratory/BUTTERFLY_RESOURCE_NOVELTY_FULL_DRAFT_20261017.md"),
    Path("docs/exploratory/NONTRIVIAL_MACROECOLOGICAL_FINDINGS_20261015.md"),
    Path("provenance/reviewer_defenses/ALTERNATIVE_RESOURCE_NOVELTY_AUTHOR_DECISION_20261017.md"),
    Path("provenance/reviewer_defenses/protocols/BCE_CLARKE_RESOURCE_NOVELTY_SENSITIVITY_20261017.md"),
    Path("provenance/reviewer_defenses/protocols/CLIMATE_CONDITIONED_RESOURCE_ABSORPTION_20261015.md"),
    Path("provenance/reviewer_defenses/protocols/SPATIALLY_LOCAL_RESOURCE_ABSORPTION_20261015.md"),
    Path("provenance/reviewer_defenses/protocols/SPATIAL_NULL_EXTENDED_MIXING_GATE_20261015.md"),
    Path("provenance/reviewer_defenses/results/butterfly_added_cell_peer_null_v0.1.json"),
    Path("provenance/reviewer_defenses/results/butterfly_bce_added_region_identity_turnover_v0.1.json"),
    Path("provenance/reviewer_defenses/results/butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json"),
    Path("provenance/reviewer_defenses/results/butterfly_climate_matched_absorption_v0.1.json"),
    Path("provenance/reviewer_defenses/results/butterfly_exact_host_pairwise_overlap_v0.1.json"),
    Path("provenance/reviewer_defenses/results/butterfly_plant_range_structural_null_two_chains_v0.2.json"),
    Path("provenance/reviewer_defenses/results/butterfly_plant_range_structural_null_v0.1.json"),
    Path("provenance/reviewer_defenses/results/butterfly_spatial_null_longchain_qualified_v0.2.json"),
    Path("provenance/reviewer_defenses/results/butterfly_spatial_null_original_mixing_audit_v0.1.json"),
    Path("scripts/analyze_bce_resource_novelty_source_sensitivity.py"),
    Path("scripts/analyze_butterfly_climate_matched_plant_range_null.py"),
    Path("scripts/analyze_butterfly_shared_host_structural_null.py"),
    Path("scripts/analyze_butterfly_spatially_local_plant_range_null.py"),
    Path("scripts/analyze_butterfly_species_redundancy.py"),
    Path("scripts/audit_butterfly_spatial_null_longchain.py"),
    Path("scripts/build_bce_additional_wcvp_ranges.R"),
    Path("scripts/build_butterfly_resource_novelty_review_bundle.py"),
    Path("scripts/render_bce_added_resource_identity_turnover_figure.py"),
    Path("scripts/render_bce_resource_novelty_sensitivity_figure.py"),
    Path("scripts/render_butterfly_resource_ecological_context_figure.py"),
    Path("scripts/render_butterfly_resource_novelty_figure.py"),
    Path("scripts/summarize_bce_added_resource_identity_turnover.py"),
    Path("tests/test_bce_added_resource_identity_turnover.py"),
    Path("tests/test_bce_resource_novelty_sensitivity_figure.py"),
    Path("tests/test_bce_resource_novelty_source_sensitivity.py"),
    Path("tests/test_butterfly_climate_matched_plant_range_null.py"),
    Path("tests/test_butterfly_resource_ecological_context_figure.py"),
    Path("tests/test_butterfly_resource_novelty_figure.py"),
    Path("tests/test_butterfly_resource_novelty_full_draft.py"),
    Path("tests/test_butterfly_shared_host_structural_null.py"),
    Path("tests/test_butterfly_spatial_null_longchain_audit.py"),
    Path("tests/test_butterfly_spatially_local_plant_range_null.py"),
    Path("tests/test_butterfly_species_redundancy.py"),
    Path("tests/test_spatially_local_range_null_toy_mixing.py"),
    Path("scripts/build_anonymous_review_bundle.py"),
    Path("scripts/release_preflight.py"),
    Path("scripts/submission_preflight.py"),
    Path("tests/test_release_metadata.py"),
    Path("tests/test_reproducibility_snapshot_integrity.py"),
    Path("tests/test_butterfly_specialization_submission_bundle.py"),
    Path("tests/test_butterfly_specialization_manuscript_bundle.py"),
    Path("tests/test_butterfly_climate_release_postgate_recovery.py"),
    Path("tests/test_butterfly_climate_release_execution_binding.py"),
    # Exploratory host-climate point studies are outside the frozen GEB review pack.
    Path("tests/test_butterfly_73_point_entry_roster_20261027.py"),
    Path("tests/test_butterfly_73_point_tdwg_resource_20261027.py"),
    Path("tests/test_butterfly_73_worldclim_point_niche_20261027.py"),
    Path("tests/test_butterfly_niche_entry_country_20261026.py"),
    Path("tests/test_butterfly_host_point_proximity_20261026.py"),
    Path("scripts/analyze_butterfly_73_country_climate_area_null_20261028.py"),
    Path("docs/exploratory/BUTTERFLY_73_COUNTRY_LAND_CLIMATE_NULL_RESULTS_20261028.md"),
    Path("provenance/reviewer_defenses/results/butterfly_73_country_land_climate_null_20261028.json"),
    Path("scripts/audit_europe_historical_butterfly_gbif_source_20261028.py"),
    Path("scripts/audit_ukbms_historical_event_source_20261028.py"),
    Path("scripts/audit_uk_public_transect_source_20261028.py"),
    Path("scripts/audit_uk_ltmn_xls_schema_20261028.py"),
    Path("scripts/audit_uk_ltmn_historical_event_join_20261028.py"),
    Path("scripts/audit_uk_ltmn_site_crosswalk_20261028.py"),
    Path("scripts/audit_uk_ltmn_site_geography_20261028.py"),
    Path("scripts/audit_uk_ltmn_butterfly_detection_panel_20261028.py"),
    Path("tests/test_uk_ltmn_historical_event_join_20261028.py"),
    Path("tests/test_uk_ltmn_site_geography_20261028.py"),
    Path("tests/test_uk_ltmn_butterfly_detection_panel_20261028.py"),
    Path("scripts/probe_uk_era5_historical_weather_20261028.py"),
    Path("scripts/audit_uk_ltmn_ERA5_historical_baseline_20261028.py"),
    Path("tests/test_uk_ltmn_ERA5_historical_baseline_20261028.py"),
    Path("scripts/probe_uk_nasa_power_1981_2021_20261028.py"),
    Path("tests/test_uk_nasa_power_1981_2021_20261028.py"),
    Path("scripts/audit_uk_ltmn_NASA_POWER_historical_baseline_20261028.py"),
    Path("tests/test_uk_ltmn_NASA_POWER_historical_baseline_20261028.py"),
    Path("tests/test_uk_public_transect_source_20261028.py"),
    Path("tests/test_ukbms_historical_frozen_result_20261028.py"),
    Path("tests/test_ukbms_historical_event_source_20261028.py"),
    Path("tests/test_europe_historical_facets_frozen_receipt_20261028.py"),
    Path("tests/test_europe_historical_butterfly_gbif_source_20261028.py"),
    Path("tests/test_butterfly_73_country_climate_area_null_20261028.py"),
    Path("provenance/reviewer_defenses/protocols/BUTTERFLY_73_COUNTRY_LAND_CLIMATE_NULL_20261028.md"),
    Path("tests/test_butterfly_73_paired_point_climate_20261028.py"),
    Path("scripts/analyze_butterfly_73_paired_point_climate_20261028.py"),
    Path("data/exploratory/butterfly_73_worldclim_case_counts_20261027.csv"),
    Path("docs/exploratory/BUTTERFLY_73_PAIRED_POINT_CLIMATE_ATTRITION_20261028.md"),
    Path("provenance/reviewer_defenses/results/butterfly_73_paired_point_climate_20261028.json"),
    Path("tests/test_repository_layout.py"),
    Path("provenance/reviewer_defenses/GEB_V02_DECISION_MEMO.md"),
    Path("provenance/reviewer_defenses/REVIEW_RESPONSE_MAP.md"),
    # Country-grain continuation: its data and dependencies are outside this snapshot.
    Path("docs/exploratory/GBIF_COUNTRY_MACRO_READOUT_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_COUNTRY_FULL_INDEX_MACRO_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_COUNTRY_RECORD_TURNOVER_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_RECOVERY_CONTINUATION_20261020.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/country_design_audit.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/country_resource_area_audit.csv.gz"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_macro_results.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_macro_summary.png"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_macro_summary.svg"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_sensitivities.csv"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/original133_alias_exclusion_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/original135_replay_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/verification_20261020.json"),
    Path("scripts/acquire_gbif_country_facets_20261020.py"),
    Path("scripts/analyze_gbif_country_host_opportunity_20261020.py"),
    Path("scripts/build_country_area_atoms_20261020.py"),
    Path("scripts/build_country_wgsrpd_crosswalk_20261020.py"),
    Path("scripts/render_gbif_country_macro_20261020.py"),
    Path("scripts/validate_gbif_recovery_taxon_keys_20261020.py"),
    Path("tests/test_country_area_atoms_20261020.py"),
    Path("tests/test_country_wgsrpd_crosswalk_20261020.py"),
    Path("tests/test_gbif_country_facets_20261020.py"),
    Path("tests/test_gbif_country_macro_20261020.py"),
    Path("tests/test_gbif_macro_information_gate_20261020.py"),
    Path("tests/test_gbif_recovery_taxon_keys_20261020.py"),
    Path("scripts/analyze_gbif_country_geographic_influence_20261020.py"),
    Path("scripts/analyze_gbif_country_record_turnover_20261020.py"),
    Path("scripts/audit_gbif_country_record_provenance_20261020.py"),
    Path("scripts/audit_gbif_country_dataset_year_20261020.py"),
    Path("scripts/audit_gbif_country_platform_reporting_20261020.py"),
    Path("scripts/audit_gbif_positive48_subcountry_20261020.py"),
    Path("scripts/summarize_gbif_positive48_subcountry_20261020.py"),
    Path("tests/test_gbif_positive48_subcountry_summary_20261020.py"),
    Path("tests/test_gbif_positive48_subcountry_20261020.py"),
    Path("tests/test_gbif_country_platform_reporting_20261020.py"),
    Path("tests/test_gbif_country_dataset_year_20261020.py"),
    Path("scripts/audit_gbif_dataset_registry_20261020.py"),
    Path("tests/test_gbif_dataset_registry_20261020.py"),
    Path("provenance/reviewer_defenses/protocols/GBIF_DATASET_REGISTRY_FOLLOWUP_20261020.md"),
    Path("tests/test_gbif_country_record_provenance_20261020.py"),
    Path("provenance/reviewer_defenses/protocols/GBIF_COUNTRY_SOURCE_INDEPENDENCE_GATE_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_COUNTRY_DATASET_YEAR_CONFOUND_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_COUNTRY_PLATFORM_EFFORT_DENOM_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_SUBCOUNTRY_POSITIVE48_QC_20261020.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_JUVENILE_GADM1_CONDITIONAL_20261021.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_LOCAL_HOST_BRIDGE_20261021.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_ALL26_INTRODUCED_HOST_INVENTORY_20261023.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_all26_source_host_country_inventory_20261023.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_ALL26_SOURCE_HOST_INVENTORY_READOUT_20261023.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_ALL26_PLANT_GADM_ACQUISITION_20261023.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_all26_host_gadm_source_20261023.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_ALL26_HOST_GADM_READOUT_20261023.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_ALL26_HOST_IDENTITY_SPATIAL_REFERENCE_20261023.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_ALL26_WEIGHTING_SENSITIVITY_20261023.md"),
    Path("provenance/reviewer_defenses/protocols/GBIF_ALL26_LOCAL_HOST_REDUNDANCY_POSTHOC_20261024.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_local_host_redundancy_summary_20261024.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_LOCAL_HOST_REDUNDANCY_READOUT_20261024.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_juvenile_host_annotation_coverage_20261023.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_JUVENILE_HOST_METADATA_READOUT_20261024.md"),
    Path("scripts/analyze_gbif_all26_local_host_redundancy_20261024.py"),
    Path("tests/test_gbif_all26_local_host_redundancy_20261024.py"),
    Path("provenance/reviewer_defenses/protocols/GBIF_JUVENILE_HOST_METADATA_GATE_20261023.md"),
    Path("scripts/audit_gbif_juvenile_host_annotations_20261023.py"),
    Path("tests/test_gbif_juvenile_host_annotations_20261023.py"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_all26_weighting_stability_20261023.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_ALL26_WEIGHTING_STABILITY_READOUT_20261023.md"),
    Path("scripts/analyze_gbif_all26_weighting_stability_20261023.py"),
    Path("tests/test_gbif_all26_weighting_stability_20261023.py"),
    Path("scripts/analyze_gbif_all26_host_identity_reference_20261023.py"),
    Path("tests/test_gbif_all26_host_identity_reference_20261023.py"),
    Path("scripts/audit_gbif_all26_host_gadm_20261023.py"),
    Path("tests/test_gbif_all26_host_gadm_20261023.py"),
    Path("scripts/derive_gbif_all26_source_host_country_inventory_20261023.py"),
    Path("tests/test_gbif_all26_source_host_country_inventory_20261023.py"),
    Path("scripts/audit_gbif_five_host_gadm_cocoverage_20261021.py"),
    Path("scripts/analyze_gbif_five_host_gadm_reference_20261022.py"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_five_host_plant_gadm_source_20261021.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_five_host_gadm_conditioned_reference_20261022.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_FIVE_HOST_GADM_LOCAL_READOUT_20261022.md"),
    Path("tests/test_gbif_five_host_gadm_reference_20261022.py"),
    Path("tests/test_gbif_five_host_gadm_cocoverage_20261021.py"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_source_introduced_hosts_five_cases_20261021.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_SOURCE_INTRODUCED_HOSTS_FIVE_20261021.md"),
    Path("scripts/derive_gbif_positive48_introduced_hosts_20261021.py"),
    Path("tests/test_gbif_positive48_introduced_hosts_20261021.py"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_juvenile_gadm1_conditional_20261021.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_juvenile_gadm1_readout_20261021.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_JUVENILE_GADM1_READOUT_20261021.md"),
    Path("scripts/audit_gbif_juvenile_gadm1_20261021.py"),
    Path("tests/test_gbif_juvenile_gadm1_20261021.py"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_positive48_gadm_lifestage_facets_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_positive48_subcountry_stage_summary_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_POSITIVE48_SUBCOUNTRY_STAGE_READOUT_20261020.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_provenance_facet_qc_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_dataset_registry_qc_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_COUNTRY_SOURCE_PROVENANCE_READOUT_20261020.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_dataset_year_matrix_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_COUNTRY_DATASET_YEAR_READOUT_20261020.md"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_platform_reporting_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_COUNTRY_PLATFORM_REPORTING_READOUT_20261020.md"),
    Path("tests/test_gbif_country_record_turnover_20261020.py"),
    Path("tests/test_gbif_country_geographic_influence_20261020.py"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/geographic_dependence_posthoc_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_record_turnover_summary_20261020.json"),
    Path("provenance/reviewer_defenses/results/gbif_country_20261020/GBIF_COUNTRY_RECORD_TURNOVER_READOUT_20261020.md"),
}

IDENTITY_REPLACEMENTS = {
    "zuizui0223": "anonymous-user",
    "ZHANG": "[BLINDED]",
    "Ruiqi": "[BLINDED]",
}

PROVENANCE_KEYS = {
    "workflow",
    "workflow_run_id",
    "run_id",
    "artifact_id",
    "artifact_name",
    "artifact_digest",
    "head_sha",
    "workflow_head_sha",
    "render_workflow_run_id",
    "descriptor_workflow_run_id",
}

TEXT_SUFFIXES = {
    ".py", ".md", ".json", ".toml", ".txt", ".yml", ".yaml", ".cff", ".tsv"
}

RAW_PRESERVE_PREFIXES = (
    "data/external/",
    "data/frozen/",
)

ANON_README = """# Butterfly specialization ecology — anonymous review snapshot

This archive is the anonymized data-and-code review supplement intended to be
uploaded directly with the manuscript during double-anonymous peer review. It is
a review snapshot of the analysis package for a global butterfly-specialization
manuscript.

It contains the scientific code, frozen protocols, de-identified result receipts, the v0.2 supplementary robustness tables,
the exact LepTraits input used by the reconstruction, and byte-exact source inputs
needed to regenerate the three main manuscript figures and the climate supplementary figure. Identifying repository metadata,
Git history, author metadata, cover letters, title pages, release metadata, and
internal hosting/run identifiers are intentionally omitted or redacted for review.

## Quick checks

Python 3.11+ is required; the reference CI uses Python 3.12.

    python -m pip install --upgrade pip
    python -m pip install -e ".[test,figure]"
    pytest -q

The two WCVP/HOSTS sidecar scripts use base R only. Full environment notes are in
docs/REPRODUCIBILITY.md.

## Regenerate the manuscript figures

    python scripts/render_butterfly_specialization_v02_figures.py \
      --anthropogenic-csv data/frozen/figure_sources/anthropogenic_species_metrics.csv \
      --matched-null-json provenance/reviewer_defenses/results/butterfly_resource_expansion_matched_null_v0.2.json \
            --hostbias-null-json provenance/reviewer_defenses/results/butterfly_resource_expansion_hostbias_null_v0.1.json \
      --plant-prominence-json provenance/reviewer_defenses/results/butterfly_host_plant_prominence_expansion_v0.1.json \
      --host-concentration-json provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json \
      --occurrence-csv data/frozen/figure_sources/occurrence_resource_validation_species.csv \
      --occurrence-null-json provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json \
      --occurrence-species-robustness-json provenance/reviewer_defenses/results/butterfly_occurrence_species_robustness_v0.1.json \
      --ceiling-json provenance/reviewer_defenses/results/butterfly_expansion_ceiling_sensitivity_v0.1.json \
      --regional-json provenance/reviewer_defenses/results/butterfly_regional_robustness_v0.1.json \
      --climate-csv data/frozen/figure_sources/climate_distance_sensitivity_species.csv \
      --climate-effect-json provenance/reviewer_defenses/results/butterfly_climate_effect_size_v0.1.json \
      --output-dir results/manuscript-figures

The review snapshot intentionally does not contain author names or a public
repository URL. A public, citable archive is supplied only in the final public
version after double-anonymous review.
"""


def sanitize_text(text: str) -> str:
    for old, new in IDENTITY_REPLACEMENTS.items():
        text = text.replace(old, new)
    text = text.replace("https://github.com/anonymous-user/chocho", "[ANONYMIZED_REPOSITORY]")
    text = text.replace("https://github.com/anonymous-user/TTF", "[ANONYMIZED_PRECURSOR_REPOSITORY]")
    return text


def redact_json(value):
    if isinstance(value, dict):
        out = {}
        for key, item in value.items():
            if key in PROVENANCE_KEYS:
                continue
            out[key] = redact_json(item)
        return out
    if isinstance(value, list):
        return [redact_json(item) for item in value]
    if isinstance(value, str):
        return sanitize_text(value)
    return value


def should_preserve_raw(rel: Path) -> bool:
    s = rel.as_posix()
    return any(s.startswith(prefix) for prefix in RAW_PRESERVE_PREFIXES)


def copy_one(source: Path, destination: Path, rel: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if should_preserve_raw(rel):
        shutil.copyfile(source, destination)
        return

    suffix = source.suffix.lower()
    if suffix == ".json":
        try:
            payload = json.loads(source.read_text(encoding="utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            shutil.copyfile(source, destination)
            return
        destination.write_text(
            json.dumps(redact_json(payload), indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        return

    if suffix in TEXT_SUFFIXES or source.name == "pyproject.toml":
        try:
            text = source.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            shutil.copyfile(source, destination)
            return
        destination.write_text(sanitize_text(text), encoding="utf-8")
        return

    shutil.copyfile(source, destination)


def stage_snapshot(stage: Path) -> None:
    for item in COPY_PATHS:
        source = ROOT / item
        if not source.exists():
            raise FileNotFoundError(source)
        if source.is_file():
            rel = source.relative_to(ROOT)
            if rel in EXCLUDE_RELATIVE:
                continue
            copy_one(source, stage / rel, rel)
            continue

        for source_file in sorted(p for p in source.rglob("*") if p.is_file()):
            rel = source_file.relative_to(ROOT)
            if rel in EXCLUDE_RELATIVE:
                continue
            if "__pycache__" in rel.parts or rel.suffix in {".pyc", ".pyo"}:
                continue
            copy_one(source_file, stage / rel, rel)

    (stage / "README.md").write_text(ANON_README, encoding="utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(stage: Path) -> None:
    rows = []
    for path in sorted(p for p in stage.rglob("*") if p.is_file()):
        rel = path.relative_to(stage).as_posix()
        if rel == "SHA256SUMS":
            continue
        rows.append(f"{sha256(path)}  {rel}")
    (stage / "SHA256SUMS").write_text("\n".join(rows) + "\n", encoding="utf-8")


def scan_identity(stage: Path) -> list[str]:
    failures: list[str] = []
    banned = ("zuizui0223", "ruiqi", "zhang.ruiqi")
    for path in sorted(p for p in stage.rglob("*") if p.is_file()):
        if should_preserve_raw(path.relative_to(stage)):
            continue
        try:
            text = path.read_text(encoding="utf-8").lower()
        except UnicodeDecodeError:
            continue
        for token in banned:
            if token in text:
                failures.append(f"{path.relative_to(stage)} contains {token!r}")
    return failures


def deterministic_zip(stage: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in sorted(p for p in stage.rglob("*") if p.is_file()):
            rel = path.relative_to(stage).as_posix()
            info = zipfile.ZipInfo(rel, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            zf.writestr(info, path.read_bytes())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--keep-stage", type=Path, default=None)
    args = ap.parse_args()

    tmp = None
    if args.keep_stage:
        stage = args.keep_stage
        if stage.exists():
            shutil.rmtree(stage)
        stage.mkdir(parents=True)
    else:
        tmp = tempfile.TemporaryDirectory(prefix="chocho-anonymous-review-")
        stage = Path(tmp.name)

    try:
        stage_snapshot(stage)
        failures = scan_identity(stage)
        if failures:
            raise RuntimeError("anonymous snapshot identity scan failed:\n" + "\n".join(failures))
        write_manifest(stage)
        deterministic_zip(stage, args.output)
        print(json.dumps({
            "output": str(args.output),
            "files": sum(1 for p in stage.rglob("*") if p.is_file()),
            "sha256": sha256(args.output),
        }, indent=2))
    finally:
        if tmp is not None:
            tmp.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

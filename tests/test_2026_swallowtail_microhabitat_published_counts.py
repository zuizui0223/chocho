"""Arithmetic consistency tests of printed 2026 article tables.

Tests do not claim original ramet survey access, biological effects or model fit.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/audit_2026_swallowtail_microhabitat_published_counts.py"
spec=importlib.util.spec_from_file_location("published_counts",SOURCE)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_four_site_and_month_margins_point_to_215():
    z=mod.audit()
    assert z["published_methods_total_quadrats"]==255
    assert z["published_site_total_quadrats"]==215
    assert z["published_table_header_total_quadrats"]==255
    assert z["source_based_plausible_quadrat_denominator"]==215
    assert z["site_category_and_month_category_sums_agree"] is True
    assert z["site_totals_agree_with_published_table_rows"] is True
    assert z["article_quadrat_header_is_internally_consistent"] is False
    assert z["inferred_internal_115_not_155_S_only_count"]["discrepancy"]==40
    assert z["denominator_requires_author_confirmation"] is True

def test_separate_ramets_and_quadrats_without_pseudoreplication():
    z=mod.audit()
    assert z["original_ramet_total"]==875
    assert z["ramets_with_both_species"]==5
    assert z["quadrat_both_count_from_all_source_breakdowns"]==11
    assert z["ramet_and_quadrat_statistics_not_independent_samples"] is True
    assert z["ramet_overlap_is_not_causal_competition"] is True

def test_never_upgrade_published_tables_to_authentic_raw_data():
    z=mod.audit()
    assert z["source_kind"]=="PUBLICATION_PRINTED_TABLES_ONLY_NOT_ORIGINAL_RECORDS"
    assert z["original_raw_field_data_accessible_in_this_audit"] is False
    assert z["new_ecological_effect_estimated"] is False
    assert z["geographic_transport_to_kyoto_A_debilis_not_identified"] is True
    assert z["GEB_PR38_untouched"] is True

def test_check_detects_corrected_paper_scenario_without_inventing_it():
    """Hypothetical test edit; only source-verified 2026 counts belong in report."""
    original=mod.PUBLICATION
    edited={**original,"declared_total_quadrats_in_methods":215,
            "table1_quadrat_labels":{**original["table1_quadrat_labels"],"S_only":115}}
    z=mod.audit(edited)
    assert z["article_quadrat_header_is_internally_consistent"] is True
    assert z["published_site_total_quadrats"]==215

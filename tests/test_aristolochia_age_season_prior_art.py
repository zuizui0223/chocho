"""Synthetic publication-only parser fixtures, not original experimental rows."""
from __future__ import annotations
import importlib.util
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/audit_aristolochia_age_season_prior_art.py"
spec=importlib.util.spec_from_file_location("agepaper",P)
assert spec and spec.loader
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)

TEXT=("10.3389/fpls.2023.1145363 Sericinus montela Spodoptera exigua "
      "Aristolochia contorta 1st-year July September C/N "
      "aristolochic acid 1 aristolochic acid 2 "
      "six days 24°C leaf area 50% leaves replaced "
      "83.80 7.61 2.26 8.952 270.727")

def test_fulltext_xml_without_individual_data_inference():
    xml=("<article><front><article-meta><article-title>"
         +TEXT+"</article-title></article-meta></front></article>").encode()
    x=s.audit_text(xml,"https://example.org/article.xml")
    assert x["status"]=="AUTHOR_SOURCE_ARTICLE_VERIFIED"
    assert x["has_original_rowlevel_per_plant_data"] is False
    assert x["source_kind"].endswith("NOT_original_individual_records")
    assert all(x["matched_checks"].values())

def test_html_parser_and_incomplete_source_fail_closed():
    html=("<html><body><main>"+TEXT+"</main></body></html>").encode()
    x=s.audit_text(html,"https://example.org/article.html")
    assert x["status"]=="AUTHOR_SOURCE_ARTICLE_VERIFIED"
    bad=s.audit_text(b"<html><body>some generic abstract, no quantitative paper</body></html>","https://example.org/incomplete.html")
    assert bad["status"]=="RETRIEVED_SOURCE_NOT_FULLY_MATCHED"
    assert bad["matched_checks"]["correct_DOI"] is False
    assert bad["matched_checks"]["reported_temperature_24_C"] is False

def test_study_claims_are_previous_2023_publication_not_new_butterfly_effect():
    assert s.DOI=="10.3389/fpls.2023.1145363"
    assert len(s.URLS)==2
    assert any("europepmc" in url for url in s.URLS)

"""Source-independent tests for exact 5-case introduced-host cohort selection."""
from __future__ import annotations
import importlib.util
from pathlib import Path
spec=importlib.util.spec_from_file_location("source_host",Path(__file__).resolve().parents[1]/"scripts/derive_gbif_positive48_introduced_hosts_20261021.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_botanical_accepted_name_classification_is_not_fuzzy():
    assert m.name_status(["Asclepias curassavica"])=="EXACT_SOURCE_BINOMIAL_GBIF_RESOLUTION_REQUIRED"
    assert m.name_status(["Citrus reticulata Blanco"])=="NOT_EXACT_WCVP_BINOMIAL"
    assert m.name_status(["Citrus x aurantium"])=="NOT_EXACT_WCVP_BINOMIAL"
    assert m.name_status(["Asclepias curassavica","Asclepias physocarpa"])=="AMBIGUOUS_OR_UNRESOLVED_ACCEPTED_NAME"
    assert m.name_status([])=="AMBIGUOUS_OR_UNRESOLVED_ACCEPTED_NAME"

def test_source_five_case_cohort_is_fixed():
    assert len(m.FIXED_CASES)==5
    assert {v[0] for v in m.FIXED_CASES}=={"Anteos clorinde","Colobura dirce","Papilio demoleus"}
    assert m.FIXED_CASES=={
      ("Anteos clorinde","CU"),("Colobura dirce","VE"),
      ("Papilio demoleus","CU"),("Papilio demoleus","DO"),("Papilio demoleus","HT")}

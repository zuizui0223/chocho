"""Basic no-network checks for safe 48-unit GBIF source replay."""
import importlib.util
from pathlib import Path
import pytest
spec=importlib.util.spec_from_file_location("positive48",Path(__file__).resolve().parents[1]/"scripts/summarize_gbif_positive48_subcountry_20261020.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def test_category_sum_and_absent_life_stage():
    rows=[{"record_count":3,"facets":{}},
          {"record_count":3,"facets":{"life_stage":{"categories":{"ADULT":2},
             "assigned_records":2,"unassigned_records":1}}}]
    assert m.category_totals(rows,"life_stage")=={"ADULT":2}
    assert m.category_totals(rows,"establishment_means")=={}
    rows[1]["facets"]["life_stage"]["assigned_records"]=3
    with pytest.raises(ValueError,match="facet totals"):
        m.category_totals(rows,"life_stage")
def test_strict_completed_and_source_gates():
    with pytest.raises(ValueError):
        m.summarise({"status":"INCOMPLETE"})
    with pytest.raises(ValueError,match="hashes drifted"):
        m.summarise({"status":"COMPLETE_48","pilot":False,"source":{}})

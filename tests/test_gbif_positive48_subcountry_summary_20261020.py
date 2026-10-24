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

def test_archived_48_positive_gadm_stage_source_replays():
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    raw=json.loads((root/"provenance/reviewer_defenses/results/gbif_country_20261020/gbif_positive48_gadm_lifestage_facets_20261020.json").read_text(encoding="utf8"))
    saved=json.loads((root/"provenance/reviewer_defenses/results/gbif_country_20261020/gbif_positive48_subcountry_stage_summary_20261020.json").read_text(encoding="utf8"))
    z=m.summarise(raw)
    assert z["positive_species_country_cells"]==48
    assert z["indexed_records"]==saved["indexed_records"]==756
    assert z["gadm1_assigned_records"]==saved["gadm1_assigned_records"]==736
    assert z["life_stage_raw_labels"]==saved["life_stage_raw_labels"]
    assert z["establishment_means_annotated_records"]==0
    assert saved["juvenile_label_count"]==39
    assert saved["juvenile_label_cells"]==11
    assert saved["juvenile_label_species"]=={
       "Papilio demoleus":27,"Colobura dirce":11,"Anteos clorinde":1}
    assert saved["host_use_or_establishment_verified"] is False

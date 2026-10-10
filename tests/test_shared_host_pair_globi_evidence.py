from __future__ import annotations
import importlib.util
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/audit_shared_host_pair_globi_evidence.py"
spec=importlib.util.spec_from_file_location("host_pair",SCRIPT)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_strict_stage_and_source_gates():
    base={"source_taxon_name":"Junonia coenia",
          "target_taxon_name":"Plantago lanceolata",
          "source_specimen_life_stage":"larva",
          "event_date":"2023-07-01","latitude":"26.1","longitude":"-80.2",
          "study_citation":"independent source","study_external_id":"x",
          "source_specimen_occurrence_id":"https://example.org/original"}
    counts, ids=mod.count_source_rows([base],"Junonia coenia")
    assert counts["independent_provenance"]==1
    assert len(ids)==1
    bad=[{**base,"source_specimen_life_stage":"adult"},{**base,"target_taxon_name":"Plantago major"},
         {**base,"study_citation":"globalbioticinteractions/hosts","study_external_id":""}]
    counts,_=mod.count_source_rows(bad,"Junonia coenia")
    assert counts["independent_provenance"]==0

def test_undated_or_unmapped_not_verified():
    base={"latitude":"", "longitude":"", "event_date":"2020"}
    assert not mod.valid_geo_and_year(base)
    assert not mod.valid_geo_and_year({"latitude":"0","longitude":"0","event_date":"2020"})
    assert mod.valid_geo_and_year({"latitude":"26.0","longitude":"-80.0","event_date":"2020"})

def test_only_species_without_substitution():
    assert mod.SPECIES==["Junonia coenia","Anartia jatrophae"]
    assert not mod.exact_species("Junonia zonalis","Junonia coenia")

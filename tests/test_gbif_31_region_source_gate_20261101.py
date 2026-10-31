"""Test GBIF archival schema logic without external network or raw GPS data."""
import importlib.util
import io
from pathlib import Path
import zipfile

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/audit_gbif_31_butterfly_wgsrpd3_source.py"


def load():
    pytest.importorskip("shapely")
    spec=importlib.util.spec_from_file_location("gbif31_source_audit",SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_canonical_species_from_source_with_lowercase_fields():
    m=load()
    assert m.normalize_species({"species":"Vanessa cardui"})=="vanessa cardui"
    assert m.normalize_species({"genus":"Danaus","specificepithet":"plexippus"})=="danaus plexippus"
    assert m.normalize_species({"scientificname":"Aglais urticae (Linnaeus, 1758)"})=="aglais urticae"
    assert m.normalize_species({"species":""}) is None
    assert m.EXPECTED_RETURNED==53144
    assert m.SOURCE_DOI=="10.15468/dl.pp5nc9"


def test_archive_table_selection_does_not_confuse_citation_with_records(tmp_path):
    m=load()
    path=tmp_path/"toy.zip"
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("citation.txt","Citation to GBIF download")
        header="gbifID\tdecimalLatitude\tdecimalLongitude\tspecies\tdatasetKey\n"
        z.writestr("occurrence.txt",header+("1\t1.0\t2.0\tVanessa cardui\ta\n"*40))
    with zipfile.ZipFile(path) as z:
        (name,delim,columns),diagnostics=m.find_occurrence_member(z)
    assert name=="occurrence.txt"
    assert delim=="\t"
    assert "decimalLatitude" in columns
    assert len(diagnostics)==1

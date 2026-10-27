"""Offline source contract checks for public Dryad insect arrival and establishment."""
import io
import sys
import zipfile
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_original_insect_arrival_dryad_schema_20261027 import (
    DOIS, EXPECTED, schema_rows,sha,acquire,SCHEMA
)


def test_2023_transport_and_establishment_are_not_same_species_response():
    assert set(DOIS)=={"insect_flows_2023","lepidoptera_establishment_2022"}
    assert {"transport_flows_updated.csv","establishment_flows.csv"} <=EXPECTED["insect_flows_2023"]
    assert "S5_list_of_non-native_species.txt" in EXPECTED["lepidoptera_establishment_2022"]
    assert "S9_plant-commodities.txt" in EXPECTED["lepidoptera_establishment_2022"]


def test_csv_and_tab_delimited_source_grains_are_detected_separately():
    data="family\tHS06\tHS12\nLycaenidae\t1\t4\n".encode()
    report=schema_rows("S9_plant-commodities.txt",data)
    assert report["source_header_column_count"]==3
    assert report["source_nonempty_lines"]==2
    assert len(report["first_5_raw_data_line_previews"])==1
    csvdata=b"source,years,destination,class,richness\nNZL,1990s,JPN,Plant products,12\n"
    r=schema_rows("transport_flows_updated.csv",csvdata)
    assert r["first_line_columns"]==["source","years","destination","class","richness"]


def test_dryad_dataset_archive_member_integrity_and_no_taxon_fabrication(monkeypatch):
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,"w") as z:
        z.writestr("source/transport_flows_updated.csv","source,destination,richness\nUSA,JPN,2\n")
        z.writestr("source/establishment_flows.csv","native_region,destination,richness\nNA,JPN,1\n")
        z.writestr("source/README.md","Source docs")
    def mock_request(url):
        return buf.getvalue(),url,"application/zip"
    import audit_original_insect_arrival_dryad_schema_20261027 as src
    monkeypatch.setattr(src,"get_bytes",mock_request)
    result=acquire("doi:10.5061/dryad.dv41ns243")
    assert result["dataset_doi"]=="doi:10.5061/dryad.dv41ns243"
    assert result["file_names"]==["README.md","establishment_flows.csv","transport_flows_updated.csv"]
    assert len(result["public_source_archive_SHA256"])==64


def test_dryad_zip_rejects_same_basename_and_path_traversal(monkeypatch):
    import audit_original_insect_arrival_dryad_schema_20261027 as src
    def zipdata(items):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            for name in items:
                z.writestr(name,"test\n")
        return buf.getvalue()
    monkeypatch.setattr(src,"get_bytes",lambda url:(zipdata(["../README.md"]),url,"application/zip"))
    with pytest.raises(ValueError,match="path traversal"):
        acquire("doi:test")
    monkeypatch.setattr(src,"get_bytes",lambda url:(
        zipdata(["a/S5_list_of_non-native_species.txt","b/S5_list_of_non-native_species.txt"]),
        url,"application/zip"
    ))
    with pytest.raises(ValueError,match="Ambiguous"):
        acquire("doi:test")

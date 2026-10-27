"""No-network Zenodo public establishment roster provenance, no false absences."""
import csv
import hashlib
import json
import sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_zenodo_2024_nonnative_insect_establishment_source_20261027 as z


def test_2024_open_original_occurrence_not_2021_old_version():
    assert z.RECORD_ID==11519262
    assert z.API.endswith("/11519262")
    assert "occurrence_table.csv" in z.REQUIRED
    assert "attribute_table.csv" in z.REQUIRED


def test_real_table_header_and_source_rows_are_preserved_without_binary_zeros():
    original=b"taxon,region,date\nPieris rapae,NZ,1930\nThymelicus lineola,CA,\n"
    parsed=z.parse_csv("occurrence_table.csv",original)
    assert parsed["original_header_columns"]==["taxon","region","date"]
    assert parsed["original_source_data_rows"]==2
    assert parsed["source_unlabeled_or_empty_value_counts"]=={"date":1}
    assert parsed["original_source_sha256"]==hashlib.sha256(original).hexdigest()


def test_exact_public_metadata_and_5_tables(monkeypatch):
    names=sorted(z.REQUIRED)
    values={name:b"species,region\nPieris rapae,NZ\n" for name in names}
    meta={
        "id":11519262,
        "metadata":{
            "title":"International Non-native Insect Establishment Data",
            "version":"0.2.1",
            "doi":"10.5281/zenodo.11519262",
            "publication_date":"2024-06-07"
        },
        "files":[{"key":name,"size":len(values[name]),
                  "checksum":"md5:"+hashlib.md5(values[name]).hexdigest(),
                  "links":{"self":"https://zenodo.test/"+name}} for name in names]
    }
    def get(url,max_bytes=9_000_000):
        if url==z.API:
            return json.dumps(meta).encode(),url,"application/json"
        filename=url.split("/")[-1]
        return values[filename],url,"text/csv"
    monkeypatch.setattr(z,"get",get)
    x=z.fetch()
    assert x["original_version"]=="0.2.1"
    assert len(x["required_source_files"])==5
    assert x["required_source_files"]["occurrence_table.csv"]["original_source_data_rows"]==1
    assert any("NOT proven" in s for s in x["scientific_inference_stop"])


def test_corrupt_source_member_md5_fails_closed(monkeypatch):
    names=sorted(z.REQUIRED)
    data=b"species,region\nPieris rapae,NZ\n"
    meta={"id":z.RECORD_ID,"files":[{"key":name,"size":len(data),
        "checksum":"md5:"+"0"*32,"links":{"self":"https://test/"+name}}
        for name in names]}
    monkeypatch.setattr(z,"get",lambda url,max_bytes=9_000_000:(
        (json.dumps(meta).encode() if url==z.API else data),url,"text/plain"))
    with pytest.raises(ValueError,match="md5 mismatch"):
        z.fetch()


def test_source_missing_species_table_is_not_allowed(monkeypatch):
    meta={"id":z.RECORD_ID,"files":[]}
    monkeypatch.setattr(z,"get",lambda url,max_bytes=9_000_000:(
        json.dumps(meta).encode(),url,"application/json"))
    with pytest.raises(ValueError,match="metadata missing"):
        z.fetch()

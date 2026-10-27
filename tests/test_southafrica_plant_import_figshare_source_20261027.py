"""No-network original Figshare source collection and sample-status QA."""
from __future__ import annotations
import json,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_southafrica_plant_import_figshare_source_20261027 as f


def test_true_original_2021_figshare_collection_not_wrong_2022_article():
    assert f.PUBLIC_COLLECTION_ID==5180681
    assert f.PUBLIC_API=="https://api.figshare.com/v2/collections/5180681/articles"
    assert f.EXPECTED_FILES==("List of contaminants","Metadata of contaminants")


def test_semicolon_sample_vocabulary_preserves_source_NA_and_zero():
    data=("catalogNumber;sample.source;commodity.type;Lepidoptera;Hemiptera\n"
          "1;audit;Seeds;NA;0\n2;intervention;Plants;Pieris rapae;NA\n").encode()
    r=f.field_values("List of contaminants.csv",data)
    assert r["original_delimiter"]==";"
    assert r["source_row_count"]==2
    assert r["column_count"]==5
    x=r["selected_original_raw_column_vocabulary"]["Lepidoptera"]
    assert sorted((y["value"],y["n"]) for y in x["top8_values"])==[
       ("NA",1),("Pieris rapae",1)
    ]


def fixture():
    c1="catalogNumber;Lepidoptera\n1;NA\n".encode()
    c2="Taxon;Order;status\nPieris rapae;Lepidoptera;absent\n".encode()
    collection=[{"id":10001,"title":"Original 2021 samples"}]
    detail={"id":10001,"title":"Original 2021 source tables","files":[
        {"id":111,"name":"List of contaminants on SA plant imports 1994-2019.csv",
         "size":len(c1),"download_url":"https://figshare.test/111"},
        {"id":112,"name":"Metadata of contaminants on SA plant imports 1994-2019.csv",
         "size":len(c2),"download_url":"https://figshare.test/112"},
    ]}
    return c1,c2,collection,detail


def test_true_collection_article_then_two_original_raw_csvs(monkeypatch):
    c1,c2,collection,detail=fixture()
    def mock_get(url):
        if url==f.PUBLIC_API:return json.dumps(collection).encode(),url,"application/json"
        if url.endswith("/10001"):return json.dumps(detail).encode(),url,"application/json"
        if url.endswith("/111"):return c1,url,"text/csv"
        if url.endswith("/112"):return c2,url,"text/csv"
        raise AssertionError("Unexpected external URL "+url)
    monkeypatch.setattr(f,"source_get",mock_get)
    result=f.run()
    assert result["figshare_collection_ID"]==5180681
    assert result["original_public_dataset_DOI"]=="10.6084/m9.figshare.c.5180681"
    assert len(result["original_files"])==2
    assert [x["source_row_count"] for x in result["original_files"]]==[1,1]
    assert {x["original_source_article_ID"] for x in result["original_files"]}=={10001}


def test_original_data_bytes_must_match_verified_figshare_file_size(monkeypatch):
    c1,c2,collection,detail=fixture()
    detail["files"][0]["size"]=999
    def mock_get(url):
        if url==f.PUBLIC_API:return json.dumps(collection).encode(),url,"application/json"
        if url.endswith("/10001"):return json.dumps(detail).encode(),url,"application/json"
        if url.endswith("/111"):return c1,url,"text/csv"
        if url.endswith("/112"):return c2,url,"text/csv"
        raise AssertionError("Unexpected")
    monkeypatch.setattr(f,"source_get",mock_get)
    with pytest.raises(ValueError,match="size mismatch"):
        f.run()


def test_wrong_2022_reanalysis_article_not_accepted_as_2021_source(monkeypatch):
    collection=[{"id":14393723}]
    detail={"id":14393723,"title":"2022 follow-up","files":[
      {"id":77,"name":"data - models Saccaggi et al 2021 Ecol Appl.csv",
       "size":100,"download_url":"https://figshare.test/77"}]}
    def mock_get(url):
        if url==f.PUBLIC_API:return json.dumps(collection).encode(),url,"application/json"
        if url.endswith("/14393723"):return json.dumps(detail).encode(),url,"application/json"
        raise AssertionError("Must reject wrong original raw CSV")
    monkeypatch.setattr(f,"source_get",mock_get)
    with pytest.raises(ValueError,match="expects two CSV"):
        f.run()

"""No-network original Figshare source API/schema checks before scientific commodity rates."""
from __future__ import annotations

import io,json,sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_southafrica_plant_import_figshare_source_20261027 as f


def test_figshare_DOI_is_exact_original_anonymized_dataset():
    assert f.PUBLIC_ARTICLE_ID==14393723
    assert f.PUBLIC_API=="https://api.figshare.com/v2/articles/14393723"
    assert f.EXPECTED_FILES==("List of contaminants","Metadata of contaminants")


def test_semicolon_original_sample_structure_and_missing_test_semantics():
    data=("catalogNumber;sample.source;commodity.type;Lepidoptera;Hemiptera\n"
          "1;audit;Seeds;NA;0\n2;intervention;Plants;Pieris rapae;NA\n").encode()
    result=f.field_values("List of contaminants on SA plant imports 1994-2019.csv",data)
    assert result["original_delimiter"]==";"
    assert result["source_row_count"]==2
    assert result["column_count"]==5
    x=result["selected_original_raw_column_vocabulary"]["Lepidoptera"]
    assert x["distinct_nonempty_values"]==2
    assert sorted((y["value"],y["n"]) for y in x["top8_values"])==[
        ("NA",1),("Pieris rapae",1)
    ]


def test_source_original_metadata_exact_file_names_and_sizes(monkeypatch):
    c1=b'catalogNumber;Lepidoptera\\n1;NA\\n'.replace(b'\\n',b'\n')
    c2=b'Taxon;Order;status\\nPieris rapae;Lepidoptera;absent\\n'.replace(b'\\n',b'\n')
    meta={
       "id":14393723,"title":"Contaminant organisms recorded on plant products",
       "files":[
           {"id":111,"name":"List of contaminants on SA plant imports 1994-2019.csv",
            "size":len(c1),"download_url":"https://figshare.test/111"},
           {"id":112,"name":"Metadata of contaminants on SA plant imports 1994-2019.csv",
            "size":len(c2),"download_url":"https://figshare.test/112"},
       ]
    }
    def mock_get(url):
        if url==f.PUBLIC_API:return json.dumps(meta).encode(),url,"application/json"
        if url.endswith("/111"):return c1,url,"text/csv"
        if url.endswith("/112"):return c2,url,"text/csv"
        raise AssertionError("Unexpected external data fetch")
    monkeypatch.setattr(f,"source_get",mock_get)
    result=f.run()
    assert result["figshare_article_ID"]==14393723
    assert len(result["original_files"])==2
    assert [x["source_row_count"] for x in result["original_files"]]==[1,1]


def test_mismatched_origin_file_size_never_silently_accepted(monkeypatch):
    meta={
        "id":14393723,"title":"Test",
        "files":[
            {"id":11,"name":"List of contaminants.csv",
             "size":999,"download_url":"https://figshare.test/11"},
            {"id":12,"name":"Metadata of contaminants.csv",
             "size":5,"download_url":"https://figshare.test/12"},
        ]
    }
    monkeypatch.setattr(f,"source_get",lambda url:(
        json.dumps(meta).encode() if url==f.PUBLIC_API else b"header;x\n",
        url,"application/json"
    ))
    with pytest.raises(ValueError,match="size mismatch"):
        f.run()

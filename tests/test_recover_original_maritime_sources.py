import pytest,io,zipfile,hashlib
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from recover_original_GLSN2018_Excel_source import (
    FILES,github_blob_sha,sheet_tables
)
from recover_original_UNCTAD_PLSCIQ_legacy_archive import (
    WPI,original_wpi,HISTORIC_SHA,country_iso2,clean,parse_value,Q_START,Q_STOP
)

def test_original_2018_github_file_blobs_pinned():
    assert FILES["original_GLNS_ports_2018.xlsx"]["sha"]=="5486d8c204a7c49b375a14762fe5cb36c7877d96"
    assert FILES["original_GLNS_edges_2018.xlsx"]["sha"]=="18098cab432ddbcd7e0a7091fc1c4cedbda40fb7"
    assert github_blob_sha(b"hello")==hashlib.sha1(b"blob 5\x00hello").hexdigest()

def test_stdlib_xlsx_parser_handles_original_excel_shared_strings():
    workbook='<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="TestPorts" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rels='<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>'
    shared='<?xml version="1.0"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><si><t>Port</t></si><si><t>COL</t></si></sst>'
    sheet='<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="inlineStr"><is><t>Latitude</t></is></c></row><row r="2"><c r="A2" t="s"><v>1</v></c><c r="B2"><v>6.1</v></c></row></sheetData></worksheet>'
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,"w") as archive:
        for path,raw in [("xl/workbook.xml",workbook),("xl/_rels/workbook.xml.rels",rels),
                         ("xl/sharedStrings.xml",shared),("xl/worksheets/sheet1.xml",sheet)]:
            archive.writestr(path,raw)
    items=sheet_tables(buf.getvalue())
    assert len(items)==1 and items[0]["rows"]==[["Port","Latitude"],["COL","6.1"]]

def test_original_WPI_source_and_legacy_PLSCIQ_vintage():
    if not WPI.is_file():
        pytest.skip("Original WPI coordinates intentionally omitted from anonymous package")
    a=original_wpi()
    assert len(a)==194
    assert len({x["island_id"] for x in a})==51
    assert sum(x["included"].lower()=="true" for x in a)==182
    assert len(HISTORIC_SHA)==3 and Q_START=="2006-Q1" and Q_STOP=="2021-Q3"
    assert country_iso2("Japan")=="JP"
    assert clean("Santa Cruz de Tenerife")=="SANTACRUZDETENERIFE"
    assert parse_value("NA") is None
    assert parse_value("12.5")==12.5

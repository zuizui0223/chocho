"""Public IMMP original-data eligibility audit: schema is not cohort survival."""
from __future__ import annotations

import importlib.util
import zipfile
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"scripts/audit_external_immp_lifestage_schema.py"
spec=importlib.util.spec_from_file_location("public_immp_schema",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def synthetic_xlsx(path:Path) -> None:
    # The XML is a 1-sheet OOXML fixture with column names and a dummy
    # next-row observation; no biological data should be read or exported.
    ns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    pr="http://schemas.openxmlformats.org/package/2006/relationships"
    strings=["Site","Visit date","Milkweed species","Egg count","Larval instar","Adult count"]
    with zipfile.ZipFile(path,"w") as out:
        out.writestr("xl/workbook.xml",
           f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>'
           '<sheet name="Act2" sheetId="1" r:id="rId1"/></sheets></workbook>')
        out.writestr("xl/_rels/workbook.xml.rels",
           f'<Relationships xmlns="{pr}"><Relationship '
           'Id="rId1" Type="sheet" Target="worksheets/sheet1.xml"/>'
           '</Relationships>')
        out.writestr("xl/sharedStrings.xml",
           f'<sst xmlns="{ns}">'+''.join(
               f'<si><t>{v}</t></si>' for v in strings)+'</sst>')
        columns=[f'<c r="{chr(65+i)}1" t="s"><v>{i}</v></c>'
                 for i in range(len(strings))]
        out.writestr("xl/worksheets/sheet1.xml",
            f'<worksheet xmlns="{ns}"><sheetData><row r="1">'
            +''.join(columns)+'</row><row r="2">'
            '<c r="A2"><v>1234</v></c></row></sheetData></worksheet>')


def test_stage_headers_are_a_data_eligibility_check_only(tmp_path):
    p=tmp_path/"sample.xlsx"
    synthetic_xlsx(p)
    sheets=m.sheet_summaries(p)
    assert [x["worksheet"] for x in sheets]==["Act2"]
    assert sheets[0]["nonempty_row_count_upper_bound"]==2
    d=m.eligibility(sheets)
    assert d["stage_count_columns_exist_in_same_worksheet"] is True
    assert d["paired_native_versus_introduced_species_at_same_site"]=="NOT_IDENTIFIABLE_FROM_COLUMN_NAMES"
    assert d["cohort_linked_egg_to_recruit_survival"]=="NOT_IDENTIFIABLE_FROM_COLUMN_NAMES"
    assert d["test_of_negative_butterfly_recruitment"]=="NOT_PERFORMED"


def test_missing_stages_or_host_precludes_pretending_a_true_survival_test():
    d=m.eligibility([{
        "worksheet":"Act1",
        "candidate_variable_classes":{
            "egg":[],"juvenile":[],"adult":["adult count"],
            "milkweed":["milkweed"],"site":["site"],"survey_date":["date"],"effort":[]
        }
    }])
    assert d["stage_count_columns_exist_in_same_worksheet"] is False
    assert d["cohort_linked_egg_to_recruit_survival"]=="NOT_IDENTIFIABLE_FROM_COLUMN_NAMES"


def test_source_version_and_license_are_locked():
    assert m.SRC_DOI=="10.5061/dryad.mkkwh718x"
    assert m.SRC_VERSION==6
    assert m.EXPECTED_SHA256=="d9ed67778e96643ea4414da8971340b9c769a0a0aaaf70871700ecce36016a3d"


def test_invalid_workbook_sheet_relationship_fails(tmp_path):
    p=tmp_path/"bad.xlsx"
    synthetic_xlsx(p)
    with zipfile.ZipFile(p) as z:
        paths={x:z.read(x) for x in z.namelist()}
    paths["xl/_rels/workbook.xml.rels"]=paths["xl/_rels/workbook.xml.rels"].replace(
        b"worksheets/sheet1.xml",b"worksheets/not_there.xml"
    )
    with zipfile.ZipFile(p,"w") as z:
        for path,data in paths.items():
            z.writestr(path,data)
    with pytest.raises(ValueError,match="relationships"):
        m.sheet_summaries(p)

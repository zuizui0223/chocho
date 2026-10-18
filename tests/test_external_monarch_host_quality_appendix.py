"""Conservative quality-host source schema audit: not a fitness estimate."""
from __future__ import annotations
import importlib.util
import zipfile
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
# The official double-anonymous submission snapshot intentionally excludes
# exploratory independently sourced host-quality audit scripts.
if not (ROOT/"scripts/audit_external_monarch_host_quality_appendix.py").is_file():
    pytest.skip("Exploratory host-quality audit excluded from frozen original GEB snapshot", allow_module_level=True)
import sys
sys.path.insert(0,str(ROOT/"scripts"))
P=ROOT/"scripts/audit_external_monarch_host_quality_appendix.py"
spec=importlib.util.spec_from_file_location("quality_s1",P)
assert spec and spec.loader
x=importlib.util.module_from_spec(spec)
spec.loader.exec_module(x)


def make_docx(path):
    ns="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    table='<w:tbl>'
    for plant,cat in [("Asclepias curassavica","H3"),("Asclepias tuberosa","H2"),
                      ("Araujia sericifera","L3"),("Nerium oleander","N")]:
        table+=(
            '<w:tr><w:tc><w:p><w:r><w:t>'+plant+'</w:t></w:r></w:p></w:tc>'
            '<w:tc><w:p><w:r><w:t>'+cat+'</w:t></w:r></w:p></w:tc></w:tr>'
        )
    table+='</w:tbl>'
    body=f'<w:document xmlns:w="{ns}"><w:body>{table}</w:body></w:document>'
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("word/document.xml",body)


def test_study_table_inspection_returns_no_unsupported_biology(tmp_path):
    p=tmp_path/"source.docx"
    make_docx(p)
    tabs=x.tables_in_docx(p)
    assert len(tabs)==1
    result=x.audit_tables(tabs)
    assert result["tables"][0]["rows"]==4
    assert result["tables"][0]["binomial_like_rows"]==4
    assert result["tables"][0]["single_status_code_rows"]==4
    assert result["eligibility"]["observed_larval_recruitment"]=="NOT_MEASURED"
    assert result["eligibility"]["independent_40_HOSTS_accepted_taxa_crosswalk"]=="NOT_PERFORMED"
    assert result["eligibility"]["exact_127_status_row_table_present"] is False


def test_code_ambiguity_rejected():
    assert x.code_in_row(["Host status H2 or L3"] ) is None
    assert x.code_in_row(["Asclepias curassavica","H3"])=="H3"
    assert x.code_in_row(["Asclepias curassavica","X9"]) is None


def test_missing_docx_word_document_is_rejected(tmp_path):
    p=tmp_path/"invalid.docx"
    with zipfile.ZipFile(p,"w") as z:z.writestr("README.txt","hello")
    with pytest.raises(ValueError,match="OOXML"):
        x.tables_in_docx(p)

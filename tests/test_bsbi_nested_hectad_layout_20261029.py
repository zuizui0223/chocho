"""Do not leak source hectad rows; check nested ZIP source integrity."""
import io,sys,zipfile
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_nested_hectad_layout_20261029 import inspect_nested,NESTED

def make_nested(path="folder/data.csv"):
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:
        z.writestr(path,"taxon,hectad,DC2,DC3\nA alpha,SP01,1,0\n")
    outer=io.BytesIO()
    with zipfile.ZipFile(outer,"w") as z:z.writestr(NESTED,b.getvalue())
    return outer.getvalue()

def test_nested_csv_early_date_period_and_privacy():
    src=inspect_nested(make_nested())
    assert src["original_nested_member_count"]==1
    assert src["nested_original_member_schema_sample"][0]["row_count_no_original_records_exported"]==1
    assert src["nested_original_member_schema_sample"][0]["plausibly_dated"]
    assert "SP01" not in str(src) and "A alpha" not in str(src)

def test_nested_archive_not_dated_plant_row_substitute():
    with pytest.raises(ValueError,match="unsafe"):
        inspect_nested(make_nested("../visits.csv"))

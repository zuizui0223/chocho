"""Source-blind schema/size/privacy regression tests for new monarch study."""
import io
import sys
from pathlib import Path
from zipfile import ZipFile
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_monarch_annotated_host_source_20261028 import inspect_zip


def test_schema_and_row_count(tmp_path):
    path = tmp_path / "public.zip"
    with ZipFile(path, "w") as z:
        z.writestr("folder/study.csv", "observation,host,cultivated,date\n1,A,yes,2020-05-01\n2,B,no,2021-12-01\n")
    result = inspect_zip(path)
    assert result["status"] == "PUBLIC_SHARED_ARCHIVE_TABULAR_SOURCE_IDENTIFIED"
    assert result["source_members"][0]["row_count"] == 2
    assert result["source_members"][0]["columns"] == ["observation", "host", "cultivated", "date"]
    assert "2020-05-01" not in str(result)
    assert result["total_tabular_rows_all_files"] == 2


def test_unknown_non_zip_fails_closed(tmp_path):
    path = tmp_path / "error.html"
    path.write_text("<html>Dropbox blocked</html>")
    with pytest.raises(Exception):
        inspect_zip(path)


def test_zip_traversal_fails_closed(tmp_path):
    path = tmp_path / "public.zip"
    with ZipFile(path, "w") as z:
        z.writestr("../bad.csv", "a,b\n1,2\n")
    with pytest.raises(ValueError, match="Unsafe ZIP member"):
        inspect_zip(path)

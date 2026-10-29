import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_chocho_author_UK_host_identity_20261029 import audit

def test_missing_author_source_no_fabrication(tmp_path):
    p=tmp_path/"SOURCE_MANIFEST.json";p.write_text("{}")
    res=audit(tmp_path)
    assert res["status"]=="NO_ACCEPTED_HOST_IDENTITY_TABLE_IN_ORIGINAL_AUTHORIZED_ARCHIVE"
    assert res["host_file_count"]==0

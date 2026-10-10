"""Audits source provenance only: no biological effects should be claimed."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_realized_host_window_source.py"
spec = importlib.util.spec_from_file_location("host_source_audit", SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_html_access_does_not_pass_as_original_data():
    class FakeResponse:
        status = 200
        headers = {"Content-Type": "text/html"}
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def read(self, _):
            return b"<html>login required</html>"

    with patch.object(module, "urlopen", return_value=FakeResponse()):
        result = module.fetch_source("phenology", "https://example.org/file")
    assert result["status"] == "SOURCE_ACCESS_BLOCKED"
    assert "sha256" not in result


def test_source_protocol_stops_without_outcome_information():
    protocol = json.loads((ROOT / "docs/exploratory/REALIZED_HOST_WINDOW_PROTOCOL_V01.json").read_text(encoding="utf-8"))
    assert protocol["independent_source"]["row_schema_verified"] is False
    assert protocol["raw_source_gate"]["fail_closed"] is True
    assert protocol["primary_model"]["biological_endpoint"].startswith("observed egg")
    assert "main" in protocol["submission_isolation"]


def test_provenance_record_remains_noncausal(tmp_path):
    with patch.object(module, "fetch_source", return_value={"status": "SOURCE_ACCESS_BLOCKED"}):
        result = module.audit(tmp_path / "receipt.json")
    assert result["source_bytes_complete"] is False
    assert result["row_level_ecological_effect_estimated"] is False
    assert result["main_manuscript_changed"] is False

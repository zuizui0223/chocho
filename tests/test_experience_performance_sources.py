"""No inference if downloadable archive is absent or not an .xlsx binary."""
from __future__ import annotations

import importlib.util
import io
import json
from pathlib import Path
from unittest.mock import Mock

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "scripts/audit_experience_performance_sources.py"
spec = importlib.util.spec_from_file_location("source_audit", path)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class FakeResponse:
    def __init__(self, raw):
        self.raw = raw
        self.status = 200
        self.headers = {"Content-Type":"application/octet-stream"}
    def __enter__(self): return self
    def __exit__(self,*args): return False
    def read(self,size): return self.raw


def test_error_page_cannot_be_called_xlsx():
    f = mod.audit_file("trial.xlsx", "https://example.org/file",
                       opener=lambda *a, **k: FakeResponse(b"<html>Access denied</html>"*600))
    assert f["status"] == "SOURCE_ACCESS_BLOCKED"
    assert "sha256" not in f


def test_binary_source_gate_checks_real_magic_and_bounds():
    f = mod.audit_file("trial.xlsx","https://example.org/file",
                       opener=lambda *a, **k: FakeResponse(b"PK\x03\x04"+b"1"*9000))
    assert f["status"] == "SOURCE_BYTES_VERIFIED_SCHEMA_UNCHECKED"
    assert f["byte_count"] == 9004


def test_hypothesis_does_not_call_native_learning_plant_globalization():
    protocol=json.loads((ROOT/"docs/exploratory/EXPERIENCE_PERFORMANCE_SOURCE_GATE_PROTOCOL_V01.json").read_text(encoding="utf-8"))
    assert "Both" in protocol["published_prior_art"]["botanical_origin"]
    assert protocol["source_access_gate"]["requires_both_files"] is True
    assert protocol["analysis_gate"]["stop"][-1].startswith("No introduced-host")

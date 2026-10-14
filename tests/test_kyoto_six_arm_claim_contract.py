"""Frozen scientific hypothesis mapping independent of synthetic test outcomes."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts"
sys.path.insert(0, str(SCRIPT))
P = SCRIPT / "analyze_kyoto_six_arm_identity_food.py"
SPEC = importlib.util.spec_from_file_location("six_arm_claim_contract", P)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)
DOC = ROOT / "docs/exploratory/KYOTO_SIX_ARM_PRIMARY_SECONDARY_CONTRACT_V02.json"


def test_contract_is_strictly_coded_before_biological_outcomes():
    data = json.loads(DOC.read_text(encoding="utf-8"))
    assert data["status"] == "PRE_FIELD_DESIGN_REVISION_NOT_PUBLIC_PREREGISTRATION"
    assert data["experiment_executed"] is False
    assert data["primary"]["primary_selected_before_any_new_field_outcomes"]
    assert not data["planning"]["pilot_block_variance_observed"]
    assert data["source_host_plant"]["native_in_Japan"] is True
    assert data["main_GEB_manuscript_unchanged"] is True
    assert data["randomization_unit"].startswith("WHOLE")


def test_primary_and_secondary_signed_coefficients_match_frozen_code():
    data = json.loads(DOC.read_text(encoding="utf-8"))
    assert mod.PRIMARY == data["primary"]["id"] == "S_competitor_x_food"
    assert set(mod.ARMS) == set(data["treatments"])
    assert all(
        tuple(original) == mod.ARMS[arm]
        for arm, original in data["treatments"].items()
    )
    assert mod.COEFFS[mod.PRIMARY] == data["primary"]["coefficients"]
    secondary = data["key_secondary"]
    assert secondary["id"] != mod.PRIMARY
    assert mod.COEFFS[secondary["id"]] == secondary["coefficients"]
    assert set(mod.COEFFS) == (
        {data["primary"]["id"], secondary["id"]} |
        set(data["other_secondary"])
    )
    assert data["original_focal_cohort_cannot_be_augmented_postassignment"]


def test_pupal_fate_accounting_cannot_claim_nonemergence_is_death():
    data = json.loads(DOC.read_text(encoding="utf-8"))
    statuses = data["outcomes"]
    assert set(statuses) == mod.FATE_VALUES
    assert statuses["diapause_alive"].endswith("not death")
    assert statuses["pupa_alive_pending_emergence"].endswith("not death")
    assert "diapause_alive" not in mod.FAILURE
    assert "pupa_alive_pending_emergence" in mod.UNRESOLVED

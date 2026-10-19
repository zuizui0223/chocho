"""Test statistical reliability gates used for post-hoc geographic-null MCMC."""
from __future__ import annotations

import importlib.util
import random
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
source=ROOT/"scripts/audit_butterfly_spatial_null_longchain.py"
spec=importlib.util.spec_from_file_location("spatial_audit",source)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_mixing_gate_accepts_independent_fixed_truth_toy_draws() -> None:
    rng=random.Random(617)
    independent=[rng.gauss(50,1) for _ in range(99)]
    check=m.mixing_check(independent)
    assert check["lag1"] <= 0.5
    assert check["approximate_ess"] >= 30
    assert check["pass"]
    assert check["n"]==99


def test_severely_autocorrelated_chain_fails() -> None:
    x=[0.0]
    rng=random.Random(771)
    for _ in range(98):
        x.append(0.95*x[-1]+rng.gauss(0,.2))
    result=m.mixing_check(x)
    assert not result["pass"]
    assert result["lag1"] > 0.5 or result["approximate_ess"] < 30


def test_trend_fails_quarter_drift_check_even_if_many_draws() -> None:
    rng=random.Random(554)
    x=[0.1*i+rng.gauss(0,0.05) for i in range(99)]
    result=m.mixing_check(x)
    assert result["late_minus_early"]>0.5*result["full_sd"]
    assert not result["drift_gate"]


def test_constant_or_too_short_worlds_do_not_claim_ess() -> None:
    assert m.integrated_acf_ess([3]*99)==(None,None,None)
    assert not m.mixing_check([3]*99)["pass"]
    assert m.integrated_acf_ess([2,1,2])==(None,None,None)


def test_failure_to_produce_two_evaluable_chains_abstains() -> None:
    failed=m.qualify({"status":"NOT_EVALUABLE","chains":[
        {"status":"NOT_EVALUABLE_SPATIAL_MIXING","seed":20261017}
    ]})
    assert failed["status"]=="NOT_EVALUABLE_MCMC_CALIBRATION"


def test_protocol_matches_code_and_blocks_parameter_fishing() -> None:
    assert m.CALIPER_KM==1500
    assert m.SEEDS==[20261017,20261018]
    assert m.DRAW_COUNT==99
    assert m.BURN_FACTOR==30
    assert m.SPACING_FACTOR==20
    text=(ROOT/"provenance/reviewer_defenses/protocols/SPATIAL_NULL_EXTENDED_MIXING_GATE_20261015.md").read_text(encoding="utf-8")
    for value in ("20261017","20261018","1,500-km","99","30×E","20×E","NOT_EVALUABLE_MCMC_CALIBRATION"):
        assert value in text

"""Unit guards for butterfly-first trait access analysis (no made-up ecology)."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

# scikit-learn is optional in the frozen independent GEB reviewer test environment.
# This exploratory trait test still runs fully in its dedicated source-locked CI.
pytest.importorskip("sklearn")

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/"scripts/run_butterfly_trait_marine_access_20261108.py"
spec=importlib.util.spec_from_file_location("traitpilot",S)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def test_frozen_trait_classification_and_unclassified_gap():
    assert mod.disturbance_class("Disturbance-associated (strong)")==1
    assert mod.disturbance_class("Disturbance-avoidant (strong)")==-1
    assert mod.disturbance_class("No evidence for disturbance associations")==0
    assert mod.disturbance_class("NA") is None

def test_jaccard_undefined_union_is_not_zero_similarity():
    assert mod.jac(set(),set()) is None
    assert mod.jac({"a"},{"b"})==0
    assert mod.jac({"a","b"},{"b","c"})==1/3

def test_logloss_brier_input_arrays():
    losses=mod.probability_loss(np.array([0,1]),np.array([0.2,0.8]))
    assert abs(losses["brier"]-0.04)<1e-12
    assert losses["logloss"]>0

def test_model_structure_preserves_necessary_main_effects():
    assert "min_log_uniform_effective_resistance" in mod.MODELS["W_SEA"]
    assert "wing_by_sea" in mod.MODELS["W_SEA"]
    assert "min_log_uniform_effective_resistance" in mod.MODELS["D_SEA"]
    assert "wing_by_port" in mod.MODELS["W_PORT"]
    assert "wing_by_latitude" in mod.MODELS["W_LAT"]

def test_exact_pinned_trait_source_checks():
    assert mod.git_blob((ROOT/mod.TRAITS[0]).read_bytes())==mod.TRAITS[1]
    assert mod.git_blob((ROOT/mod.RESOURCES[0]).read_bytes())==mod.RESOURCES[1]

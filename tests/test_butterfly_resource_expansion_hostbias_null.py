from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_butterfly_resource_expansion_hostbias_null.py"
SPEC = importlib.util.spec_from_file_location("hostbias_null", SCRIPT)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_hostbias_plan_preserves_count_and_excludes_observed_hosts():
    observed = {"h1": "F", "h2": "F"}
    family_pool = {"F": ("h1", "h2", "c1", "c2", "c3")}
    native = {
        "h1": frozenset({"A", "B"}),
        "h2": frozenset({"A", "B", "C"}),
        "c1": frozenset({"A", "B"}),
        "c2": frozenset({"A", "B", "C"}),
        "c3": frozenset({"A", "B", "C", "D"}),
    }
    consumers = {
        "h1": frozenset({"Butterfly focal"}),
        "h2": frozenset({"Butterfly focal", "Other one"}),
        "c1": frozenset({"Other one"}),
        "c2": frozenset({"Other one", "Other two"}),
        "c3": frozenset({"Other one", "Other two", "Other three"}),
    }

    plan = mod.build_sampling_plan(
        "Butterfly focal",
        observed,
        family_pool,
        native,
        consumers,
        bandwidth=0.20,
        usage_exponent=1.0,
    )
    assert len(plan) == 1
    assert set(plan[0]["alternatives"]) == {"c1", "c2", "c3"}
    assert len(plan[0]["targets"]) == 2

    for model in ("native_range", "usage", "native_range_usage"):
        sampled, native_difference, sampled_usage = mod.sample_from_plan(
            plan,
            np.random.default_rng(123),
            model,
        )
        assert len(sampled) == 2
        assert len(set(sampled)) == 2
        assert not (set(sampled) & set(observed))
        assert len(native_difference) == 2
        assert len(sampled_usage) == 2


def test_hostbias_portfolio_metrics_uses_union_and_added_units():
    native = {
        "h1": frozenset({"A", "B"}),
        "h2": frozenset({"B", "C"}),
    }
    contemporary = {
        "h1": frozenset({"A", "B", "D"}),
        "h2": frozenset({"B", "C", "E"}),
    }
    metrics = mod.portfolio_metrics(("h1", "h2"), native, contemporary)
    assert metrics["native_units"] == 3
    assert metrics["contemporary_units"] == 5
    assert metrics["added_units"] == 2
    assert metrics["log_expansion"] > 0


def test_usage_weights_are_other_lepidoptera_counts_plus_one():
    observed = {"h1": "F"}
    family_pool = {"F": ("h1", "c1", "c2")}
    native = {
        "h1": frozenset({"A"}),
        "c1": frozenset({"A"}),
        "c2": frozenset({"A"}),
    }
    consumers = {
        "h1": frozenset({"Butterfly focal"}),
        "c1": frozenset({"Other one"}),
        "c2": frozenset({"Other one", "Other two", "Other three"}),
    }
    plan = mod.build_sampling_plan(
        "Butterfly focal",
        observed,
        family_pool,
        native,
        consumers,
        bandwidth=0.20,
        usage_exponent=1.0,
    )
    weights = plan[0]["targets"][0]["weights"]["usage"]
    lookup = dict(zip(plan[0]["alternatives"], weights))
    assert lookup["c1"] == 2.0
    assert lookup["c2"] == 4.0

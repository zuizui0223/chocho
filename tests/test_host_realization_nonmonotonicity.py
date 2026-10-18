"""Source-independent mathematical controls for ecological monotonicity debate.

These are deterministic identities and hypothetical examples, NOT findings
of host-driven butterfly extinction, interaction impairment, or immigration.
"""
import math
import random

import pytest

from butterfly_specialization_ecology.host_realization_nonmonotonicity import (
    SeasonalEggExposure,
    additions_without_displacement,
    illustrative_cases,
    mapped_resource_cells,
    multi_year_growth,
    seasonally_realized_recruitment,
    source_sink_breeder_diversion,
)


def test_structural_host_geography_monotone_even_when_hosts_overlap():
    x=mapped_resource_cells(["native", "both"],["both","introduced"])
    assert x["native_cells"] == 2
    assert x["after_host_addition_cells"] == 3
    assert x["unique_new_cells"] == 1
    assert x["new_cell_keys"] == ["introduced"]
    assert x["structural_monotonicity"] is True


def test_absent_diversion_or_native_source_loss_prevents_negative_effect():
    for i in range(100):
        x=additions_without_displacement(i,12.2)
        assert x["after"] >= x["before"] and x["difference"]==12.2
    unchanged=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.95,0.02,0.0)
    ])
    assert unchanged["introduced_minus_native_recruits"]==0
    assert unchanged["effective_reproduction_drops"] is False


def test_lower_recruitment_requires_diverted_eggs_and_quality_deficit():
    x=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.9,0.1,0.3)
    ])
    assert x["baseline_recruits"]==90
    assert x["mixed_host_recruits"]==pytest.approx(66)
    assert x["introduced_minus_native_recruits"]==pytest.approx(-24)
    assert x["effective_reproduction_drops"] is True
    y=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.1,0.9,0.3)
    ])
    assert y["introduced_minus_native_recruits"]==pytest.approx(24)
    assert y["effective_reproduction_drops"] is False


def test_true_new_local_patch_cannot_have_negative_per_egg_yield():
    local=seasonally_realized_recruitment([
        SeasonalEggExposure(150,0.0,0.15,1.0)
    ])
    assert local["introduced_minus_native_recruits"]==22.5
    assert not local["effective_reproduction_drops"]


def test_analytical_delta_agrees_with_direct_counterfactual_over_parameter_grid():
    rng=random.Random(20261018)
    for _ in range(999):
        e=rng.random()*1000
        n=rng.random()
        i=rng.random()
        q=rng.random()
        ans=seasonally_realized_recruitment([
            SeasonalEggExposure(e,n,i,q)
        ])
        assert ans["introduced_minus_native_recruits"]==pytest.approx(
            e*q*(i-n),abs=1e-10
        )
        assert ans["mixed_host_recruits"]>=-1e-10
        assert ans["mixed_host_recruits"]<=e+1e-10


def test_seasonal_performance_can_switch_sign_across_windows():
    spring=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.5,0.8,0.6)
    ])
    autumn=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.7,0.2,0.6)
    ])
    year=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.5,0.8,0.6),
        SeasonalEggExposure(100,0.7,0.2,0.6),
    ])
    assert spring["introduced_minus_native_recruits"]==pytest.approx(18)
    assert autumn["introduced_minus_native_recruits"]==pytest.approx(-30)
    assert year["introduced_minus_native_recruits"]==pytest.approx(-12)
    assert year["effective_reproduction_drops"] is True


def test_productive_new_patch_can_still_be_landscape_sink_under_breeder_diversion():
    sink=source_sink_breeder_diversion(2.0,0.6,0.4,0.9)
    assert sink["new_patch_can_be_locally_productive"] is True
    assert sink["after_total_per_original_breeder"]==pytest.approx(1.416)
    assert sink["global_recruitment_drops"] is True
    nondisturb=source_sink_breeder_diversion(2.0,0.6,0.0,0.9)
    assert nondisturb["after_total_per_original_breeder"]==2.0
    highquality=source_sink_breeder_diversion(2.0,4.0,0.4,1.0)
    assert highquality["global_recruitment_drops"] is False


def test_temporal_growth_can_decline_even_if_arithmetic_mean_increases():
    native=multi_year_growth([0.8,0.8])
    introduced=multi_year_growth([1.28,0.326])
    assert introduced["arithmetic_mean_multiplier"] > native["arithmetic_mean_multiplier"]
    assert introduced["geometric_mean_multiplier"] < native["geometric_mean_multiplier"]
    assert introduced["geometric_mean_multiplier"]==pytest.approx(math.sqrt(1.28*.326))


def test_examples_label_is_categorical_no_species_or_pvalues():
    d=illustrative_cases()
    assert d["status"]=="MATHEMATICAL_THOUGHT_EXPERIMENT_NOT_EMPIRICAL_RESULT"
    assert d["species_observations"]==0
    assert "not fitted butterfly" in d["warning"]


@pytest.mark.parametrize("kw",[
    {"eggs":math.nan,"native_recruit_yield":.6,"introduced_recruit_yield":.3,"introduced_share":.8},
    {"eggs":-1,"native_recruit_yield":.6,"introduced_recruit_yield":.3,"introduced_share":.8},
    {"eggs":3,"native_recruit_yield":1.1,"introduced_recruit_yield":.3,"introduced_share":.8},
    {"eggs":3,"native_recruit_yield":.5,"introduced_recruit_yield":.3,"introduced_share":-0.1},
])
def test_invalid_recruitment_parameters_rejected(kw):
    with pytest.raises(ValueError):
        SeasonalEggExposure(**kw)


def test_invalid_annual_growth_and_displacement_parameters_rejected():
    with pytest.raises(ValueError):
        multi_year_growth([.8,0.0])
    with pytest.raises(ValueError):
        source_sink_breeder_diversion(2.0,0.6,1.5)
    with pytest.raises(ValueError):
        mapped_resource_cells(["native",""],["new"])

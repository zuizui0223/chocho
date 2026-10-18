"""Minimal, falsifiable distinction between host presence and larval recruitment.

Everything here is an explicit conditional model, NOT a fitted demographic
result. Fixed native-host plant maps are monotone under host additions; egg
reallocation, local egg-to-recruit quality and dispersal can break monotonic
recruitment. Seasons and annual growth are kept separate so there is no
accidental inference from mapped host presence to butterfly success.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log
from typing import Iterable


def _nonnegative_finite(label: str, value: float) -> None:
    if not isfinite(value) or value < 0:
        raise ValueError(f"{label} must be finite and non-negative")


def _fraction(label: str, value: float) -> None:
    if not isfinite(value) or not (0 <= value <= 1):
        raise ValueError(f"{label} must lie in [0,1]")


def mapped_resource_cells(native: Iterable[str], introduced: Iterable[str]) -> dict[str, object]:
    """Structural set monotonicity; NO measured local butterfly use."""
    native_set = set(native)
    introduced_set = set(introduced)
    if any(not isinstance(x, str) or not x.strip() for x in native_set | introduced_set):
        raise ValueError("all resource region keys must be non-empty strings")
    contemporary = native_set | introduced_set
    if not native_set.issubset(contemporary):
        raise ValueError("resource union is unexpectedly nonmonotone")
    return {
        "native_cells": len(native_set),
        "after_host_addition_cells": len(contemporary),
        "unique_new_cells": len(contemporary - native_set),
        "new_cell_keys": sorted(contemporary - native_set),
        "structural_monotonicity": True,
    }


@dataclass(frozen=True)
class SeasonalEggExposure:
    """A season's egg budget and COUNTERFACTUAL reallocation, not empirical rates.

    native_recruit_yield and introduced_recruit_yield are per-egg expected
    emergence probabilities (0–1) conditioned on field-specific plant quality.
    Introduced_share is the fraction of the SAME adult egg budget that would
    otherwise be laid on the best local native option. It is a diversion
    counterfactual, not a plant-preference score inferred from butterfly GBIF.
    """
    eggs: float
    native_recruit_yield: float
    introduced_recruit_yield: float
    introduced_share: float

    def __post_init__(self) -> None:
        _nonnegative_finite("eggs", self.eggs)
        for name in ("native_recruit_yield", "introduced_recruit_yield", "introduced_share"):
            _fraction(name, getattr(self, name))


def seasonally_realized_recruitment(windows: Iterable[SeasonalEggExposure]) -> dict[str, float | bool]:
    """Analytic equality: Δ recruits = Σ eggs*q*(introduced - native)."""
    days=list(windows)
    if not days:
        raise ValueError("at least one measured or hypothetical seasonal window needed")
    base=sum(w.eggs*w.native_recruit_yield for w in days)
    new=sum(
        w.eggs*((1-w.introduced_share)*w.native_recruit_yield
                +w.introduced_share*w.introduced_recruit_yield)
        for w in days
    )
    delta=sum(
        w.eggs*w.introduced_share*(w.introduced_recruit_yield-w.native_recruit_yield)
        for w in days
    )
    if abs((new-base)-delta)>1e-8*max(1,base,new):
        raise ValueError("weighted quality/diversion identity not conserved")
    if any(w.introduced_share>0 for w in days):
        diversion_active=True
    else:
        diversion_active=False
    return {
        "baseline_recruits":base,
        "mixed_host_recruits":new,
        "introduced_minus_native_recruits":delta,
        "effective_reproduction_drops":delta < -1e-10,
        "egg_diversion_exists":diversion_active,
    }


def additions_without_displacement(
    native_recruits: float, independently_added_recruits: float
) -> dict[str, float | bool]:
    """Nonnegative added recruitment cannot harm an unchanged native source."""
    _nonnegative_finite("native_recruits",native_recruits)
    _nonnegative_finite("independently_added_recruits",independently_added_recruits)
    final=native_recruits+independently_added_recruits
    return {
        "before":native_recruits,
        "after":final,
        "difference":independently_added_recruits,
        "monotone":True,
    }


def source_sink_breeder_diversion(
    native_recruits_per_breeder:float,
    new_patch_recruits_per_breeder:float,
    diverted_breeder_share:float,
    survived_travel_share:float=1.0,
) -> dict[str,float | bool]:
    """New patches can be locally productive but reduce total sources.

    Assumption: diverted breeders would otherwise reproduce at the native
    source; there is a fixed number of breeders and no compensation.
    """
    _nonnegative_finite("native_recruits_per_breeder",native_recruits_per_breeder)
    _nonnegative_finite("new_patch_recruits_per_breeder",new_patch_recruits_per_breeder)
    _fraction("diverted_breeder_share",diverted_breeder_share)
    _fraction("survived_travel_share",survived_travel_share)
    baseline=native_recruits_per_breeder
    new=(1-diverted_breeder_share)*baseline+(
        diverted_breeder_share*survived_travel_share*new_patch_recruits_per_breeder
    )
    return {
        "baseline_total_per_original_breeder":baseline,
        "after_total_per_original_breeder":new,
        "delta_per_breeder":new-baseline,
        "new_patch_can_be_locally_productive":new_patch_recruits_per_breeder>0,
        "global_recruitment_drops":new < baseline-1e-10,
        "assumes_finite_reallocated_breeders":True,
    }


def multi_year_growth(annual_positive_multipliers: Iterable[float]) -> dict[str,float]:
    """Long-term log-growth is not arithmetic mean; no fitted variance here."""
    xs=list(annual_positive_multipliers)
    if not xs:
        raise ValueError("requires annual rates")
    for r in xs:
        if not isfinite(r) or r<=0:
            raise ValueError("annual multipliers must be finite and > 0 for logarithmic rate")
    average=sum(xs)/len(xs)
    log_mean=sum(log(r) for r in xs)/len(xs)
    geometric=exp(log_mean)
    if geometric>average+1e-10*max(1,average):
        raise ValueError("Jensen inequality violated")
    return {
        "years":len(xs),
        "arithmetic_mean_multiplier":average,
        "geometric_mean_multiplier":geometric,
        "log_long_run_growth_rate":log_mean,
    }


def illustrative_cases() -> dict[str,object]:
    """Predeclared EXAMPLES chosen to show three logical regimes, NOT data."""
    early_and_late=seasonally_realized_recruitment([
        SeasonalEggExposure(100,0.5,0.8,0.6),
        SeasonalEggExposure(100,0.7,0.2,0.6),
    ])
    pure_addition=additions_without_displacement(120,18)
    new_sink=source_sink_breeder_diversion(2.0,0.6,0.4,0.9)
    # Arithmetic mean of introduced-host annual growth is 0.803 > 0.8;
    # geometric mean sqrt(1.28*0.326) ≈0.646 < 0.8.
    source_years=multi_year_growth([0.8,0.8])
    alien_years=multi_year_growth([1.28,0.326])
    return {
        "status":"MATHEMATICAL_THOUGHT_EXPERIMENT_NOT_EMPIRICAL_RESULT",
        "species_observations":0,
        "seasonal_preference_performance":early_and_late,
        "nondiversion_addition":pure_addition,
        "locally_positive_but_landscape_negative":new_sink,
        "stable_native_two_year_growth":source_years,
        "unstable_introduced_two_year_growth":alien_years,
        "warning":"Illustrative values selected to cross sign boundaries. They are not fitted butterfly field data, null distributions or evidence for widespread ecological traps.",
    }

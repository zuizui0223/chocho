"""Analytic two-period *hypothesis generator* for temporally supplied host plants.

Not fitted to butterfly observations. Demand and conversion from unmet food to
adult emergence are assumptions. Independent of frozen GEB manuscript inputs.
"""
from __future__ import annotations

from math import exp, isfinite, log
from typing import Mapping

TwoPeriods = tuple[float, float]
Demand = Mapping[str, TwoPeriods]


def _validate(demands: Demand) -> None:
    if set(demands) != {"A", "B"}:
        raise ValueError("Exactly two species named A and B are required")
    for species, pair in demands.items():
        if len(pair) != 2 or any(not isfinite(x) or x <= 0 for x in pair):
            raise ValueError(f"Each {species} demand needs two positive finite values")


def supply_regimes(demands: Demand, shortage: float) -> dict[str, TwoPeriods]:
    """Equal total supply in every arm and no shortage for either solo species."""
    _validate(demands)
    if not isfinite(shortage) or shortage < 0:
        raise ValueError("shortage must be finite and nonnegative")
    if shortage >= min(x for pair in demands.values() for x in pair):
        raise ValueError("shortage must be less than each species-period demand")
    early = demands["A"][0] + demands["B"][0]
    late = demands["A"][1] + demands["B"][1]
    return {
        "early_shortage": (early - shortage, late),
        "distributed_shortage": (early - shortage / 2, late - shortage / 2),
        "late_shortage": (early, late - shortage),
    }


def realized_food(demands: Demand, supply: TwoPeriods) -> dict[str, TwoPeriods]:
    """Proportional sharing of contemporaneous supply, with no carryover."""
    _validate(demands)
    if len(supply) != 2 or any(not isfinite(x) or x < 0 for x in supply):
        raise ValueError("supply must be finite and nonnegative")
    fractions = tuple(
        min(1.0, supply[t] / (demands["A"][t] + demands["B"][t]))
        for t in (0, 1)
    )
    return {s: (demands[s][0] * fractions[0], demands[s][1] * fractions[1])
            for s in ("A", "B")}


def fraction_unmet(demand: TwoPeriods, received: TwoPeriods) -> float:
    if any(r < 0 or r > d + 1e-12 for d, r in zip(demand, received)):
        raise ValueError("received must lie between zero and demand")
    return sum(max(0.0, d - r) for d, r in zip(demand, received)) / sum(demand)


def reversal_interval(demands: Demand) -> dict[str, float] | None:
    """Open k_A/k_B interval where A loses more early, B loses more late."""
    _validate(demands)
    a = demands["A"][0] / sum(demands["A"])
    b = demands["B"][0] / sum(demands["B"])
    lower, upper = b / a, (1 - b) / (1 - a)
    if lower >= upper:
        return None
    return {"lower_exclusive": lower, "upper_exclusive": upper}


def evaluate(
    demands: Demand,
    shortage: float,
    *,
    baseline_emergence: Mapping[str, float] | None = None,
    shortage_sensitivity: Mapping[str, float] | None = None,
) -> dict:
    """Compute hypothetical mixed-versus-solo log-emergence penalties.

    p_i = p0_i * exp(-k_i * unmet_i). This is an illustrative monotone
    response, NOT measured survival or a fitted butterfly fitness model.
    """
    p0 = dict(baseline_emergence) if baseline_emergence is not None else {"A": 0.8, "B": 0.8}
    k = dict(shortage_sensitivity) if shortage_sensitivity is not None else {"A": 3.0, "B": 3.0}
    if set(p0) != {"A", "B"} or any(not isfinite(x) or not 0 < x <= 1 for x in p0.values()):
        raise ValueError("baseline_emergence must give valid probabilities for A,B")
    if set(k) != {"A", "B"} or any(not isfinite(x) or x < 0 for x in k.values()):
        raise ValueError("shortage_sensitivity must be nonnegative finite for A,B")
    regimes = supply_regimes(demands, shortage)
    out: dict[str, dict] = {}
    for name, supply in regimes.items():
        mixed = realized_food(demands, supply)
        unmet = {s: fraction_unmet(demands[s], mixed[s]) for s in ("A", "B")}
        solo_unmet = {
            s: fraction_unmet(demands[s], (min(demands[s][0], supply[0]),
                                           min(demands[s][1], supply[1])))
            for s in ("A", "B")
        }
        solo_probability = {s: p0[s] * exp(-k[s] * solo_unmet[s]) for s in ("A", "B")}
        mixed_probability = {s: p0[s] * exp(-k[s] * unmet[s]) for s in ("A", "B")}
        penalty = {s: -log(mixed_probability[s] / solo_probability[s]) for s in ("A", "B")}
        out[name] = {
            "supply": list(supply),
            "mixed_unmet_fraction": unmet,
            "solo_unmet_fraction": solo_unmet,
            "illustrative_solo_emergence_probability": solo_probability,
            "illustrative_mixed_emergence_probability": mixed_probability,
            "mixed_vs_solo_log_penalty": penalty,
            "relative_disadvantage_A_minus_B": penalty["A"] - penalty["B"],
        }
    early = out["early_shortage"]["relative_disadvantage_A_minus_B"]
    late = out["late_shortage"]["relative_disadvantage_A_minus_B"]
    return {
        "status": "SYNTHETIC_THEORY_ONLY_NO_EMPIRICAL_FIT",
        "demands": {s: list(demands[s]) for s in ("A", "B")},
        "seasonal_supply_total": sum(out["early_shortage"]["supply"]),
        "shortage": shortage,
        "assumptions": [
            "two calendar periods, no carryover",
            "proportional food sharing within period",
            "unchanged host quality, predators, chemistry and development",
            "illustrative exponential mapping to emergence, not biological data",
        ],
        "regimes": out,
        "relative_disadvantage_reverses": early * late < 0,
        "inversion_sensitivity_interval_kA_over_kB": reversal_interval(demands),
    }

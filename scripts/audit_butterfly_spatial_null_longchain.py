#!/usr/bin/env python3
"""Extended-chain MCMC reliability audit of the 1,500-km spatial host null.

This is a finite post-hoc numerical QUALIFICATION attempt after a pilot
whose high autocorrelation invalidated naive permutation p interpretations.
No result can be promoted if the frozen mixing/ESS/drift gates fail.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from analyze_butterfly_shared_host_structural_null import build_model, read_csv
from analyze_butterfly_climate_matched_plant_range_null import (
    EXPECTED, read_unit_climate, climate_strata, evaluate,
)
from analyze_butterfly_spatially_local_plant_range_null import (
    load_geography, run_caliper, MAIN_METRIC, SECONDARY_METRIC,
)

CLIMATE_SHA256 = "0087117e0915a0847abaae8d11762b75f35b8eac94434d50913c2781bea727e2"
GEOGRAPHY_SHA256 = "c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840"
CALIPER_KM = 1500
SEEDS = [20261017, 20261018]
DRAW_COUNT = 99
BURN_FACTOR = 30
SPACING_FACTOR = 20


def integrated_acf_ess(series: list[float]) -> tuple[float | None, float | None, float | None]:
    """Initial-positive adjacent-lag-pair ESS diagnostic, NOT exact inference."""
    n = len(series)
    if n < 16:
        return None, None, None
    x = [float(z) for z in series]
    mu = sum(x)/n
    c = [z-mu for z in x]
    den = sum(z*z for z in c)
    if den <= 1e-10:
        return None, None, None

    def rho(k: int) -> float:
        return sum(c[i]*c[i+k] for i in range(n-k))/den
    ac1 = rho(1)
    # Truncate paired autocorrelation at first nonpositive sum.
    tau = 1.0
    for lag in range(1, n//2-1, 2):
        summed = rho(lag)+rho(lag+1)
        if summed <= 0:
            break
        tau += 2*summed
    return ac1, tau, n/tau


def mixing_check(series: list[float]) -> dict[str, object]:
    n = len(series)
    rho1, tau, ess = integrated_acf_ess(series)
    sd = statistics.stdev(series) if len(series)>1 else None
    q = n//4
    early = statistics.mean(series[:q])
    late = statistics.mean(series[-q:])
    drift = late-early
    gate = (
        rho1 is not None and ess is not None and sd is not None
        and rho1 <= 0.5 and ess >= 30
        and abs(drift) <= 0.5*sd
    )
    return {
        "n":n,"lag1":rho1,"iact":tau,"approximate_ess":ess,
        "full_sd":sd,"early_quarter_mean":early,
        "late_quarter_mean":late,"late_minus_early":drift,
        "lag1_gate":rho1 is not None and rho1 <= 0.5,
        "ess_gate":ess is not None and ess >= 30,
        "drift_gate":sd is not None and abs(drift) <= 0.5*sd,
        "pass":bool(gate),
    }


def qualify(out: dict[str, object]) -> dict[str, object]:
    ch = out.get("chains",[])
    checks = []
    for c in ch:
        if c.get("status")!="EXPLORATORY_SPATIAL_NULL_EVALUABLE":
            checks.append({
                "seed":c.get("seed"),"status":c.get("status"),"pass":False,
            })
            continue
        prim = mixing_check(c["primary_null_draws"])
        sec = mixing_check(c["secondary_null_draws"])
        checks.append({
            "seed":c["seed"],
            "status":c["status"],
            "primary_A":prim,
            "secondary_U":sec,
            "pass":prim["pass"] and sec["pass"],
        })
    between = None
    if len(checks)==2 and all(c["status"]=="EXPLORATORY_SPATIAL_NULL_EVALUABLE" for c in ch):
        medians = [
            float(c["outcomes"][MAIN_METRIC]["null_median"])
            for c in ch
        ]
        sigma = max(
            float(q["primary_A"]["full_sd"]) for q in checks
        )
        between = {
            "abs_difference_primary_medians":abs(medians[0]-medians[1]),
            "maximum_within_chain_sd":sigma,
            "threshold":0.5*sigma,
            "pass":abs(medians[0]-medians[1]) <= 0.5*sigma,
        }
    ok = (
        out.get("status")=="EVALUABLE"
        and len(checks)==2
        and all(x["pass"] for x in checks)
        and between is not None and between["pass"]
    )
    return {
        "status":(
            "QUALIFIED_EXPLORATORY_MCMC_SENSITIVITY"
            if ok else "NOT_EVALUABLE_MCMC_CALIBRATION"
        ),
        "checks":checks,
        "between_chain":between,
        "interpretation":(
            "The outcome was screened by predetermined post-hoc simulation "
            "reliability gates. Passing does not convert it to prospectively "
            "tested biology; failing precludes inferential promotion."
        ),
    }


def main() -> None:
    p=argparse.ArgumentParser()
    for k in (
        "descriptors_csv","insect_host_csv","native_distribution_csv",
        "contemporary_distribution_csv","level1_json",
        "unit_climate_csv","level3_geojson","output_json",
    ):
        p.add_argument("--"+k.replace("_","-"),type=Path,required=True)
    args=p.parse_args()
    lookup=json.loads(args.level1_json.read_text(encoding="utf-8"))["units"]
    nrows=read_csv(args.native_distribution_csv)
    crows=read_csv(args.contemporary_distribution_csv)
    hrows=read_csv(args.insect_host_csv)
    model=build_model(
        read_csv(args.descriptors_csv),hrows,nrows,crows,lookup
    )
    focal=set(model["host_families"])
    codes=sorted({
        r["area_code_l3"] for r in nrows+crows
        if r["accepted_plant_name_id"] in focal
    })
    if len(codes)!=365 or len(codes)!=len(model["level1"]):
        raise ValueError("historical 365 region reconstruction drift")
    if any(lookup[u]!=model["level1"][i] for i,u in enumerate(codes)):
        raise ValueError("historical Level1 code drift")
    climate, climate_source=read_unit_climate(args.unit_climate_csv)
    positions, geo_source=load_geography(
        args.level3_geojson,codes,model["level1"]
    )
    if climate_source["source_file_sha256"]!=CLIMATE_SHA256:
        raise ValueError("previous CHELSA climate input hash changed")
    if geo_source["source_geometry_sha256"]!=GEOGRAPHY_SHA256:
        raise ValueError("previous TDWG polygon input hash changed")
    if len(positions)!=365 or not set(codes)<=set(climate):
        raise ValueError("climate/geometry resource coverage failure")
    strata, bin_info=climate_strata(
        climate,codes,model["level1"],
        "temperature_precipitation_2x2"
    )
    if len(strata)!=365:
        raise ValueError("climate-binning support lost")
    obs=evaluate(model,model["introduced_sets"])
    if any(obs[k]!=v for k,v in EXPECTED.items()):
        raise ValueError("baseline resource effect drift")
    data=run_caliper(
        model,strata,positions,CALIPER_KM,SEEDS,DRAW_COUNT,
        BURN_FACTOR,SPACING_FACTOR
    )
    qual=qualify(data)
    output={
        "schema":"chocho_butterfly_spatial_null_extended_mixing_v0.2",
        "status":qual["status"],
        "source_climate":climate_source,
        "source_geography":geo_source,
        "study":{"focal_butterflies":239,"wgsrpd3_regions":365},
        "protocol":"provenance/reviewer_defenses/protocols/SPATIAL_NULL_EXTENDED_MIXING_GATE_20261015.md",
        "frozen_protocol":{
            "primary_caliper_km":CALIPER_KM,
            "seeds":SEEDS,"draws_per_seed":DRAW_COUNT,
            "burn_factor":BURN_FACTOR,"spacing_factor":SPACING_FACTOR,
            "climate_bins":bin_info["mode"],
            "minimum_lag1_max":0.50,
            "minimum_ess":30,
            "quarter_drift_max_in_full_sd":0.5,
            "between_chain_median_difference_max_in_sd":0.5,
        },
        "observed":obs,
        "null":data,
        "qualification":qual,
        "scientific_boundary":(
            "Caveated post-hoc geographic sensitivity of mapped larval "
            "resource opportunity, not independent butterfly recruitment, "
            "colonization, or demographic effects. Nominal chain p-values "
            "are not a valid exact permutation test even if qualified."
        ),
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(
        json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    print(json.dumps({
        "schema":output["schema"],"status":output["status"],
        "observed":{k:obs[k] for k in (MAIN_METRIC,SECONDARY_METRIC)},
        "geometry":data["geometry"],
        "qualification":qual,
        "chains":[
            {
                "seed":c["seed"],
                "status":c["status"],
                "observed_absorption_minus_null":c.get("outcomes",{}).get(
                    MAIN_METRIC,{}).get("observed_minus_null_median"
                ),
                "null_absorption_median":c.get("outcomes",{}).get(
                    MAIN_METRIC,{}).get("null_median"
                ),
                "observed_unique_added_minus_null":c.get("outcomes",{}).get(
                    SECONDARY_METRIC,{}).get("observed_minus_null_median"
                ),
            } for c in data.get("chains",[])
        ],
    },indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

import analyze_butterfly_resource_expansion_hostbias_null as base


def build_focal(descriptors, pairs, family_pool, native, contemporary, consumers, bandwidth, exponent):
    focal=[]
    excluded={}
    valid={
        host
        for host, units in native.items()
        if units and units <= contemporary.get(host, frozenset())
    }
    pool={
        family: tuple(host for host in hosts if host in valid)
        for family, hosts in family_pool.items()
    }

    for species, descriptor in descriptors.items():
        if (
            descriptor["host_family_count"] <= 0
            or descriptor["native_resource_units"] <= 0
            or descriptor["resolved_host_species"] < descriptor["host_family_count"]
        ):
            continue
        host_map=pairs.get(species,{})
        if not host_map:
            excluded[species]="no_resolved_hosts"
            continue
        missing=[h for h in host_map if not native.get(h,frozenset())]
        if missing:
            excluded[species]=f"observed_hosts_without_native_footprint:{len(missing)}"
            continue
        comp=Counter(host_map.values())
        bad=[
            fam for fam,n in comp.items()
            if len([h for h in pool.get(fam,()) if h not in host_map]) < n
        ]
        if bad:
            excluded[species]="insufficient_alternative_pool:"+",".join(sorted(bad))
            continue
        observed=base.portfolio_metrics(tuple(host_map),native,contemporary)
        if observed is None:
            excluded[species]="no_native_union"
            continue
        focal.append({
            "species":species,
            "descriptor":descriptor,
            "host_map":host_map,
            "observed":observed,
            "sampling_plan":base.build_sampling_plan(
                species,host_map,pool,native,consumers,bandwidth,exponent
            ),
        })
    return focal, excluded


def observed_usage_mean(focal, consumers):
    vals=[
        len(consumers.get(host,frozenset())-{row["species"]})
        for row in focal
        for host in row["host_map"]
    ]
    return float(np.mean(vals))


def simulate_mean_usage(focal, replicates, seed_tag):
    means=[]
    for rep in range(replicates):
        vals=[]
        for row in focal:
            rng=np.random.default_rng(
                base.stable_seed(seed_tag,"native_range_usage",rep,row["species"])
            )
            _,_,sampled_usage=base.sample_from_plan(
                row["sampling_plan"],rng,"native_range_usage"
            )
            vals.extend(sampled_usage)
        means.append(float(np.mean(vals)))
    return np.asarray(means,float)


def calibrate_exponent(
    descriptors,pairs,family_pool,native,contemporary,consumers,
    bandwidth,target,replicates,iterations,seed_tag
):
    low,high=0.0,3.0
    trace=[]

    def evaluate(exp):
        focal,_=build_focal(
            descriptors,pairs,family_pool,native,contemporary,consumers,
            bandwidth,exp
        )
        vals=simulate_mean_usage(
            focal,replicates,f"{seed_tag}|cal|bw={bandwidth}|exp={exp:.8f}"
        )
        med=float(np.median(vals))
        return med,focal

    low_value,_=evaluate(low)
    high_value,_=evaluate(high)
    if not (low_value <= target <= high_value):
        raise RuntimeError(
            f"target prominence {target} not bracketed for bandwidth {bandwidth}: "
            f"{low_value}..{high_value}"
        )

    chosen_focal=None
    for _ in range(iterations):
        mid=(low+high)/2
        value,focal=evaluate(mid)
        trace.append({"exponent":mid,"median_sampled_mean_degree":value})
        chosen_focal=focal
        if value < target:
            low=mid
        else:
            high=mid

    exponent=(low+high)/2
    final_value,chosen_focal=evaluate(exponent)
    return exponent,final_value,trace,chosen_focal


def simulate_expansion(focal, replicates, seed_tag):
    total_added=[]
    mean_log=[]
    median_log=[]
    sampled_degree=[]
    native_diff=[]
    family_breadth=[
        row["descriptor"]["host_family_count"] for row in focal
    ]
    rho=[]

    for rep in range(replicates):
        metrics=[]
        usage=[]
        diffs=[]
        for row in focal:
            rng=np.random.default_rng(
                base.stable_seed(seed_tag,"native_range_usage",rep,row["species"])
            )
            sampled,d,deg=base.sample_from_plan(
                row["sampling_plan"],rng,"native_range_usage"
            )
            m=base.portfolio_metrics(sampled,NATIVE,CONTEMPORARY)
            metrics.append(m)
            usage.extend(deg)
            diffs.extend(d)
        logs=[m["log_expansion"] for m in metrics]
        total_added.append(sum(m["added_units"] for m in metrics))
        mean_log.append(float(np.mean(logs)))
        median_log.append(float(np.median(logs)))
        sampled_degree.append(float(np.mean(usage)))
        native_diff.append(float(np.median(diffs)))
        rho.append(base.spearman(family_breadth,logs))

    return {
        "total_added":total_added,
        "mean_log":mean_log,
        "median_log":median_log,
        "sampled_degree":sampled_degree,
        "native_diff":native_diff,
        "rho":rho,
    }


NATIVE={}
CONTEMPORARY={}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--protocol-json",type=Path,required=True)
    ap.add_argument("--descriptors-csv",type=Path,required=True)
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--native-distribution-csv",type=Path,required=True)
    ap.add_argument("--contemporary-distribution-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    protocol=json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema")!="chocho_butterfly_prominence_calibrated_null_protocol_v0.1":
        raise RuntimeError("unexpected calibrated-null protocol")

    descriptors=base.load_descriptors(args.descriptors_csv)
    pairs,consumers,family_pool=base.load_pairs(args.insect_host_csv)
    global NATIVE,CONTEMPORARY
    NATIVE=base.load_units(args.native_distribution_csv)
    CONTEMPORARY=base.load_units(args.contemporary_distribution_csv)

    bandwidths=[float(x) for x in protocol["native_log_bandwidths"]]
    calibration_reps=int(protocol["calibration_replicates"])
    final_reps=int(protocol["final_permutations"])
    iterations=int(protocol["bisection_iterations"])
    seed_tag=str(protocol["seed_tag"])

    # Exponent 1 focal set defines eligibility; eligibility itself is exponent-invariant.
    baseline_focal,excluded=build_focal(
        descriptors,pairs,family_pool,NATIVE,CONTEMPORARY,consumers,
        bandwidths[0],1.0
    )
    target=observed_usage_mean(baseline_focal,consumers)
    observed={
        "species":len(baseline_focal),
        "mean_other_lepidoptera_users_per_host":target,
        "total_added_units":sum(r["observed"]["added_units"] for r in baseline_focal),
        "mean_log_expansion":float(np.mean([r["observed"]["log_expansion"] for r in baseline_focal])),
        "median_log_expansion":float(np.median([r["observed"]["log_expansion"] for r in baseline_focal])),
        "rho_host_family_vs_log_expansion":base.spearman(
            [r["descriptor"]["host_family_count"] for r in baseline_focal],
            [r["observed"]["log_expansion"] for r in baseline_focal],
        ),
    }

    results=[]
    for bw in bandwidths:
        exp,cal_value,trace,focal=calibrate_exponent(
            descriptors,pairs,family_pool,NATIVE,CONTEMPORARY,consumers,
            bw,target,calibration_reps,iterations,seed_tag
        )
        sim=simulate_expansion(
            focal,final_reps,f"{seed_tag}|final|bw={bw}|exp={exp:.8f}"
        )
        rho_arr=np.asarray(sim["rho"],float)
        results.append({
            "native_log_bandwidth":bw,
            "calibrated_usage_exponent":exp,
            "target_observed_mean_degree":target,
            "calibrated_median_sampled_mean_degree":cal_value,
            "calibration_trace":trace,
            "matching_diagnostics":{
                "final_sampled_mean_degree_median":float(np.median(sim["sampled_degree"])),
                "final_sampled_mean_degree_q025":float(np.quantile(sim["sampled_degree"],0.025)),
                "final_sampled_mean_degree_q975":float(np.quantile(sim["sampled_degree"],0.975)),
                "median_absolute_log1p_native_breadth_difference":float(np.median(sim["native_diff"])),
            },
            "total_introduced_added_units":base.null_summary(
                sim["total_added"],observed["total_added_units"]
            ),
            "mean_log_expansion":base.null_summary(
                sim["mean_log"],observed["mean_log_expansion"]
            ),
            "median_log_expansion":base.null_summary(
                sim["median_log"],observed["median_log_expansion"]
            ),
            "rho_host_family_vs_log_expansion":{
                "q025":float(np.quantile(rho_arr,0.025)),
                "median":float(np.median(rho_arr)),
                "q975":float(np.quantile(rho_arr,0.975)),
                "p_two_sided":float(
                    (1+np.sum(np.abs(rho_arr)>=abs(observed["rho_host_family_vs_log_expansion"])))
                    /(len(rho_arr)+1)
                ),
            },
        })

    payload={
        "schema":"chocho_butterfly_prominence_calibrated_null_v0.1",
        "status":"POSTHOC_OUTCOME_BLIND_PROMINENCE_CALIBRATED_NULL",
        "protocol":protocol["schema"],
        "observed":observed,
        "excluded_species":excluded,
        "results":results,
        "interpretation_rule":(
            "Because the use-weight exponent is selected solely to reproduce the observed "
            "mean other-Lepidoptera consumer degree, expansion outcomes do not influence "
            "calibration. If observed expansion is not above the calibrated null across "
            "native-range bandwidths, there is no identifiable butterfly-specific "
            "host-identity excess beyond matched plant prominence and starting geography."
        ),
        "claim_boundary":protocol["claim_boundary"],
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2,sort_keys=True))


if __name__=="__main__":
    main()

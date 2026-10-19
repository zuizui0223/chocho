#!/usr/bin/env python3
"""Macroecological trait filter of potential butterfly host-resource geography.

Uses ORIGINAL frozen chocho 239-species HOSTS/WCVP resource reconstructions
and the independent 339-species LepTraits-based habitat descriptor panel.
It does NOT estimate butterfly occupancy, recruitment, or causal plant effects.

The protocol and all claims are exploratory, not confirmatory:
docs/exploratory/GLOBAL_BUTTERFLY_HABITAT_FILTER_PROTOCOL_20261019.md

Dependencies: numpy >=1.26 (already required by chocho).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESOURCE = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
TRAITS = ROOT / "data/frozen/figure_sources/s1_resource_descriptors.csv"
LEPTRAITS = ROOT / "data/external/leptraits_consensus_v1.0.csv"
N_REGIONS = 365
SEED = 20261019
STRONG_ASSOC = "Disturbance-associated (strong)"
STRONG_AVOID = "Disturbance-avoidant (strong)"
BROAD_ASSOC = (STRONG_ASSOC, "Disturbance-associated (weak)")
BROAD_AVOID = (STRONG_AVOID, "Disturbance-avoidant (weak)")


def _read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def _digest(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_inputs():
    res, habitat, lep = map(_read_csv, (RESOURCE, TRAITS, LEPTRAITS))
    assert len(res) == 239, f"unexpected eligible resource panel {len(res)}"
    assert len(habitat) == 339, f"unexpected original trait source population {len(habitat)}"
    assert len({x["species"] for x in res}) == len(res)
    assert len({x["species"] for x in habitat}) == len(habitat)
    taxon = {}
    for line in lep:
        s = str(line["Species"]).strip()
        f = str(line["Family"]).strip()
        if s and f:
            taxon[s] = f
    trait_by_name = {row["species"]: row for row in habitat}
    assert all(x["species"] in trait_by_name and x["species"] in taxon for x in res)
    rows = []
    for row in res:
        s = row["species"]
        a = float(row["native_resource_units"])
        c = float(row["contemporary_resource_units"])
        added = float(row["introduced_added_units"])
        assert 0 < a <= c <= N_REGIONS and abs((c-a)-added) < 1e-8, (s,a,c,added)
        actual_y = math.log1p(c) - math.log1p(a)
        assert abs(actual_y - float(row["log_resource_expansion"])) < 1e-10, s
        n_hosts = float(row["resolved_host_species"])
        n_families = float(row["host_family_count"])
        wing = float(row["wing_size_proxy"])
        assert n_hosts > 0 and n_families > 0 and wing > 0
        trait = trait_by_name[s]
        rows.append({
            "species": s, "genus": s.split(" ")[0], "family": taxon[s],
            "disturbance": trait["disturbance_affinity"],
            "edge": trait["edge_affinity"],
            "canopy": trait["canopy_affinity"],
            "voltinism": row["voltinism"],
            "native": a, "contemporary": c, "added": added,
            "host_n": n_hosts, "host_families": n_families, "wing": wing,
            "adequate": row["host_taxonomy_lower_bound_adequate"] == "True",
            "log_expansion": actual_y,
            "filled_fraction": (added / (N_REGIONS - a)) if a < N_REGIONS else None,
            "capacity_logit": math.log((added + 0.5) / (N_REGIONS - a - added + 0.5))
            if a < N_REGIONS else None,
        })
    assert round(sum(r["native"] for r in rows)) == 26530
    assert round(sum(r["contemporary"] for r in rows)) == 41083
    assert round(sum(r["added"] for r in rows)) == 14553
    assert len({r["family"] for r in rows}) == 6
    return rows, habitat, {
        "resource_sha256": _digest(RESOURCE),
        "trait_source_sha256": _digest(TRAITS),
        "leptraits_source_sha256": _digest(LEPTRAITS),
    }


def mean(values):
    return float(np.mean(np.asarray(values, dtype=float)))


def median(values):
    return float(np.median(np.asarray(values, dtype=float)))


def describe(rows):
    if not rows:
        return {"species": 0}
    y = [r["log_expansion"] for r in rows]
    filled = [r["filled_fraction"] for r in rows if r["filled_fraction"] is not None]
    return {
        "species": len(rows), "genera": len({r["genus"] for r in rows}),
        "butterfly_families": len({r["family"] for r in rows}),
        "expanded": sum(r["added"] > 0 for r in rows),
        "mean_log_resource_expansion": mean(y),
        "median_log_resource_expansion": median(y),
        "mean_filled_fraction": mean(filled) if filled else None,
        "median_native_regions": median([r["native"] for r in rows]),
        "median_resolved_host_species": median([r["host_n"] for r in rows]),
        "median_host_family_richness": median([r["host_families"] for r in rows]),
        "adequate_host_inventory_fraction": mean([int(r["adequate"]) for r in rows]),
    }


def edge_class(edge):
    if edge.startswith("Edge-associated"):
        return "edge-associated"
    if edge.startswith("Edge-avoidant"):
        return "edge-avoidant"
    return "other"


def _columns(row, family_levels, *, extended_edge=False):
    # Fixed transformations need no information from held-out labels or outcomes.
    n = math.log(row["native"] / 100.0)
    c = [
        1.0,
        n, n*n,
        math.log1p(row["host_n"]), math.log1p(row["host_families"]),
        math.log(row["wing"] / 5.0),
        float(row["voltinism"] == "B"), float(row["voltinism"] == "U"),
        float(not row["adequate"]),
    ]
    c.extend(float(row["family"] == k) for k in family_levels[1:])
    if extended_edge:
        e = edge_class(row["edge"])
        c.extend([float(e == "edge-associated"), float(e == "edge-avoidant")])
    return c


def model_matrices(rows, all_families, *, y_key="log_expansion", extended_edge=False):
    z = np.array([_columns(r, all_families, extended_edge=extended_edge)
                  for r in rows], dtype=float)
    # Drop invariant/redundant CONTROL columns only, e.g. the all-adequate
    # subgroup's inventory-missing dummy. Exposure is NEVER dropped.
    keep = []
    for i in range(z.shape[1]):
        indices = keep + [i]
        if np.linalg.matrix_rank(z[:,indices]) > len(keep):
            keep.append(i)
    z = z[:,keep]
    x = np.array([float(r["exposure"]) for r in rows], dtype=float)[:,None]
    y = np.array([float(r[y_key]) for r in rows], dtype=float)
    return z, np.column_stack([z, x]), y


def fit_ols(z, x, y):
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    resid = y - x @ beta
    p = int(np.linalg.matrix_rank(x))
    return beta, resid, p


def _coefficient(rows, families, *, y_key="log_expansion", extended_edge=False):
    z, full, y = model_matrices(rows, families, y_key=y_key, extended_edge=extended_edge)
    beta, residuals, p = fit_ols(z, full, y)
    if p < full.shape[1]:
        raise ValueError("collinear design; disallow unidentifiable adjusted claim")
    return float(beta[-1]), beta, residuals, full, z, y


def cluster_se(full, residuals, groups):
    n, p = full.shape
    assert n > p and np.linalg.matrix_rank(full) == p
    labels = sorted(set(groups))
    G = len(labels)
    assert G > 10
    xx_inv = np.linalg.pinv(full.T @ full, rcond=1e-10)
    meat = np.zeros((p,p))
    for group in labels:
        idx = np.array([i for i,g in enumerate(groups) if g==group])
        score = full[idx].T @ residuals[idx]
        meat += np.outer(score,score)
    correction = (G / (G-1)) * ((n-1)/(n-p))
    cov = correction * xx_inv @ meat @ xx_inv
    return float(math.sqrt(max(0.0,cov[-1,-1]))), G


def wild_genus_signflip(z, y, x_full, genera, *, seed=SEED, n=4999):
    # Frisch-Waugh-Lovell residualized-exposure coefficient, under a
    # reduced-model residual wild sign-flip null, genus as the block.
    feature = x_full[:,-1]
    rx = feature - z @ np.linalg.lstsq(z, feature, rcond=None)[0]
    ry = y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
    denom = float(rx @ rx)
    if denom <= 1e-9:
        return {"status": "UNIDENTIFIABLE_WITHIN_CONTROLS"}
    value = float((rx @ ry) / denom)
    groups = sorted(set(genera))
    members = defaultdict(list)
    for i,g in enumerate(genera):
        members[g].append(i)
    rng = np.random.default_rng(seed)
    score_by_group = np.array([
        float(np.dot(rx[members[g]],ry[members[g]]))/denom
        for g in groups
    ])
    sim = rng.choice((-1.0,1.0), size=(n,len(groups))) @ score_by_group
    p = (int(np.sum(np.abs(sim) >= abs(value)-1e-12)) + 1) / (n + 1)
    return {
        "status": "EXPLORATORY_GENUS_WILD_SIGNFLIP_REDUCED_RESIDUAL",
        "residualized_beta": value, "simulation_draws": n,
        "two_sided_tail_fraction_not_randomization_p": p,
    }


def genus_bootstrap(rows, families, *, seed=SEED, nboot=2499, y_key="log_expansion"):
    members = defaultdict(list)
    for r in rows:
        members[r["genus"]].append(r)
    genera = sorted(members)
    rng = np.random.default_rng(seed)
    arr = []
    failures = 0
    for i in range(nboot):
        chosen = rng.integers(0,len(genera),len(genera))
        sample = [r for j in chosen for r in members[genera[int(j)]]]
        if len({r["exposure"] for r in sample}) != 2:
            failures += 1
            continue
        try:
            beta, _, _, _, _, _ = _coefficient(sample,families,y_key=y_key)
            if math.isfinite(beta):
                arr.append(beta)
            else:
                failures += 1
        except ValueError:
            failures += 1
    if len(arr) < .95*nboot:
        raise RuntimeError(f"GENUS_BOOTSTRAP_NOT_IDENTIFIABLE: {len(arr)}/{nboot}")
    return {
        "n_draws_requested": nboot,
        "n_draws_valid": len(arr),
        "invalid_rank_or_missing_class": failures,
        "percentile_95ci": [float(np.quantile(arr,.025)),float(np.quantile(arr,.975))],
        "median": float(np.median(arr)),
    }


def adjusted(rows, families, *, label, y_key="log_expansion", bootstrap=True, edge=False):
    b, coeffs, resid, X, z, y = _coefficient(rows,families,y_key=y_key,extended_edge=edge)
    err,G = cluster_se(X,resid,[r["genus"] for r in rows])
    out = {
        "label": label, "n": len(rows), "genera": G, "parameter": "associated_strong_minus_avoidant_strong" if
        label.startswith("strict") else "associated_minus_avoidant",
        "outcome": y_key,
        "adjusted_coefficient": b,
        "genus_cluster_se_CR1": err,
        "approx_genus_cluster_95ci": [b-1.96*err,b+1.96*err],
        "OLS_residual_mse": float(np.mean(resid**2)),
        "design_rank": int(np.linalg.matrix_rank(X)),
        "n_covariates_including_exposure": X.shape[1],
        "wild_reduced_genus_signflip": wild_genus_signflip(z,y,X,
            [r["genus"] for r in rows],seed=SEED+len(rows)),
    }
    if bootstrap and not edge:
        out["genus_bootstrap"] = genus_bootstrap(rows,families,
             seed=SEED+len(rows)+7,y_key=y_key)
    return out


def grouped_out_of_genus_cv(rows, families, *, repeat=10, folds=5, y_key="log_expansion"):
    members = defaultdict(list)
    for row in rows:
        members[row["genus"]].append(row)
    genera = sorted(members)
    assert len(genera) > 5 and len(rows)>100
    scores = []
    for rep in range(repeat):
        rng = np.random.default_rng(SEED + 541*rep)
        order = list(rng.permutation(len(genera)))
        bins=[[] for _ in range(folds)]
        fold_sizes=[0]*folds
        for i in order:
            g=genera[int(i)]
            smallest=min(range(folds),key=lambda j:(fold_sizes[j],j))
            bins[smallest].append(g)
            fold_sizes[smallest]+=len(members[g])
        actual = np.array([r[y_key] for r in rows])
        pred0 = np.zeros(len(rows));pred1=np.zeros(len(rows))
        for fold in bins:
            test_groups=set(fold)
            ii_test=[i for i,r in enumerate(rows) if r["genus"] in test_groups]
            ii_train=[i for i,r in enumerate(rows) if r["genus"] not in test_groups]
            assert not ({rows[i]["genus"] for i in ii_train} &
                        {rows[i]["genus"] for i in ii_test})
            train=[rows[i] for i in ii_train]
            test=[rows[i] for i in ii_test]
            # For genuine held-out prediction, the exact SAME control columns
            # must be used in train and test, even when a family occurs only
            # in the held-out fold. Least squares handles absent train levels
            # via minimum-norm coefficients; no test-label information leaks.
            z_train = np.array([_columns(r,families) for r in train],dtype=float)
            z_test = np.array([_columns(r,families) for r in test],dtype=float)
            X_train = np.column_stack([z_train,[r["exposure"] for r in train]])
            X_test = np.column_stack([z_test,[r["exposure"] for r in test]])
            y_train = np.array([r[y_key] for r in train],dtype=float)
            b0=np.linalg.lstsq(z_train,y_train,rcond=None)[0]
            b1=np.linalg.lstsq(X_train,y_train,rcond=None)[0]
            pred0[ii_test] = z_test @ b0
            pred1[ii_test] = X_test @ b1
        mse0=float(np.mean((actual-pred0)**2))
        mse1=float(np.mean((actual-pred1)**2))
        mae0=float(np.mean(abs(actual-pred0)))
        mae1=float(np.mean(abs(actual-pred1)))
        scores.append({
            "repeat":rep,"MSE_baseline":mse0,"MSE_plus_disturbance":mse1,
            "MAE_baseline":mae0,"MAE_plus_disturbance":mae1,
            "relative_MSE_reduction":(mse0-mse1)/mse0,
        })
    return {
        "heldout_grouping": "butterfly_genus_5fold",
        "outcome":y_key,
        "repeats":repeat,
        "mean_MSE_baseline":mean([x["MSE_baseline"] for x in scores]),
        "mean_MSE_plus_disturbance":mean([x["MSE_plus_disturbance"] for x in scores]),
        "mean_MAE_baseline":mean([x["MAE_baseline"] for x in scores]),
        "mean_MAE_plus_disturbance":mean([x["MAE_plus_disturbance"] for x in scores]),
        "mean_relative_MSE_reduction":mean([x["relative_MSE_reduction"] for x in scores]),
        "repeats_with_MSE_improvement":sum(x["relative_MSE_reduction"]>0 for x in scores),
        "individual_repeats":scores,
    }


def taxonomy_family_signs(rows):
    by_f = defaultdict(list)
    for r in rows:
        by_f[r["family"]].append(r)
    out={}
    for family, taxa in sorted(by_f.items()):
        a=[t for t in taxa if t["exposure"]==1]
        b=[t for t in taxa if t["exposure"]==0]
        out[family]={
            "n_affiliated":len(a), "n_avoidant":len(b),
            "unadjusted_log_mean_difference":
                mean([r["log_expansion"] for r in a]) -
                mean([r["log_expansion"] for r in b]) if a and b else None,
            "usable_stratum_n_ge_3_per_class":len(a)>=3 and len(b)>=3,
        }
    return out


def within_mixed_genus(rows, *, y_key="log_expansion"):
    by_genus=defaultdict(list)
    for r in rows:
        by_genus[r["genus"]].append(r)
    both={}
    for genus,taxa in by_genus.items():
        a=[r for r in taxa if r["exposure"]==1]
        b=[r for r in taxa if r["exposure"]==0]
        if a and b:
            both[genus]={
                "affiliated_n":len(a),
                "avoidant_n":len(b),
                "mean_log_difference":mean([r[y_key] for r in a]) -
                    mean([r[y_key] for r in b]),
            }
    return {
        "mixed_genera_n":len(both),
        "mixed_genus_species_n":sum(x["affiliated_n"]+x["avoidant_n"] for x in both.values()),
        "unweighted_genus_mean_within_difference": mean([x["mean_log_difference"] for x in both.values()]) if both else None,
        "genera_positive":sum(x["mean_log_difference"]>0 for x in both.values()),
        "genera_negative":sum(x["mean_log_difference"]<0 for x in both.values()),
        "genera": dict(sorted(both.items())),
        "warning":"Seven small mixed genera are not seven geographic replications; this is a conservative descriptive phylogenetic check, not a new population effect."
    }


def annotate(rows, *, broad=False):
    if broad:
        chosen=[r.copy() for r in rows if r["disturbance"] in BROAD_ASSOC+BROAD_AVOID]
        for r in chosen:
            r["exposure"]=int(r["disturbance"] in BROAD_ASSOC)
    else:
        chosen=[r.copy() for r in rows if r["disturbance"] in (STRONG_ASSOC,STRONG_AVOID)]
        for r in chosen:
            r["exposure"]=int(r["disturbance"]==STRONG_ASSOC)
    return chosen


def selection_audit(rows,trait_source):
    available={r["species"] for r in rows}
    categories=sorted({x["disturbance_affinity"] for x in trait_source})
    return {
        "all_trait_candidate_species": len(trait_source),
        "resource_eligible_species": len(rows),
        "by_disturbance": {
            cat: {
                "trait_candidates": sum(x["disturbance_affinity"]==cat for x in trait_source),
                "in_239_resource_panel": sum(r["disturbance"]==cat for r in rows),
                "retained_share": (
                    sum(r["disturbance"]==cat for r in rows) /
                    sum(x["disturbance_affinity"]==cat for x in trait_source)
                ),
            } for cat in categories
        },
        "bias_boundary": "Eligibility is disproportionately lower for disturbance-avoidant groups: missing response is not zero, not missing at random. No inference to globally all butterfly species."
    }


def main_analysis():
    all_rows, trait_source, digests = read_inputs()
    strong = annotate(all_rows)
    assert len(strong)==183
    assert sum(x["exposure"] for x in strong)==150
    assert len(strong)-sum(x["exposure"] for x in strong)==33
    full_families=sorted({r["family"] for r in all_rows})
    unadjusted_mean_diff=mean([r["log_expansion"] for r in strong if r["exposure"]==1])-mean(
        [r["log_expansion"] for r in strong if r["exposure"]==0])
    main=adjusted(strong,full_families,label="strict_strong",bootstrap=True)
    edge=adjusted(strong,full_families,label="strict_plus_edge",bootstrap=False,edge=True)
    ceiling_rows=[r for r in strong if r["capacity_logit"] is not None]
    ceiling=adjusted(ceiling_rows,full_families,label="strict_capacity_logit",
                     y_key="capacity_logit",bootstrap=True)
    adequate=[r for r in strong if r["adequate"]]
    quality=adjusted(adequate,full_families,label="strict_host_adequate_only",bootstrap=True)
    quality_capacity=adjusted(adequate,full_families,label="strict_capacity_adequate_only",
                              y_key="capacity_logit",bootstrap=True)
    capacity_edge=adjusted(ceiling_rows,full_families,label="strict_capacity_plus_edge",
                           y_key="capacity_logit",bootstrap=False,edge=True)
    broad=annotate(all_rows,broad=True)
    broadmod=adjusted(broad,full_families,label="broad_weak_plus_strong",bootstrap=True)
    broad_capacity=adjusted([r for r in broad if r["capacity_logit"] is not None],
                        full_families,label="broad_capacity_logit",
                        y_key="capacity_logit",bootstrap=False)
    # Disturbance groups differ in their taxonomic sampling; genus-level contrasts
    # are important, but do not become independent tests of global treatment.
    mixed=within_mixed_genus(strong)
    mixed_capacity=within_mixed_genus(ceiling_rows,y_key="capacity_logit")
    fam=taxonomy_family_signs(strong)
    cv=grouped_out_of_genus_cv(strong,full_families)
    cv_capacity=grouped_out_of_genus_cv(ceiling_rows,full_families,y_key="capacity_logit")
    genus_sizes=Counter(r["genus"] for r in strong)
    largest=sorted(genus_sizes.items(),key=lambda p:(-p[1],p[0]))[0]
    no_largest=[r for r in strong if r["genus"] != largest[0]]
    leave_largest=adjusted(no_largest,full_families,label="strict_excluding_largest_genus",
                            bootstrap=False)
    leave_largest_capacity=adjusted(no_largest,full_families,
                        label="strict_capacity_excluding_largest_genus",
                        y_key="capacity_logit",bootstrap=False)
    outcome={
        "schema":"chocho_global_butterfly_disturbance_macrofilter_v0.1",
        "status":"EXPLORATORY_MACRO_TRAIT_ASSOCIATION_NOT_OBSERVED_BUTTERFLY_DISTRIBUTION",
        "source_digests":digests,
        "source_traits":"LepTraits-based 339 species descriptor; exact taxon matched",
        "units":"BUTTERFLY_SPECIES, not observed populations or adult abundance",
        "scope":{"resource_species":239,"trait_candidates":339,
                 "region_grid_WGSRPD3":N_REGIONS,
                 "butterfly_families":len(full_families),
                 "resource_original_total_native":26530,
                 "resource_contemporary_total":41083,
                 "resource_introduced_added_total":14553},
        "population_selection":selection_audit(all_rows,trait_source),
        "strict_groups":{
            "exposure_class":STRONG_ASSOC,
            "comparison_class":STRONG_AVOID,
            "n_total":len(strong),
            "affiliated":describe([r for r in strong if r["exposure"]==1]),
            "avoidant":describe([r for r in strong if r["exposure"]==0]),
            "unadjusted_mean_log_difference":unadjusted_mean_diff,
            "unadjusted_effect_ratio_of_C_over_N_ratios":math.exp(unadjusted_mean_diff)
        },
        "strict_primary_adjusted":main,
        "strict_plus_edge_adjusted":edge,
        "strict_capacity_filled_logit":ceiling,
        "strict_host_adequate_only":quality,
        "strict_capacity_adequate_only":quality_capacity,
        "strict_capacity_plus_edge":capacity_edge,
        "broad_weak_plus_strong":broadmod,
        "broad_capacity_logit":broad_capacity,
        "within_genus_descriptive":mixed,
        "within_genus_capacity_descriptive":mixed_capacity,
        "taxonomy_family_stratification":fam,
        "leave_largest_genus": {
            "largest_genus":largest[0], "excluded_species":largest[1],
            "adjusted":leave_largest,
            "capacity_logit_adjusted":leave_largest_capacity
        },
        "out_of_genus_cv":cv,
        "out_of_genus_cv_capacity_logit":cv_capacity,
        "negative_gates":{
            "population_recruitment_measured":False,
            "actual_butterfly_occupancy_change_measured":False,
            "historical_plant_introduction_time_series_measured":False,
            "geographic_species_native_continent_adjusted":False,
            "habitat_trait_classification_experimentally_validated":False,
            "subnational_observer_effort_adjusted":False,
            "qualitative_habitat_descriptors_may_be_ascertainment_biased":True,
            "test_was_designed_after_descriptive_data_review":True,
            "not_confirmatory":True
        }
    }
    return outcome


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-json",type=Path,required=True)
    args=parser.parse_args()
    ans=main_analysis()
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(ans,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":ans["status"],
        "population_selection":ans["population_selection"]["by_disturbance"],
        "strong_unadjusted_log_contrast":ans["strict_groups"]["unadjusted_mean_log_difference"],
        "strong_adjusted_beta":ans["strict_primary_adjusted"]["adjusted_coefficient"],
        "strong_genus_ci":ans["strict_primary_adjusted"]["approx_genus_cluster_95ci"],
        "strong_boot_ci":ans["strict_primary_adjusted"]["genus_bootstrap"]["percentile_95ci"],
        "strong_wild_tail":ans["strict_primary_adjusted"]["wild_reduced_genus_signflip"]["two_sided_tail_fraction_not_randomization_p"],
        "strict_adequate_beta":ans["strict_host_adequate_only"]["adjusted_coefficient"],
        "strict_adequate_capacity_beta":ans["strict_capacity_adequate_only"]["adjusted_coefficient"],
        "strict_ceiling_logit_beta":ans["strict_capacity_filled_logit"]["adjusted_coefficient"],
        "strict_edge_adjusted_beta":ans["strict_plus_edge_adjusted"]["adjusted_coefficient"],
        "broad_beta":ans["broad_weak_plus_strong"]["adjusted_coefficient"],
        "full_sensitivity_effect_sizes":{
            k:{
                "n":ans[k]["n"],
                "beta":ans[k]["adjusted_coefficient"],
                "genus_cluster_ci":ans[k]["approx_genus_cluster_95ci"],
                "genus_bootstrap_ci":ans[k].get("genus_bootstrap",{}).get("percentile_95ci"),
                "wild_null_tail":ans[k]["wild_reduced_genus_signflip"].get(
                    "two_sided_tail_fraction_not_randomization_p")
            } for k in (
                "strict_primary_adjusted","strict_plus_edge_adjusted",
                "strict_capacity_filled_logit","strict_host_adequate_only",
                "strict_capacity_adequate_only","strict_capacity_plus_edge",
                "broad_weak_plus_strong","broad_capacity_logit")
        },
        "taxonomy_family_stratification":ans["taxonomy_family_stratification"],
        "leave_largest_genus":ans["leave_largest_genus"],
        "mixed_genus":ans["within_genus_descriptive"],
        "mixed_genus_capacity":ans["within_genus_capacity_descriptive"],
        "out_of_genus_cv": {
            k:v for k,v in ans["out_of_genus_cv"].items() if k!="individual_repeats"
        },
        "out_of_genus_cv_capacity_logit":{
            k:v for k,v in ans["out_of_genus_cv_capacity_logit"].items() if k!="individual_repeats"
        }
    },sort_keys=True))

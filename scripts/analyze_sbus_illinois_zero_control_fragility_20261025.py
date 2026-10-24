#!/usr/bin/env python3
"""Post-result challenge of Illinois host-sharing lag association.

Two fixed stress checks:
(1) remove mechanically composite site-year peer mean,
(2) keep only focal adult counts explicitly positive in both source years.
Never call the pilot coefficient causal competition.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from analyze_sbus_illinois_niche_counts_20261025 import (
    PROGRAMS, calc_panel, frozen_focal, group_within_matrix,
    host_sets, load_source, read_events, regression
)
ROOT = Path(__file__).resolve().parents[1]
EXPECTED_PRIMARY = 9178
EXPECTED_CONTRAST = -0.28350714681067035


def audit_model(obs, remove_site_mean=False):
    X,y,groups = group_within_matrix(obs)
    if remove_site_mean:
        # Controls columns: self, sharing peers, nonsharing peers, effort diff,
        # site-year community mean, year dummies. Drop only the composite.
        X=np.delete(X, 4, axis=1)
    n,k=X.shape
    sites=sorted({r["site"] for r in obs})
    if n < 300 or len(sites) < 50 or n<=k+100 or np.linalg.matrix_rank(X)!=k:
        return {"status":"INELIGIBLE_SAMPLE_SIZE_OR_RANK",
                "n_transitions":n,"n_sites":len(sites),"rank":int(np.linalg.matrix_rank(X)),
                "covariates":k}
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@beta
    siteidx=defaultdict(list)
    for i,r in enumerate(obs):
        siteidx[r["site"]].append(i)
    bread=np.linalg.pinv(X.T@X)
    meat=np.zeros((k,k))
    for ix in siteidx.values():
        z=X[ix].T @ resid[ix]
        meat+=np.outer(z,z)
    correction=len(siteidx)/(len(siteidx)-1)*(n-1)/(n-k)
    V=correction*bread@meat@bread
    stderr=np.sqrt(np.maximum(0,np.diag(V)))
    contrast=float(beta[1]-beta[2])
    contrastSE=math.sqrt(max(0,V[1,1]+V[2,2]-2*V[1,2]))
    return {
        "status":"OBSERVATIONAL_ASSOCIATION_ONLY",
        "n_transitions":n,
        "n_sites":len(siteidx),
        "n_species_site_clusters":len(groups),
        "drop_composite_site_year_mean":remove_site_mean,
        "beta_lag_shared":float(beta[1]),
        "shared_sitecluster_95_normal":[float(beta[1]-1.96*stderr[1]),
                                        float(beta[1]+1.96*stderr[1])],
        "beta_lag_nonshared":float(beta[2]),
        "contrast_shared_minus_nonshared":contrast,
        "contrast_sitecluster_95_normal":[contrast-1.96*contrastSE, contrast+1.96*contrastSE],
        "not_causal_or_population_density":True,
    }


def run(hostcsv:Path):
    focal=frozen_focal()
    hosts=host_sets(hostcsv,focal)
    with tempfile.TemporaryDirectory() as d:
        binary=Path(d)/"original_us_butterflies.zip"
        load_source(binary)
        events,counts=read_events(binary,focal)
    report, obs=calc_panel(PROGRAMS[0],events,focal,hosts)
    if report["site_species_transitions_with_both_comparator_groups"] != EXPECTED_PRIMARY:
        raise ValueError("The fixed original Illinois source panel drifted")
    baseline=regression(obs)
    if not math.isclose(baseline["contrast_shared_minus_nonshared"],
                        EXPECTED_CONTRAST, rel_tol=0, abs_tol=1e-9):
        raise ValueError("Original Illinois observed contrast drifted")
    positive=[o for o in obs if o["focal_positive_both_years"]]
    summary={
        "schema":"chocho_us_sbus_illinois_host_sharing_posthoc_fragility_v01",
        "status":"POSTHOC_OBSERVATIONAL_SENSITIVITY_COMPLETE_NO_COMPETITION_CAUSAL_ESTIMATE",
        "source_doi":"10.6084/m9.figshare.27934602.v1",
        "original_observations":len(obs),
        "positive_focal_both_year_observations":len(positive),
        "focal_positive_both_year_fraction":len(positive)/len(obs),
        "original_source_implied_zeros":report["site_year_species_missing_imputed_zeros"],
        "original_source_explicit_species_site_year":report["site_year_species_observed_entries"],
        "original_within_model":{
            "beta_shared":baseline["beta_shared_host_peer_lag"],
            "shared_SE":baseline["SE_shared_sitecluster"],
            "beta_nonshared":baseline["beta_nonshared_peer_lag"],
            "contrast":baseline["contrast_shared_minus_nonshared"],
            "contrast_CI":baseline["contrast_normal_95_interval"],
        },
        "leave_out_site_mean_control":audit_model(obs,remove_site_mean=True),
        "positive_focal_in_both_years_original_control_set":audit_model(positive),
        "positive_focal_in_both_years_no_site_mean":audit_model(positive,remove_site_mean=True),
        "limitations":[
            "Focal positivity conditioning selects survivors with observation, and cannot restore undocumented survey zeros",
            "Peers still use imputed zeros, so no complete elimination of source nonreporting risk",
            "All regressions remain observational adult counts per recorded survey hour",
            "Shared host documented in global HOSTS literature, not locally observed larvae",
            "No fine-scale weather, plant biomass, abundance calibration or phylogenetic structure controls",
            "Network degree-preserving null required before even a source-conditional niche association claim",
        ],
    }
    return summary


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--accepted-host-csv",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    result=run(a.accepted_host_csv)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "sample_original":result["original_observations"],
        "positive_both":result["positive_focal_both_year_observations"],
        "variants":{k:result[k] for k in (
            "leave_out_site_mean_control",
            "positive_focal_in_both_years_original_control_set",
            "positive_focal_in_both_years_no_site_mean",
        )},
    },ensure_ascii=False,indent=2))

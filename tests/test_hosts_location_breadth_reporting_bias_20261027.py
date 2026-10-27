"""No-network source-independent controls for genus-bootstrap HOSTS documentation models."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_hosts_location_breadth_reporting_bias_20261027 import (
    regression, design, bootstrap_by_genus, logit_inverse, get_source_rows,
)


def test_logistic_fit_checks_target_sample_and_numerics():
    xx=np.array([
        [1.,-3.,-1.,0.],
        [1.,-2.,-1.,0.],
        [1.,-1.,-1.,0.],
        [1.,0.,0.,0.],
        [1.,1.,1.,0.],
        [1.,2.,1.,0.],
        [1.,3.,1.,0.],
        [1.,4.,2.,0.],
    ])
    y=np.array([0.,0.,0.,0.,1.,1.,1.,1.])
    beta=regression(y,xx)
    assert np.isfinite(beta).all()
    assert beta[1]>0.
    assert logit_inverse(xx[-1]@beta)>logit_inverse(xx[0]@beta)
    with pytest.raises(ValueError,match="few positives"):
        regression(np.zeros(len(y)),xx)


def fixture(n=60):
    r=[]
    for i in range(n):
        g=i//2
        r.append({
            "species":f"Genus{g} species{i}",
            "original_source_host_family_count":(i%5)+1,
            "original_source_resolved_host_species":(i%5)+2,
            "original_source_introduced_added_resource_units":(i%15)+1,
            "source_unique_accepted_host_location_triples":(i%9)+1,
            "source_mapped_introduced_triples":int(i%4==0),
            "source_mapped_native_triples":int(i%3==0),
            "source_WCVP_unlisted_triples":int(i%5==0),
        })
    return r


def test_design_has_exact_three_nuisance_or_interest_predictors():
    r=fixture()
    x,(mu,sd)=design(r)
    assert x.shape==(60,4)
    assert np.allclose(x[:,0],1.)
    assert np.allclose(x[:,1:].mean(axis=0),0.,atol=1e-12)
    assert np.allclose(x[:,1:].std(axis=0,ddof=1),1.,atol=1e-12)
    z,_=design(r,feature="species_richness")
    assert z.shape==x.shape
    with pytest.raises(ValueError,match="Unfrozen"):
        design(r,feature="host_stage")


def test_source_reporting_controls_use_same_genus_resamples_for_negative_control():
    r=fixture(80)
    out=bootstrap_by_genus(r,draws=49,seed=10)
    assert out["sample_size"]==80
    assert out["source_genus_clusters"]==40
    assert out["bootstrap_draws"]==49
    assert out["outcome_models"]["introduced_plant_region_host_report"][
        "source_positive_butterflies"]==20
    assert out["outcome_models"][
        "introduced_plant_region_host_report__minus__native_plant_region_host_report"
    ]["bootstrap"]["valid_draws"]>=44


def test_original_239_source_integrity_is_gate_not_imputed_positive():
    bad={"schema":"other","result":{"primary_and_sensitivity":{"0.95":{}}}}
    with pytest.raises(ValueError,match="Not exact WCVP"):
        get_source_rows(bad)

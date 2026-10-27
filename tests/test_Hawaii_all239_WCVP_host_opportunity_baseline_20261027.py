"""No-network original 239 Hawaii source saturation and artificial rich-host reference."""
import math,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_Hawaii_all239_WCVP_host_opportunity_baseline_20261027 import (
    stratum,exact_hypergeom,exact_bin_preserving_reference,SCHEMA
)


def test_original_host_list_bin_frozen_before_ecological_comparison():
    assert [stratum(x) for x in (1,5,6,20,21,50,51,104)]==[
        "1-5","1-5","6-20","6-20","21-50","21-50","51+","51+"]
    with pytest.raises(ValueError,match="nonempty"):
        stratum(0)


def test_source_conditional_hypergeom_margins():
    z=exact_hypergeom(4,2,2)
    assert math.isclose(z[0],1/6)
    assert math.isclose(z[1],4/6)
    assert math.isclose(z[2],1/6)
    with pytest.raises(ValueError,match="Impossible"):
        exact_hypergeom(4,5,1)


def test_positive_ten_coverage_vs_original_source_population_not_arrival_negative():
    species=[
        {"species":f"Butterfly{i}",
         "source_accepted_host_count":2 if i<10 else 15,
         "hawaii_WCVP_introduced_wild_host_count":1 if i%3 else 0,
        } for i in range(30)
    ]
    positive=[x for x in species if x["species"] in {"Butterfly1","Butterfly2",
                                                     "Butterfly11","Butterfly12"}]
    x=exact_bin_preserving_reference(species,positive)
    assert x["source_Zenodo_10_original_positive_aliens"]==4
    assert x["source_observed_10_with_Hawaii_wild_introduced_host"]==3
    assert x["source_all239_host_breadth_bin_matched_artificial_expected"]>=0
    assert len(x["original_source_bins"])==4
    assert x["source_artificial_95pct_interval"][0] <= \
           x["source_artificial_95pct_interval"][1]
    assert x["reference_is_noncausal_missing_insect_arrival_exposure"]


def test_no_support_positive_kept_with_missing_not_ecological_absence():
    species=[
        {"species":f"Butterfly{i}",
         "source_accepted_host_count":i+1,
         "hawaii_WCVP_introduced_wild_host_count":0,
        } for i in range(12)
    ]
    positive=[species[0],species[8]]
    result=exact_bin_preserving_reference(species,positive)
    assert result["source_observed_10_with_Hawaii_wild_introduced_host"]==0
    assert result["source_all239_host_breadth_bin_matched_artificial_expected"]==0
    assert result["source_artificial_95pct_interval"]==[0,0]

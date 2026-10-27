"""Offline synthetic tests for exact source-only native/introduced pair overlap."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_hosts_native_introduced_pair_overlap_null_20261027 import (
    exact_hypergeom, exact_overlap_reference, family_strata,
    compare_documentation,
)


def test_hypergeom_margins_exact_probability():
    d=exact_hypergeom(4,2,2)
    assert sorted(d)==[0,1,2]
    assert math.isclose(d[0],1/6)
    assert math.isclose(d[1],4/6)
    assert math.isclose(d[2],1/6)
    assert math.isclose(sum(d.values()),1.0)
    assert exact_hypergeom(5,0,2)=={0:1.}
    assert exact_hypergeom(1,1,1)=={1:1.}
    with pytest.raises(ValueError):
        exact_hypergeom(3,4,1)
    with pytest.raises(ValueError):
        exact_hypergeom(3,1,4)


def test_exact_convolution_and_artificial_tails():
    ps=[
        {"species":"GenusA a","U":4,"N":2,"T":2,"O":1},
        {"species":"GenusB b","U":4,"N":2,"T":2,"O":1},
    ]
    v=exact_overlap_reference(ps,2)
    assert math.isclose(v["expected_shared_pairs_source_reference"],2)
    assert math.isclose(v["null_variance"],2/3)
    assert v["source_reference_95pct_interval"]==[0,4]
    assert math.isclose(
        v["source_mathematical_left_tail_not_biological_p"],
        v["source_mathematical_right_tail_not_biological_p"])
    assert v["permutable_strata"]==2
    assert v["source_recording_information_gate"]== \
        "SOURCE_REFERENCE_TOO_SPARSE_FOR_INTERPRETATION"


def test_family_condition_cannot_invent_permutable_host_identity():
    U={"pA","pB","pC","pD"}
    native={"pA","pC"}
    introduced={"pA","pD"}
    fam={"pA":["A"],"pB":["A"],"pC":["B"],"pD":["B"]}
    strata=family_strata(U,fam,native,introduced)
    assert len(strata)==2
    assert sum(x["O"] for x in strata)==1
    z=exact_overlap_reference([{**s,"species":"F butterfly"} for s in strata],1)
    assert z["source_reference_95pct_interval"]==[0,2]
    assert math.isclose(z["expected_shared_pairs_source_reference"],1)
    with pytest.raises(ValueError,match="ambiguous"):
        family_strata(U,{**fam,"pA":["A","B"]},native,introduced)


def test_no_native_documentation_does_not_mean_new_ecological_host():
    sets={"S a":{
        "all":{"1","2","3"},
        "NATIVE_CONFIRMED":set(),
        "INTRODUCED_CONFIRMED":{"1","2"},
        "WCVP_UNLISTED":{"3"},
    }}
    x=compare_documentation(sets,{},"INTRODUCED_CONFIRMED")
    assert x["source_margins"]["source_target_documented_species_host_pairs"]==2
    assert x["source_margins"]["source_native_and_target_documented_pairs"]==0
    assert x["within_butterfly_uniform_host_list_null"][
        "expected_shared_pairs_source_reference"]==0
    assert x["within_butterfly_uniform_host_list_null"][
        "source_recording_information_gate"]=="SOURCE_REFERENCE_TOO_SPARSE_FOR_INTERPRETATION"
    assert "adopted" not in str(x["source_margins"])


def test_source_family_and_unlisted_controls_share_original_universe():
    sets={
        "A butterfly":{
            "all":{"p1","p2","p3","p4"},
            "NATIVE_CONFIRMED":{"p1","p3"},
            "INTRODUCED_CONFIRMED":{"p1","p2"},
            "WCVP_UNLISTED":{"p3","p4"},
        },
        "B butterfly":{
            "all":{"p5","p6"},
            "NATIVE_CONFIRMED":{"p5"},
            "INTRODUCED_CONFIRMED":{"p6"},
            "WCVP_UNLISTED":{"p5"},
        }
    }
    fam={"p1":["X"],"p2":["X"],"p3":["Y"],"p4":["Y"],
         "p5":["X"],"p6":["X"]}
    intro=compare_documentation(sets,fam,"INTRODUCED_CONFIRMED")
    unlist=compare_documentation(sets,fam,"WCVP_UNLISTED")
    assert intro["source_margins"]["source_native_and_target_documented_pairs"]==1
    assert unlist["source_margins"]["source_native_and_target_documented_pairs"]==2
    assert intro["source_margins"]["accepted_host_universe_species_pair_sum"]==6
    assert unlist["source_margins"]["accepted_host_universe_species_pair_sum"]==6
    assert intro["within_butterfly_uniform_host_list_null"][
        "expected_shared_pairs_source_reference"]==1.5
    assert not intro["within_butterfly_uniform_host_list_null"][
        "conditional_model_is_artificial_not_biological_randomization"] is False

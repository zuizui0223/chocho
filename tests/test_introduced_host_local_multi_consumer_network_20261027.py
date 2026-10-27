"""No-network source independent tests for introduced-host co-reporting networks."""
import sys
from pathlib import Path
import math
from collections import Counter

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_introduced_host_local_multi_consumer_network_20261027 import (
    row_bin, hypergeom, exact_fixed_margin_source_reference, make_strata,
    summarise_status, classify_nodes, source_locality_concentration,
)

def node(pid, loc, family, status, degree, raw=None, different_genera=False):
    r=raw if raw is not None else degree
    return {
        "plant_id":pid,"original_Location":loc,"plant_family":family,
        "source_plant_region_origin":status,"raw_original_HOSTS_rows":r,
        "report_intensity_bin":row_bin(r),
        "distinct_original_239_butterfly_consumers":degree,
        "documented_at_least_two_butterfly_species":degree>=2,
        "documented_at_least_two_butterfly_genera":different_genera,
    }


def test_fixed_source_host_intensity_categories():
    assert [row_bin(v) for v in (1,2,3,5,6,150)] == [
        "1","2","3-5","3-5","6+","6+"]
    with pytest.raises(ValueError,match="No indexed"):
        row_bin(0)


def test_exact_origin_shared_reference_probabilities():
    p=hypergeom(4,2,2)
    assert p.keys()=={0,1,2}
    assert math.isclose(p[0],1/6) and math.isclose(p[1],4/6)
    assert math.isclose(p[2],1/6)
    with pytest.raises(ValueError,match="Invalid"):
        hypergeom(2,3,1)
    with pytest.raises(ValueError,match="Invalid"):
        hypergeom(0,0,0)
    strata=[
        {"location":"A","family":"Fabaceae","N":4,"target":2,
         "shared":2,"target_shared":1},
        {"location":"B","family":"Poaceae","N":4,"target":2,
         "shared":2,"target_shared":1},
    ]
    v=exact_fixed_margin_source_reference(strata,"INTRODUCED_CONFIRMED")
    assert v["source_observed_target_shared_nodes"]==2
    assert math.isclose(v["source_exact_conditional_expected_target_shared_nodes"],2)
    assert math.isclose(v["source_exact_conditional_variance"],2/3)
    assert v["source_permutable_location_count"]==2
    assert v["information_gate"]=="INSUFFICIENT_MATCHED_ORIGIN_GEOGRAPHY_INFORMATION"
    assert v["source_central_95pct_reference_interval"]==[0,4]
    assert v["source_conditional_permutable_support"]==[0,4]


def test_source_single_class_places_do_not_create_spurious_matched_effects():
    n=[
        node("a","IslandA","Fabaceae","INTRODUCED_CONFIRMED",2),
        node("b","IslandA","Fabaceae","INTRODUCED_CONFIRMED",1),
        node("c","IslandB","Fabaceae","NATIVE_CONFIRMED",2),
        node("d","IslandB","Fabaceae","NATIVE_CONFIRMED",1),
    ]
    s=make_strata(n,"INTRODUCED_CONFIRMED")
    x=exact_fixed_margin_source_reference(s,"INTRODUCED_CONFIRMED")
    assert x["source_matched_location_family_strata"]==0
    assert x["source_matched_target_plant_nodes"]==0
    assert x["source_central_95pct_reference_interval"]==[0,0]
    assert x["information_gate"]=="INSUFFICIENT_MATCHED_ORIGIN_GEOGRAPHY_INFORMATION"


def test_family_and_original_raw_report_count_matched_controls():
    n=[
        node("a","Japan","Fabaceae","INTRODUCED_CONFIRMED",3,raw=4),
        node("b","Japan","Fabaceae","INTRODUCED_CONFIRMED",1,raw=1),
        node("c","Japan","Fabaceae","NATIVE_CONFIRMED",2,raw=3),
        node("d","Japan","Fabaceae","NATIVE_CONFIRMED",1,raw=1),
        node("e","Japan","Poaceae","INTRODUCED_CONFIRMED",2,raw=2),
        node("f","Japan","Poaceae","NATIVE_CONFIRMED",2,raw=4),
        node("g","Taiwan","Fabaceae","WCVP_UNLISTED",2,raw=3),
        node("h","Taiwan","Fabaceae","NATIVE_CONFIRMED",1,raw=4),
    ]
    first=make_strata(n,"INTRODUCED_CONFIRMED",report_bin=False)
    strict=make_strata(n,"INTRODUCED_CONFIRMED",report_bin=True)
    assert sum(0<s["target"]<s["N"] for s in first)==2
    assert sum(0<s["target"]<s["N"] for s in strict)==2
    primary=exact_fixed_margin_source_reference(first,"INTRODUCED_CONFIRMED")
    assert primary["source_matched_target_plant_nodes"]==3
    assert primary["source_observed_target_shared_nodes"]==2
    source_bin=exact_fixed_margin_source_reference(strict,"INTRODUCED_CONFIRMED")
    assert source_bin["source_matched_target_plant_nodes"]==2
    control=exact_fixed_margin_source_reference(
        make_strata(n,"WCVP_UNLISTED"),"WCVP_UNLISTED")
    assert control["source_matched_target_plant_nodes"]==1
    assert control["source_observed_target_shared_nodes"]==1


def test_wcvp_native_unlisted_status_by_region_not_guessed_from_country():
    data=[
      node("x","Hawaii","Fabaceae","ORIGIN_PENDING",2),
      node("y","Hawaii","Fabaceae","ORIGIN_PENDING",1),
      node("z","Hawaii","Fabaceae","ORIGIN_PENDING",1),
    ]
    atoms={"US":[(("HAW",),100)]}
    metadata={"US":{"country_area_km2":100,"guard_coverage":1,
                    "is_geo_safe":True,"area_identity_checks":"PASS"}}
    x,edgecounts=classify_nodes(data,.95,atoms,metadata,{"HAW"},
        native={"x":{"HAW"}},
        current={"x":{"HAW"},"y":{"HAW"}})
    assert [v["source_plant_region_origin"] for v in x]==[
        "NATIVE_CONFIRMED","INTRODUCED_CONFIRMED","WCVP_UNLISTED"]
    assert edgecounts=={"NATIVE_CONFIRMED":2,
                        "INTRODUCED_CONFIRMED":1,"WCVP_UNLISTED":1} or \
           sum(edgecounts.values())==4


def test_node_rate_is_distinct_from_index_edge_rate():
    n=[
      node("a","Taiwan","Fabaceae","INTRODUCED_CONFIRMED",5,raw=8,
           different_genera=True),
      node("b","Taiwan","Fabaceae","INTRODUCED_CONFIRMED",1),
      node("c","Taiwan","Poaceae","NATIVE_CONFIRMED",1),
      node("d","Taiwan","Poaceae","NATIVE_CONFIRMED",2,
           different_genera=False)
    ]
    x=summarise_status(n)
    introduced=x["INTRODUCED_CONFIRMED"]
    assert introduced["accepted_plant_Location_nodes"]==2
    assert introduced["documented_butterfly_Location_host_edges"]==6
    assert introduced["shared_plant_nodes_2plus_butterflies"]==1
    assert math.isclose(introduced["source_fraction_plant_nodes_shared_by_two_butterflies"],.5)
    assert math.isclose(introduced["source_fraction_edges_on_shared_plant_nodes"],5/6)
    assert introduced["cross_genus_shared_plant_nodes"]==1
    assert len(source_locality_concentration(n,"INTRODUCED_CONFIRMED"))==1

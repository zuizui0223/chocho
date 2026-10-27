"""Synthetic independent spatial controls for source-coherent host-resource point maps."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_butterfly_73_point_tdwg_resource_20261027 import (
    RegionLocator, classify_point, fetch_and_measure, group_report,
)


def feature(name, lo, hi, south=0, north=1):
    return {
        "type":"Feature","properties":{"LEVEL3_COD":name},
        "geometry":{"type":"Polygon","coordinates":[[
            [lo,south],[hi,south],[hi,north],[lo,north],[lo,south]
        ]]}
    }


def fixture():
    codes=["AAA","BBB","CCC","DDD"]
    features=[
        feature("AAA",0,1),
        feature("BBB",1,2),
        feature("CCC",0,1,1,2),
        feature("DDD",0,1,2,3)
    ]
    bysp={"Test butterfly":{
        "native":(1<<2),
        "added":(1<<0)|(1<<1),
        "masks":{.05:1<<0,.10:0}
    }}
    return RegionLocator(features,codes),bysp


def test_polygon_index_classes_not_country_labels():
    loc, bysp=fixture()
    assert classify_point(.5,.5,"Test butterfly",loc,bysp)== \
        "introduced_analog_resource_region"
    assert classify_point(.5,1.5,"Test butterfly",loc,bysp)== \
        "introduced_novel_resource_region"
    assert classify_point(1.5,.5,"Test butterfly",loc,bysp)== \
        "native_resource_region"
    assert classify_point(2.5,.5,"Test butterfly",loc,bysp)== \
        "no_documented_added_resource_region"
    assert classify_point(5,5,"Test butterfly",loc,bysp)== "unmapped_point"
    assert classify_point(.5,1.0,"Test butterfly",loc,bysp)== \
        "ambiguous_tdwg_region"
    assert classify_point(.5,.5,"Test butterfly",loc,bysp,q=.10)== \
        "introduced_novel_resource_region"


def test_polygons_missing_regions_do_not_fall_back_to_country():
    with pytest.raises(ValueError,match="incomplete"):
        RegionLocator([feature("AAA",0,1)],["AAA","BBB"])


class DummyReader:
    def __init__(self, docs):
        self.docs=docs
        self.calls=0
        self.retries=0
    def get(self,url):
        self.calls+=1
        if not self.docs:
            raise AssertionError("Unexpected GBIF query")
        return self.docs.pop(0),"2"*64


def record(i, lat=.5, lon=.5, uncertain=250):
    return {
        "key":i,"decimalLatitude":lat,"decimalLongitude":lon,
        "year":2020,"countryCode":"CU","speciesKey":1944574,
        "basisOfRecord":"HUMAN_OBSERVATION",
        "coordinateUncertaintyInMeters":uncertain
    }


def test_fixed_country_novel_class_does_not_imply_point_novel_class():
    loc, bysp=fixture()
    case={
        "species":"Test butterfly","country":"CU",
        "resource_class":"novel","gbif_accepted_species_key":1944574,
        "original_gbif_record_count":2,
        "added_host_country_area_fraction":.4,
        "analog_host_country_area_fraction":.1,
        "novel_host_country_area_fraction":.3
    }
    fetcher=DummyReader([{"count":2,"results":[record(1),record(2,lon=1.5)]}])
    x=fetch_and_measure(case,fetcher,100,loc,bysp)
    assert x["full_source_status"]=="COMPLETE_SOURCE_POINT_TO_WGSRPD3_MAPPING"
    assert x["coordinate_valid_index_records"]==2
    assert x["regional_climate_point_classes"]["0.05"]=={
        "introduced_analog_resource_region":1,
        "introduced_novel_resource_region":1}
    assert x["regional_climate_point_classes"]["0.1"]=={
        "introduced_novel_resource_region":2}
    assert x["coordinate_strict_le1000m_records"]==2
    assert x["source_count_delta"]==0
    a=group_report([x],q="0.05")
    assert a["coordinate_valid_index_records"]==2
    assert a["species_country_equal_case_mean_introduced_resource_overlap"]==1


def test_incomplete_source_not_turned_into_zero_overlap():
    loc,bysp=fixture()
    case={"species":"Test butterfly","country":"CU",
          "resource_class":"novel","gbif_accepted_species_key":1944574,
          "original_gbif_record_count":1}
    bad=fetch_and_measure(case,DummyReader([
        {"count":2,"results":[record(1)]}]),100,loc,bysp)
    assert bad["full_source_status"]=="INCOMPLETE_SOURCE_POINT_MAPPING"
    report=group_report([bad])
    assert report["cases_source_complete"]==0
    assert report["coordinate_valid_index_records"]==0
    assert report["point_record_weighted_within_any_introduced_host_region"] is None

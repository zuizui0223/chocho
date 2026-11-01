"""No-network source tests for full 14.6M-record GBIF 369-region aggregation.

Do not downsample the real source silently or infer species absences from zeros.
Only synthetic points are used in these tests; GBIF raw coordinates are never
put in repository fixtures.
"""
from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import sys
import zipfile

import numpy as np
import pytest


ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"scripts/aggregate_gbif_global224_wgsrpd3_annual.py"


def module():
    pytest.importorskip("pandas")
    pytest.importorskip("shapely")
    spec=importlib.util.spec_from_file_location("global224_annual",SRC)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def toy_shapes():
    from shapely.geometry import box
    from shapely.strtree import STRtree
    return STRtree([box(0,0,2,2),box(2,0,4,2)]),["AAA","BBB"]


def test_vectorized_point_source_assignments_are_deterministic():
    m=module()
    tree,codes=toy_shapes()
    actual,overlaps=m.assign_xy_to_369(
        np.array([-1,1,2,3,3,0],float),np.array([1,1,1,1,92,0],float),
        tree)
    assert codes==["AAA","BBB"]
    assert actual.tolist()==[-1,0,0,1,-1,0]
    assert overlaps==1
    with pytest.raises(ValueError,match="dimensions"):
        m.assign_xy_to_369([0,1],[0],tree)


def toy_member():
    f=io.StringIO()
    w=csv.writer(f,delimiter="\t",lineterminator="\n")
    w.writerow(["gbifID","species","speciesKey","datasetKey","eventDate",
                "year","decimalLatitude","decimalLongitude",
                "coordinateUncertaintyInMeters","basisOfRecord",
                "occurrenceStatus"])
    w.writerow([1,"Genus alpha","101","d1","2020-01-01T00:00:00",2020,1,1,100,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([2,"Genus alpha","101","d1","2020-01-01T10:00:00",2020,1,1,200,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([3,"Genus beta","102","d2","2020-01-02",2020,1,3,1200,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([4,"Genus alpha","101","d2","2021-05-05",2021,1,3,800,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([5,"Genus alpha","101","d1","2019-01-01",2019,1,-5,100,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([6,"Genus alpha","101","d1","2020-01-01",2020,99,1,100,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([7,"Other species","999","d1","2020-01-01",2020,1,1,100,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([8,"Genus alpha","101","d1","2020-01-01",2020,1,1,100,
                "HUMAN_OBSERVATION","ABSENT"])
    w.writerow([9,"Genus beta","999","d2","2020-01-04",2020,1,3,100,
                "HUMAN_OBSERVATION","PRESENT"])
    w.writerow([10,"Genus beta","102","d2","2016-01-01",2020,1,3,100,
                "HUMAN_OBSERVATION","PRESENT"])
    return f.getvalue()


def test_source_csv_uses_species_keys_and_provider_date_unique_not_record_count(tmp_path):
    m=module()
    path=tmp_path/"toy.zip"
    with zipfile.ZipFile(path,"w",zipfile.ZIP_DEFLATED) as z:
        z.writestr("metadata.xml","test")
        z.writestr("occurrence.txt",toy_member()+"\n"*2500)
    with zipfile.ZipFile(path) as z:
        (member,sep,columns),info=m.archive_table(z)
        assert member=="occurrence.txt"
        assert sep=="\t"
        assert "specieskey" in columns
        with z.open(member) as stream:
            state=m.process_chunks(
                stream,sep,columns,*toy_shapes(),
                {"101":"Genus alpha","102":"Genus beta"},
                {"genus alpha":"Genus alpha","genus beta":"Genus beta"},
                1000,max_source_rows=10)
    q=state["source_counts"]
    assert q["rows"]==10
    assert q["candidate_records_mapped_to_369"]==6
    assert q["candidate_name_only_match"]==1
    assert q["not_a_frozen_candidate"]==1
    assert q["candidate_key_and_name_disagree"]==0
    assert q["source_occurrence_status_nonpresent"]==1
    assert q["candidate_valid_rows_outside_official_369_regions"]==1
    assert q["mapped_records_without_valid_dataset_date"]==1
    assert q["ambiguous_multiple_level3_polygon_hits"]==0
    assert state["record_cells"][("Genus alpha","AAA",2020)]==2
    assert len(state["provider_days"][("Genus alpha","AAA",2020)])==1
    assert state["record_cells"][("Genus beta","BBB",2020)]==3
    assert len(state["provider_days"][("Genus beta","BBB",2020)])==2
    assert state["precision_cells"][("Genus beta","BBB",2020)]==2
    assert m.accepted_iso_date("2022-01-01",2020,{}) is None
    assert m.accepted_iso_date("2020-01-02",2020,{})=="2020-01-02"


def test_source_panel_canonical_frozen_list_and_protocol():
    m=module()
    bykey,byname,resource,disposition=m.candidate_panel()
    assert len(bykey)==len(byname)==225
    assert len(resource)==239
    assert sum(disposition.values())==225
    assert m.EXPECTED_ZIP_BYTES==1415705617
    assert m.EXPECTED_RECORDS==14637152
    assert m.DOI=="10.15468/dl.umqmsf"

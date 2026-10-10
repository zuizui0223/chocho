"""Synthetic audit tests: 4 endpoint plants ≠ 4 timepoints."""
from __future__ import annotations
import csv
import importlib.util
import io
from pathlib import Path

PATH=Path(__file__).resolve().parents[1]/"scripts/audit_kyoto_competition_identifiability.py"
spec=importlib.util.spec_from_file_location("kyoto_design",PATH)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def build():
    l=io.StringIO()
    p=io.StringIO()
    lr=csv.DictWriter(l,fieldnames=["plot.no","treat.no","s.density","a.density","date","day","total.s","total.a","cumul.sp","cumul.ap"])
    pr=csv.DictWriter(p,fieldnames=["plot.no","treat.no","s.density","a.density","defoliation.1","defoliation.2","defoliation.3","defoliation.4"])
    lr.writeheader()
    pr.writeheader()
    cage=0
    for total in (4,8,12):
        for numerator in (0,1,2,3,4):
            if numerator==0:sr,ar=0,total
            elif numerator==4:sr,ar=total,0
            else:sr,ar=total*numerator//4,total-total*numerator//4
            for repeat in range(2):
                cage+=1
                pr.writerow({"plot.no":cage,"treat.no":cage,"s.density":sr,"a.density":ar,
                             "defoliation.1":10,"defoliation.2":20,"defoliation.3":30,"defoliation.4":40})
                for day in range(23):
                    lr.writerow({"plot.no":cage,"treat.no":cage,"s.density":sr,"a.density":ar,
                                 "date":"2018-07-01","day":day,
                                 "total.s":sr,"total.a":ar,"cumul.sp":0,"cumul.ap":0})
    return l.getvalue().encode(),p.getvalue().encode()

def test_four_plant_columns_cannot_be_time_series():
    larval,plants=build()
    result=mod.audit(larval,plants)
    assert result["n_larval_longitudinal_rows"]==690
    assert result["n_plant_final_rows"]==30
    assert result["n_independent_cages"]==30
    assert result["replicates_per_density_combination"]==2
    assert result["plant_variables"]["n_temporal_plant_biomass_measurements"]==0
    assert result["plant_variables"]["independently_manipulated_host_regrowth_timing"] is False
    assert result["outcome_variables"]["competitor_contact_or_behavior_measured"] is False

def test_exact_density_replication_is_required():
    larval,plants=build()
    text=plants.decode().replace("1,1,0,4,10,20,30,40","1,1,4,0,10,20,30,40")
    try:
        mod.audit(larval,text.encode())
    except ValueError as err:
        assert "treatment disagrees" in str(err)
    else:
        raise AssertionError("plant data swapped treatment not rejected")

def test_duplicate_day_cannot_create_extra_replicate():
    larval,plants=build()
    rows=larval.decode().splitlines()
    first=rows[1].split(",")
    second=rows[2].split(",")
    second[5]=first[5] # modify the exact day column, never string-replace treatment codes
    rows[2]=",".join(second)
    try:
        mod.audit(("\n".join(rows)+"\n").encode(),plants)
    except ValueError as err:
        assert "duplicated day" in str(err)
    else:
        raise AssertionError("duplicated day treated as new timepoint")

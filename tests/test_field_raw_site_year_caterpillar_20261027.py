"""Offline raw standardized field survey ID/site-year normalization guards."""
from __future__ import annotations
import csv,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_field_raw_site_year_caterpillar_20261027 import (
    raw_filter,site_year_matched,aggregate_matched, read_inferred, read_origin,
)


def test_raw_long_format_arthropods_not_independent_survey_events(tmp_path):
    header=["ID","Name","PlantFK","sciName","Species","julianday",
            "Longitude","Latitude","WetLeaves","Year","Group","Quantity",
            "ObservationMethod"]
    rows=[
        ["101","site A","Plant 1","Quercus alba","Quercus alba",
         "160","-82","35","FALSE","2022","caterpillar","3","Visual"],
        ["101","site A","Plant 1","Quercus alba","Quercus alba",
         "160","-82","35","FALSE","2022","spider","2","Visual"],
        ["102","site A","Plant 2","Rosa rugosa","Rosa rugosa",
         "160","-82","35","FALSE","2022","beetle","3","Visual"],
        ["103","site A","Plant 2","Rosa rugosa","Rosa rugosa",
         "160","-82","35","FALSE","2022","caterpillar","0","Visual"],
        ["104","Coweeta - BB","Plant 1","Quercus alba","Quercus alba",
         "160","-82","35","FALSE","2022","caterpillar","1","Visual"],
        ["105","site A","Plant 1","Quercus alba","Quercus alba",
         "160","-82","35","TRUE","2022","caterpillar","1","Visual"],
        ["106","site A","Plant 1","Quercus alba","Quercus alba",
         "300","-82","35","FALSE","2022","caterpillar","1","Visual"],
    ]
    f=tmp_path/"raw.csv"
    with f.open("w",newline="") as o:
        w=csv.writer(o);w.writerow(header);w.writerows(rows)
    out,qc,stats,columns=raw_filter(f,{},{
        "Quercus alba":"native","Rosa rugosa":"alien"
    },{"Quercus alba":2})
    assert stats["survey_IDs_final_quality_gate"]==3
    assert stats["source_caterpillar_positive_survey_IDs"]==1
    assert stats["accepted_239_linked_survey_IDs"]==1
    assert stats["USDA_native_survey_IDs"]==1
    assert stats["USDA_alien_survey_IDs"]==2
    assert len(out)==3
    assert qc["Coweeta_source_rows_excluded"]==1
    assert qc["wet_or_missing_leaf_condition"]==1
    assert qc["out_of_JuneJuly_eastNA_or_year"]==1
    assert out[0]["plant"]=="Quercus alba"


def test_both_origin_site_year_gate_no_false_competition_effect():
    rows=[]
    for site in ("S1","S2"):
        for year in (2021,2022):
            for i in range(8):
                rows.append({
                    "site":site,"year":year,"week":22,
                    "origin":"native" if i<4 else "alien",
                    "plant":"Quercus alba" if i<4 else "Rosa rugosa",
                    "taxon_status":"SPECIES_BINOMIAL",
                    "caterpillar":int(i%3==0),
                    "host_239":1 if i<4 else None,
                })
    strict=site_year_matched(rows,min_survey=5)
    assert strict==[]
    relaxed=site_year_matched(rows,min_survey=4)
    assert len(relaxed)==4
    totals=aggregate_matched(relaxed)
    assert totals["matched_site_year_or_week_strata"]==4
    assert totals["independent_field_sites_with_matching_both_plant_origins"]==2
    assert totals["source_site_year_paired_evidence_information_gate"]== \
        "INSUFFICIENT_MATCHED_SITES_OR_YEARS"
    assert "descriptive_Mantel_Haenszel_native_vs_alien_larval_presence_odds_ratio" not in totals


def test_plant_inference_source_confidence_not_genus_imputation(tmp_path):
    a=tmp_path/"inferred.csv"
    a.write_text("PlantFK,InferredSciName,NameConfidence\n"
                 "p1,Rosa rugosa,3\np2,Quercus alba,1\n")
    inf=read_inferred(a)
    assert inf["p1"]==("Rosa rugosa",3.)
    assert inf["p2"]==("Quercus alba",1.)
    b=tmp_path/"origin.csv"
    b.write_text("scientificName,nativeStatus,plantOrigin\n"
                 "Quercus alba,L48(N),native\n"
                 "Rosa rugosa,L48(I),alien\n")
    origins=read_origin(b)
    assert origins=={"Quercus alba":"native","Rosa rugosa":"alien"}


def test_matched_source_239_filter_is_not_populated_from_no_documented_host_list():
    rows=[]
    for i in range(10):
        rows.append({
            "site":"A","year":2020,"week":25,
            "origin":"alien" if i>=5 else "native",
            "plant":"Rosa rugosa" if i>=5 else "Quercus alba",
            "taxon_status":"SPECIES_BINOMIAL",
            "caterpillar":int(i%2==0),
            "host_239":None if i>=5 else 1
        })
    assert len(site_year_matched(rows,min_survey=5))==1
    assert not site_year_matched(rows,min_survey=5,source_links=True)
    assert len(site_year_matched(rows,min_survey=2,week=True))==1

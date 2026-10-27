"""Offline controls for accepted butterfly host links vs published woody field larvae."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_field_caterpillar_host_realization_filter_20261027 as m


def field_origin_fixture():
    rows=[]
    for i in range(100):
        genus="Zanthoxylum" if i==85 else "Genus"+chr(65+i//3)+chr(65+i%3)
        name=genus if i==85 else genus+" speciosa"
        origin="alien" if i>=86 else "native"
        n=40
        # Avoid nonexact source percentage quantization.
        pct=5 if origin=="native" else 2.5
        rows.append({
            "Plant species":name,
            "Family":"Fabaceae",
            "USDA Native Status":"L48(I,N)CAN(N)" if i==85 else "L48(N)",
            "Origin":origin,
            "n_Surveys":str(n),"n_Branches":"3","n_Sites":"2",
            "Pct with caterpillars":str(pct),
        })
    return rows


def test_source_binomial_is_not_species_level_when_only_genus():
    assert m.normalize_binomial("Quercus alba")=="Quercus alba"
    assert m.normalize_binomial("  Quercus   alba ")=="Quercus alba"
    assert m.normalize_binomial("Zanthoxylum") is None
    assert m.normalize_binomial("Quercus sp.") is None
    assert m.normalize_binomial("Quercus spp.") is None
    assert m.normalize_binomial("Quercus alba extra") is None


def test_frozen_table_100_plants_including_one_genus_but_no_origin_missing(tmp_path,monkeypatch):
    data=field_origin_fixture()
    f=tmp_path/"field.csv"
    with f.open("w",newline="",encoding="utf-8") as stream:
        out=csv.DictWriter(stream,fieldnames=list(data[0]))
        out.writeheader();out.writerows(data)
    monkeypatch.setattr(m,"EXTERNAL_BLOB",m.git_blob_sha(f.read_bytes()))
    plants,meta=m.original_survey_plant_panel(f)
    assert len(plants)==100
    assert meta["original_100_plants_by_USDA_origin"]=={"native":86,"alien":14}
    assert meta["original100_genus_only_taxa"]==1
    assert meta["original100_species_level_rows"]==99
    name=next(x for x in plants if x["field_plant_binomial"]=="Zanthoxylum")
    assert name["original_USDA_origin"]=="native"
    assert name["field_plant_taxonomic_resolution"]=="GENUS_ONLY_UNSCORABLE"
    out,status=m.join(plants,{"GenusAA speciosa":{"p1"}},{"p1":{"Butterfly x"}},set(),{})
    assert status["EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH"]==1
    assert status["GENUS_ONLY_UNSCORABLE"]==1
    assert any(x["original239_accepted_HOSTS_match_status"]==
               "GENUS_ONLY_UNSCORABLE" for x in out)
    assert m.model_records(out,"original239",include_zero=True)[0].__len__()==99


def test_source_taxonomic_ambiguous_accepted_plant_id_does_not_invent_single_host():
    rows=[
        {"field_plant_binomial":"Acer palmatum",
         "field_plant_taxonomic_resolution":"SPECIES_BINOMIAL"},
        {"field_plant_binomial":"Quercus alba",
         "field_plant_taxonomic_resolution":"SPECIES_BINOMIAL"}
    ]
    for r in rows:
        r.update({
            "original_USDA_origin":"alien",
            "field_plant_genus":r["field_plant_binomial"].split()[0],
            "field_original_caterpillar_pct":0.,
            "field_original_n_surveys":50,
            "field_reconstructed_positive_survey_events":0,
        })
    out,qc=m.join(rows,{
        "Acer palmatum":{"p1","p2"},"Quercus alba":{"p3"}
    },{
        "p1":{"a b"},"p2":{"c d"},"p3":{"e f"}
    },{"Acer palmatum"},{})
    assert qc["AMBIGUOUS_SOURCE_ACCEPTED_PLANT_MATCH"]==1
    assert out[0]["documented_original_239_butterfly_host_species_count"] is None
    assert qc["EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH"]==1
    assert out[1]["documented_original_239_butterfly_host_species_count"]==1


def test_original_surveys_are_not_independent_species_and_binomial_n_rejected():
    x=np.column_stack([np.ones(10),np.arange(10)/10,np.arange(10)%2])
    y=np.array([1,0,1,0,2,1,2,1,3,2],float)
    n=np.array([10]*10,float)
    result=m.solve_logit(x,y,n)
    assert np.all(np.isfinite(result))
    assert m.expit(x@result).min()>0
    with pytest.raises(ValueError,match="Invalid original rounded"):
        m.solve_logit(x,[11]+list(y[1:]),n)
    with pytest.raises(ValueError,match="rank insufficient"):
        m.solve_logit(np.ones((10,3)),y,n)


def test_source_reporting_model_declines_poor_exact_host_taxonomic_overlap():
    data=[]
    for i in range(20):
        r={
            "field_plant_binomial":f"Plant{i} testii",
            "field_plant_genus":f"Plant{i}",
            "field_plant_taxonomic_resolution":"SPECIES_BINOMIAL",
            "original_USDA_origin":"alien" if i<2 else "native",
            "original239_accepted_HOSTS_match_status":"EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH",
            "documented_original_239_butterfly_host_species_count":i+1,
            "source_broader_original_HOSTS_lepidopteran_consumers":i+1,
            "field_original_n_surveys":100,
            "field_reconstructed_positive_survey_events":2 if i<2 else 8,
            "field_original_caterpillar_pct":2. if i<2 else 8.
        }
        data.append(r)
    r=m.conditional_model(data,"original239",boot=9)
    assert r["prespecified_information_gate"]=="INSUFFICIENT_EXACT_MATCH_NATIVE_ALIEN_SUPPORT"
    assert "source_host_count_per_SD_adjusted_odds_ratio" not in r
    assert r["alien_plants"]==2


def test_equal_plants_vs_survey_weights_and_plant_genus_out_of_sample_validation():
    data=[]
    for i in range(36):
        alien=i>=24
        orig=(1+(i%6))*1.0
        rate=(.02 if alien else .08)+(i%4)*.004
        n=100
        data.append({
            "field_plant_binomial":f"Plant{chr(65+i)} veranii",
            "field_plant_genus":f"Plant{chr(65+i)}",
            "field_plant_taxonomic_resolution":"SPECIES_BINOMIAL",
            "original_USDA_origin":"alien" if alien else "native",
            "original239_accepted_HOSTS_match_status":"EXACT_SOURCE_239_ACCEPTED_PLANT_MATCH",
            "documented_original_239_butterfly_host_species_count":orig,
            "source_broader_original_HOSTS_lepidopteran_consumers":orig*2,
            "field_original_n_surveys":n,
            "field_reconstructed_positive_survey_events":round(n*rate),
            "field_original_caterpillar_pct":100*rate,
        })
    z=m.conditional_model(data,"original239",boot=49)
    assert z["prespecified_information_gate"]=="SOURCE_LINKED_FIELD_MODEL_INFORMATION_OK"
    assert z["leave_one_plant_genus_out"]["validated_heldout_genera"]>=25
    assert z["plant_genus_bootstrap_499"]["successful_draws"]>=44
    assert math.isfinite(z["source_host_count_per_SD_adjusted_odds_ratio"])
    assert z["leave_one_plant_genus_out"][
        "incremental_log_predictive_gain_per_1000_surveys"] is not None

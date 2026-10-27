"""Source-independent site-year host portfolio incremental prediction controls."""
from __future__ import annotations
import math,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_field_site_year_HOSTS_portfolio_predictiveness_20261027 import (
    aggregate_groups,assess_source_field_link_counts,
    build_X,fit_site_year_model,blocked_leave_site_out
)


def synthetic_source_surveys():
    result=[]
    plants=[
        (f"Nativeplant{chr(65+i)} alba","native",1+i%4)
        for i in range(6)
    ]+[
        (f"Alienplant{chr(65+i)} verna","alien",1+i%5)
        for i in range(6)
    ]
    for site in range(12):
        for year in (2020,2021):
            for plant,origin,consumer in plants:
                for i in range(6):
                    # Any caterpillar on branch all-Lepidoptera presence.
                    value=int((i+site+consumer)%11 < (
                         2 if origin=="native" else 1))
                    result.append({
                        "site":f"Original Field Site {site}",
                        "year":year,"week":25,
                        "plant":plant,"origin":origin,
                        "host_239":consumer,
                        "caterpillar":value,
                        "taxon_status":"SPECIES_BINOMIAL"
                    })
    return result


def test_matched_original_source_sites_retained_only_when_siteyear_pair_support():
    raw=synthetic_source_surveys()
    matched=[
        {"site":f"Original Field Site {s}",
         "site_year_key":(f"Original Field Site {s}",year)}
        for s in range(12) for year in (2020,2021)
    ]
    grouped=aggregate_groups(raw,matched)
    assert len(grouped)==12*2*12
    qc=assess_source_field_link_counts(grouped)
    assert qc["source_site_year_strata"]==24
    assert qc["source_field_sites"]==12
    assert qc["original_239_linked_field_native_plant_species"]==6
    assert qc["original_239_linked_field_alien_plant_species"]==6
    assert qc["model_information_gate"]=="INSUFFICIENT_MATCHED_BUTTERFLY_HOST_SOURCE_FIELD_SUPPORT"
    # Only 12 plants / genera <15, so numerical model is intentionally not
    # promoted even with many repeated branch surveys.
    out=fit_site_year_model(grouped)
    assert out["status"]=="NO_MODEL_INSUFFICIENT_MATCHED_SOURCE_SUPPORT"


def test_model_can_have_multiple_plants_and_source_holdout_sites():
    raw=synthetic_source_surveys()
    # 4 further native and 4 alien woody plant species to clear 15-plant gate.
    for site in range(12):
        for year in (2020,2021):
            for i in range(8):
                for j in range(6):
                    orig="native" if i<4 else "alien"
                    sp=f"Extra{chr(65+i)} corymbosa"
                    raw.append({
                        "site":f"Original Field Site {site}","year":year,"week":25,
                        "plant":sp,"origin":orig,"host_239":1+(i%3),
                        "caterpillar":int((site+i+j)%12<2 if orig=="native"
                                          else (site+i+j)%12<1),
                        "taxon_status":"SPECIES_BINOMIAL"
                    })
    matched=[
        {"site":f"Original Field Site {s}",
         "site_year_key":(f"Original Field Site {s}",year)}
        for s in range(12) for year in (2020,2021)
    ]
    grouped=aggregate_groups(raw,matched)
    q=assess_source_field_link_counts(grouped)
    assert q["model_information_gate"]=="SOURCE_LINKED_PLANT_AND_SITE_YEAR_SUPPORT"
    assert q["original_239_HOSTS_exact_linked_field_plants"]==20
    X,y,n,names,standard=build_X(grouped,with_site_year=True)
    assert len(names)==24
    assert X.shape[1]==26
    assert len(y)==len(n)==len(grouped)
    assert np.linalg.matrix_rank(X)==X.shape[1]
    model=fit_site_year_model(grouped)
    assert model["status"]=="DESCRIPTIVE_SOURCE_SITE_YEAR_FIXED_EFFECT_MODEL_ONLY"
    assert math.isfinite(model["original_239_HOSTS_consumers_perSD_adjusted_odds_ratio"])
    assert model["leave_one_entire_survey_site_out"]["sites_validated"]==12
    assert math.isfinite(model["leave_one_entire_survey_site_out"][
        "incremental_log_score_per_1000_out_of_site_surveys"])


def test_unlinked_source_239_plants_never_force_false_zero_accepted_host_count():
    raw=synthetic_source_surveys()
    raw.append({
        "site":"Original Field Site 0","year":2020,"week":25,
        "plant":"Unknown plant","origin":"alien",
        "host_239":None,"caterpillar":1,
        "taxon_status":"GENUS_ONLY"
    })
    matched=[{"site":"Original Field Site 0",
              "site_year_key":("Original Field Site 0",2020)}]
    grouped=aggregate_groups(raw,matched)
    assert all(x["source_239_butterfly_consumers"] is not None for x in grouped)
    assert not any(x["plant"]=="Unknown plant" for x in grouped)

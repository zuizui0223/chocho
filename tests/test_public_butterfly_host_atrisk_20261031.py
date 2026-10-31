"""Regression gates for butterfly-first port/foodplant island index analysis.

The butterfly index is a current GBIF *indexing* cohort, not historical
immigration, occupancy or breeding. The 35 subset includes only six named
documented HOSTS plants and is prior-exploratory selected. A no-index indicator
can never be relabelled as foodplant biological absence.

Separate GEB anonymous bundle is intentionally redacted. Only skip if it
identifiably lacks these exploratory inputs; full source repo must replay.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import csv
import json
import math
import pytest

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/public_butterfly_island_risk_20261031"
NEEDED=["full17_query_receipt.json","full17_butterfly_island_225_atrisk.csv",
        "public_35host_named_plant_join.json",
        "public_35host_butterfly_x_17island_risk.csv"]


def sources():
    if "_anonymous_stage" in ROOT.parts or not all((D/p).exists() for p in NEEDED):
        pytest.skip("Independent blinded GEB publication excludes the separate exploratory island source datasets")
    source=json.loads((D/NEEDED[0]).read_text(encoding="utf8"))
    matched=json.loads((D/NEEDED[2]).read_text(encoding="utf8"))
    def table(name):
        with (D/name).open(encoding="utf8",newline="") as f:
            return list(csv.DictReader(f))
    return source,table(NEEDED[1]),matched,table(NEEDED[3])


def test_full225_by17_islands_no_future_based_species_selection():
    source,rows,matched,_=sources()
    assert source["status"]=="PUBLIC_REINDEXED_OCCURRENCE_HISTORY_NOT_TRUE_COLONIZATION_OR_HOST_MEDIATION"
    assert source["species_candidates_frozen_before_future_225"]==225
    assert source["all_17_original_islands"]==17
    assert len(source["islands_tested"])==17
    assert len(rows)==225*17==source["risk_cells"]==3825
    assert len({r["speciesKey"] for r in rows})==225
    assert len({r["landmass_id"] for r in rows})==17
    assert len({(r["speciesKey"],r["landmass_id"]) for r in rows})==3825
    assert all(r["source_not_future_selected"]=="True" for r in rows)
    assert all(r["index_not_true_absence_or_colonization"]=="True" for r in rows)
    first=sum(int(r["first_indexed_in_2011_2017"]) for r in rows)
    holdout=sum(int(r["first_indexed_in_2018_2022"]) for r in rows)
    eligible_first=sum(int(r["eligible_2011_2017"]) for r in rows)
    eligible_holdout=sum(int(r["eligible_2018_2022"]) for r in rows)
    assert (eligible_first,first,eligible_holdout,holdout)==(3584,80,3504,18)
    assert (first,holdout,eligible_first,eligible_holdout)==(
        source["first_indexed_2011_2017"],source["first_indexed_2018_2022"],
        source["eligible_2011_2017"],source["eligible_2018_2022"])
    assert all(int(r["eligible_2011_2017"]) == (int(r["pre2001_2010_indexed_records"])==0) for r in rows)
    assert all(int(r["eligible_2018_2022"]) == (
       int(r["pre2001_2010_indexed_records"])==0 and
       int(r["future2011_2017_indexed_records"])==0) for r in rows)
    assert all(int(r["first_indexed_in_2011_2017"]) == (
       int(r["eligible_2011_2017"])==1 and
       int(r["future2011_2017_indexed_records"])>0) for r in rows)
    assert all(int(r["first_indexed_in_2018_2022"]) == (
       int(r["eligible_2018_2022"])==1 and
       int(r["future2018_2022_indexed_records"])>0) for r in rows)
    assert all(r["index_not_true_absence_or_colonization"]=="True" for r in rows)


def test_original_island_area_and_butterfly_taxon_source_envelope():
    source,_,_,_=sources()
    filter_=source["butterfly_taxon_filter"]
    assert filter_["method"]=="REPEATED_GBIF_FAMILY_TAXONKEY_OR_FILTER"
    assert len(filter_["family_taxonkeys"])==6
    assert set(filter_["family_taxonkeys"])=={
        "Nymphalidae","Pieridae","Lycaenidae","Papilionidae","Hesperiidae","Riodinidae"}
    assert len(set(filter_["family_taxonkeys"].values()))==6
    assert source["gshhg_shape_sha256"]=="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
    geometry=source["geometry_audit"]
    assert len(geometry)==17
    assert sum(x["query_fragments"] for x in geometry)==24
    assert all(x["status"]=="AREA_PRESERVED_COUNTERCLOCKWISE_DISJOINT_FRAGMENT_GEOMETRY"
              for x in geometry)
    assert all(.97<=x["simplified_area_retained"]<=1.03 for x in geometry)
    assert all(x["max_fragment_WKT_chars"]<=2300 for x in geometry)
    assert next(x for x in geometry if x["island_id"]=="110")["query_fragments"]==8
    receipts=source["GBIF_index_query_receipts"]
    assert len(receipts)==17*3
    assert all(r["geometry_fragments"]==len(r["fragment_receipts"]) for r in receipts)
    assert all(r["all_butterfly_indexed_occurrences_disjoint_polygon_sum"]>=0 for r in receipts)
    assert any("GBIF" in s and "index" in s for s in source["limits"])


def test_42_confirmed_HOSTS_larval_links_joined_to_35_butterflies_without_new_host_invention():
    source,all225,result,join=sources()
    assert result["status"]=="EXPLORATORY_NAME_SPECIFIC_HISTORICAL_HOST_BUTTERFLY_INDEX_JOIN_NOT_CAUSAL"
    assert (result["original_225_predeclared_species"],result["host_species_butterflies"],
            result["source_TRUE_Host_edges"],result["foodplant_species_names"],
            result["candidate_species_island_cells"],result["islands"])==(225,35,42,6,595,17)
    assert len(join)==595
    assert len({(x["speciesKey"],x["landmass_id"]) for x in join})==595
    assert len({x["species"] for x in join})==35
    original_by={(x["speciesKey"],x["landmass_id"]):x for x in all225}
    for x in join:
        original=original_by[(x["speciesKey"],x["landmass_id"])]
        assert all(x[k]==original[k] for k in (
            "pre2001_2010_indexed_records",
            "future2011_2017_indexed_records","future2018_2022_indexed_records",
            "eligible_2011_2017","first_indexed_in_2011_2017",
            "eligible_2018_2022","first_indexed_in_2018_2022"))
        assert x["not_selected_35_on_future_first_index"]=="True"
        assert x["unrecorded_named_foodplant_not_biological_absence"]=="True"
        assert all(int(x[k]) in (0,1) for k in (
           "pre2011_port_compatible_speciesKey_recorded",
           "pre2011_offport_compatible_speciesKey_recorded",
           "pre2011_any_compatible_speciesKey_recorded"))
        assert int(x["pre2011_any_compatible_speciesKey_recorded"])==int(
           bool(int(x["pre2011_port_compatible_speciesKey_recorded"]) or
                int(x["pre2011_offport_compatible_speciesKey_recorded"])))
    assert len(result["taxonomy"])==6
    assert all(x["source_match_type"]=="EXACT" and x["source_rank"]=="SPECIES"
               and x["source_confidence"]==100 and
               x["name_applies_to_historic_island_alien_status"] is False for x in result["taxonomy"])


def test_historical_host_vs_future_index_rates_and_island_fixed_margins_replay():
    _,_,result,rows=sources()
    for period,expect in (("2011_2017",(505,33,20,0,37,2)),
                          ("2018_2022",(472,6,20,0,35,0))):
        eligible="eligible_"+period
        outcome="first_indexed_in_"+period
        source=[x for x in rows if x[eligible]=="1"]
        counts=result["periods"][period]
        assert (len(source),sum(int(x[outcome]) for x in source))==expect[:2]
        assert (counts["risk_cells"],counts["new_indexed_butterfly_island_pairs"])==expect[:2]
        for label,key,expo,positive in (
           ("port","pre2011_port_compatible_speciesKey_recorded",expect[2],expect[3]),
           ("off","pre2011_offport_compatible_speciesKey_recorded",expect[4],expect[5])):
            observed=[x for x in source if x[key]=="1"]
            assert len(observed)==expo
            assert sum(int(x[outcome]) for x in observed)==positive
            block=counts["by_early_host_presence"][key]
            assert block["1"]["risk_cells"]==expo
            assert block["1"]["first_indexed"]==positive
            assert block["0"]["risk_cells"]==expect[0]-expo
            assert block["0"]["first_indexed"]==expect[1]-positive
            by=defaultdict(list)
            for x in source:by[x["landmass_id"]].append(x)
            null=counts["island_margin_null_"+key]
            target=sum(
                sum(int(x[key]) for x in sub)
                *sum(int(x[outcome]) for x in sub)/len(sub)
                for sub in by.values())
            assert math.isclose(target,null["conditional_island_margin_expected"],
                                rel_tol=0,abs_tol=1e-11)
            assert null["observed_coindexed_new_records"]==positive
            assert null["uncertainty_is_not_corrected_for_repeated_species_or_recording_bias"]
    assert result["source_direct_HOSTS_butterflies_are_historically_posthoc_selected"] is True
    assert any("not botanical absence" in s for s in result["limits"])
    assert any("ship" in s for s in result["limits"])

"""SYNTHETIC source parser and unit-of-randomization tests; no real monarch data in fixtures."""
from __future__ import annotations

import csv
import importlib.util
import io
import math
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"scripts/reanalyze_original_monarch_plot_host_warming_adults.py"
SPEC=importlib.util.spec_from_file_location("monarch_plots",P)
assert SPEC and SPEC.loader
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def fake_rows(*, one_unknown=False, nonzero_main_effects=False):
    rows=[]
    for plot in range(1,31):
        temp="ambient" if plot<=15 else "elevated"
        for host in m.HOST:
            for oe in m.INOC:
                plant=f"P{plot}-{host}-{oe}"
                for ind in (1,2):
                    # Fully fabricated 2x2x2 values; randomized temperature
                    # at plot, NOT individual larva.
                    survival="1"
                    if host=="tropical" and temp=="elevated" and ind==1:
                        survival="0"
                    elif nonzero_main_effects and oe=="infected" and ind==2:
                        survival="0"
                    if one_unknown and plot==1 and host=="tropical" and oe=="control" and ind==2:
                        survival="NA"
                    rows.append({
                        "ID":f"{plant}-L{ind}","Plant_ID":plant,
                        "PlotNum":f"G{plot:02d}","Temp":temp,
                        "Milkweed":host,"OE_treatment":oe,"Lineage":"L1",
                        "Notes":"synthetic not source original",
                        "Surv_adult":survival,"Surv_pupa":survival,
                    })
    return rows


def encode(rows):
    f=io.StringIO()
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return f.getvalue().encode()


def load_fake(**kwargs):
    return m.load_original(encode(fake_rows(**kwargs)),require_blob=False)


def test_source_validates_true_30_plot_clusters_not_240_independent_temps():
    rows,receipt=load_fake()
    assert len(rows)==240
    assert receipt["n_original_larvae"]==240
    assert receipt["n_distinct_experimental_plots"]==30
    assert receipt["n_distinct_independent_host_plants"]==120
    assert receipt["larvae_per_host_plant"]==2
    assert receipt["plots_by_temperature"]=={"ambient":15,"elevated":15}
    assert receipt["original_temperature_is_plot_assignment"]
    assert receipt["original_status"]["both_observed"]==240


def test_synthetic_two_milkweed_by_temperature_risk_contrast_and_missing_bounds():
    row=fake_rows(one_unknown=True)
    lower=m.plot_observed_adults(row,unknown_success=False)
    upper=m.plot_observed_adults(row,unknown_success=True)
    a=m.contrast_by_plot(lower,draws=300,seed=1234)
    b=m.contrast_by_plot(upper,draws=300,seed=1234)
    assert a["observed_plot_paired_host_difference_ambient"]==pytest.approx(-1/60)
    assert b["observed_plot_paired_host_difference_ambient"]==pytest.approx(0)
    assert a["observed_plot_paired_host_difference_elevated"]==pytest.approx(-.5)
    assert b["observed_plot_paired_host_difference_elevated"]==pytest.approx(-.5)
    assert a["observed_temperature_by_species_adult_emergence_risk_interaction"]==pytest.approx(-.5+1/60)
    assert b["observed_temperature_by_species_adult_emergence_risk_interaction"]==pytest.approx(-.5)
    assert a["bootstrap_resampling_unit"].startswith("original field plot")
    assert a["source_temperature_plot_replicates_each"]==15
    assert len(a["leave_one_plot_out_values"])==30


def test_fresh_source_is_checked_by_git_blob_identity():
    raw=encode(fake_rows())
    with pytest.raises(ValueError,match="Git blob"):
        m.load_original(raw)
    assert m._git_blob(raw)!=m.ORIGINAL_BLOB


def test_original_plots_must_have_all_four_host_oe_plants():
    rows=fake_rows()
    rows[0]["OE_treatment"]="control"
    with pytest.raises(ValueError,match="eight 30-larva cells mismatch"):
        m.load_original(encode(rows),require_blob=False)


def test_original_plant_is_two_nested_larvae_not_two_hosts():
    rows=fake_rows()
    rows[0]["Plant_ID"]="P1-unrelated"
    with pytest.raises(ValueError,match="120 plants"):
        m.load_original(encode(rows),require_blob=False)


def test_original_temperature_cannot_vary_within_plot():
    rows=fake_rows()
    rows[0]["Temp"]="elevated"
    with pytest.raises(ValueError,match="eight 30-larva cells mismatch"):
        m.load_original(encode(rows),require_blob=False)


def test_invalid_source_outcomes_or_adult_before_pupation_fail_closed():
    rows=fake_rows()
    rows[0]["Surv_adult"]="99"
    with pytest.raises(ValueError,match="indicator"):
        m.load_original(encode(rows),require_blob=False)
    rows=fake_rows()
    rows[0]["Surv_pupa"]="0"
    with pytest.raises(ValueError,match="Recorded original adult with no original pupation"):
        m.load_original(encode(rows),require_blob=False)


def test_missing_original_larva_cannot_be_omitted_for_complete_case_fitness():
    rows=fake_rows()
    rows.pop()
    with pytest.raises(ValueError,match="Exactly 240"):
        m.load_original(encode(rows),require_blob=False)
    rows=fake_rows()
    rows[1]["ID"]=rows[0]["ID"]
    with pytest.raises(ValueError,match="duplicate original individual"):
        m.load_original(encode(rows),require_blob=False)


def test_host_temperature_delta_stays_zero_for_pure_infection_main_effect():
    a=m.plot_observed_adults(fake_rows(),unknown_success=False)
    b=m.plot_observed_adults(fake_rows(nonzero_main_effects=True),unknown_success=False)
    delta_a=m.contrast_by_plot(a,draws=300,seed=1)
    delta_b=m.contrast_by_plot(b,draws=300,seed=1)
    assert delta_a["observed_temperature_by_species_adult_emergence_risk_interaction"]==pytest.approx(
        delta_b["observed_temperature_by_species_adult_emergence_risk_interaction"]
    )


def test_unknown_original_outcome_is_not_mortality_or_population_causal_exotic_effect():
    raw=encode(fake_rows(one_unknown=True))
    result=m.run(raw,draws=300,seed=77,require_blob=False)
    assert result["source_units"]["original_status"]["ambiguous_stage_fate"]==1
    assert result["unknown_original_larvae_never_imputed_as_observed"]
    assert result["botanical_origin_not_randomized_across_species"]
    assert result["not_flight_capable_adult_or_reproductive_recruitment"]
    assert result["no_go_for_independent_new_mechanism_paper_from_single_published_source"]
    assert result["GEB_main_and_PR38_unmodified"]
    assert "new_empirical_butterfly_result" not in result


def test_seed_reproduces_source_plot_bootstrap_without_larval_pseudoreplication():
    plots=m.plot_observed_adults(fake_rows(),unknown_success=False)
    first=m.contrast_by_plot(plots,draws=400,seed=41)
    second=m.contrast_by_plot(plots,draws=400,seed=41)
    assert first==second
    assert not any(math.isnan(x) for x in first["stratified_whole_plot_bootstrap_95_intervals"]["interaction"])
    with pytest.raises(ValueError,match="256"):
        m.contrast_by_plot(plots,draws=1)

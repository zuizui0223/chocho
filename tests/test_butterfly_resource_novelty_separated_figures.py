"""Numeric/semantic agreement of the separated manuscript Figure 1 and Fig S2."""
from __future__ import annotations

import json
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
PATH=ROOT/"scripts/render_butterfly_resource_novelty_separated_figures.py"
spec=importlib.util.spec_from_file_location("separate_novelty_figs",PATH)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def source():
    receipts=ROOT/"provenance/reviewer_defenses/results"
    main=mod.extract_results(
        mod.load(receipts/"butterfly_plant_range_structural_null_two_chains_v0.2.json"),
        mod.load(receipts/"butterfly_climate_matched_absorption_v0.1.json"),
        mod.load(receipts/"butterfly_spatial_null_longchain_qualified_v0.2.json"),
    )
    return mod.split_figure_data(main)


def test_separated_main_figure_units_and_accounting() -> None:
    main,_=source()
    assert main["schema"] == "chocho_butterfly_fig1_geographic_novelty_accounting_v0.1"
    assert main["native_butterfly_region_cells"]==26530
    assert main["contemporary_butterfly_region_cells"]==41083
    assert main["novel_butterfly_region_cells"]==14553
    assert [v["n"] for v in main["accounting"]]==[46387,14406,14553]
    assert sum(x["n"] for x in main["accounting"])==main["incidence_budget"]
    assert "not butterfly" not in main["claim_boundary"].lower() or (
        "NOT relative butterfly abundance" in main["claim_boundary"]
    )


def test_separated_supplementary_figure_has_all_two_seed_nulls() -> None:
    main,supp=source()
    assert supp["schema"]=="chocho_butterfly_figure_S2_matched_null_v0.1"
    assert supp["incidence_budget"] == main["incidence_budget"]
    assert [x["draws_per_chain"] for x in supp["controls"]]==[499,299,99]
    assert [x["a_excess"] for x in supp["controls"]]==[
        [1086,1089], [552,558], [271,276]
    ]
    assert [x["u_deficit"] for x in supp["controls"]]==[
        [831,832], [426,431], [264,267]
    ]
    assert "nor causal variance partition" in supp["claim_boundary"]
    assert "post-hoc" in supp["claim_boundary"].lower()


def test_split_refuses_nonconservative_incidence_totals() -> None:
    m,s=source()
    inp={
        "focal_butterflies":239,
        "resource_regions":365,
        "total_candidate_incidences":75346,
        "accounting":[dict(x) for x in m["accounting"]],
        "controls":[dict(x) for x in s["controls"]],
    }
    inp["accounting"][1]["n"]+=1
    with pytest.raises(ValueError,match="drift"):
        mod.split_figure_data(inp)


def test_originally_small_spatial_effect_is_not_a_main_figure_claim() -> None:
    main,supp=source()
    assert "controls" not in main
    assert "climate" not in json.dumps(main).lower()
    assert "post-hoc" in supp["claim_boundary"].lower()
    assert len(supp["controls"])==3

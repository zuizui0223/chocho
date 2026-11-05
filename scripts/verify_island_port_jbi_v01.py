#!/usr/bin/env python3
"""Independent, dependency-free claim/denominator consistency gate for island–port
JBI Research Article v0.1, using existing frozen result receipts only.

This does NOT test causation, ecological inference assumptions, references or
human peer review. No web/API calls. Never edits the underlying GEB article.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MSS=ROOT/"manuscript/island_port_flora_assembly_20261105"
def js(path):
    return json.loads((ROOT/path).read_text(encoding="utf-8"))

def near(a,b,rtol=1e-7):
    if not math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=1e-8):
        raise AssertionError(f"Numerical source mismatch: {a} vs expected {b}")

def validate():
    q=js("results/island-mainland-resource-q1/island_mainland_summary.json")
    a=js("results/island-mainland-resource-area-v02/area_adjusted_summary.json")
    ef=js("results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json")
    cv=js("results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json")
    h=js("results/three-gate-port-introduced-host-taxa/introduced_host_taxon_port_association_audit.json")
    pair=js("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json")
    hist=js("results/q4-host-year-botanical-effort/host_first_record_botanical_effort_audit.json")
    wing=js("results/q3-wing-port-index-screen-v02/wing_port_observation_screen_audit.json")
    gate=js("results/three-gate-integrated-targeted-v04/three_gate_enrichment_audit.json")
    ledger=json.loads((MSS/"EVIDENCE_LEDGER_V01.json").read_text())
    draft=(MSS/"ISLAND_PORT_MANUSCRIPT_V01.md").read_text(encoding="utf-8")
    figures=(MSS/"FIGURES_AND_BOUNDARIES_V01.md").read_text(encoding="utf-8")
    claims={r["id"]:r for r in ledger["claims"]}
    assert len(claims)==12

    island=q["groups"]["island"]
    mainland=q["groups"]["mainland"]
    assert q["panel_total"]["panel_species"]==239
    assert (island["all_regions"],island["regions_with_contemporary_opportunity"])==(126,123)
    assert (mainland["all_regions"],mainland["regions_with_contemporary_opportunity"])==(237,236)
    assert [island[z] for z in ("native_butterfly_region_units",
        "contemporary_butterfly_region_units","added_butterfly_region_units")]==[4829,10570,5741]
    assert [mainland[z] for z in ("native_butterfly_region_units",
        "contemporary_butterfly_region_units","added_butterfly_region_units")]==[21215,29901,8686]
    near(island["pooled_increase_pct_if_native_positive"],118.88589770138745)
    near(mainland["pooled_increase_pct_if_native_positive"],40.942729201037004)
    near(a["stratum_test"]["weighted_island_minus_mainland_share"],.10702159757115921)
    near(a["stratum_test"]["within_stratum_descriptive_p_two_sided"],.0002)

    models=ef["models"]
    port_all=models["all_flora_with_gbif_total_effort"]
    port_host=models["host_taxa_with_effort"]
    assert port_all["n"]==98 and port_host["n"]==98
    near(port_all["regression_beta_per_log2_port_unit"],.24447852550694815)
    near(port_all["permutation"]["p_two_sided_centered"],.0012)
    near(port_host["regression_beta_per_log2_port_unit"],.16390595789001258)
    near(port_host["permutation"]["p_two_sided_centered"],.011)

    condition=cv["model_results"]["M1eff_port_host_taxa_given_alien_flora"]
    assert condition["n"]==98
    near(condition["coefficient"],-.041526655753407464)
    near(condition["within_level1_and_area_half_label_permutation"]["p_two_centered"],.093)
    assert h["source"]["distinct_plant_island_pairs"]==5611
    assert h["source"]["global_plant_taxa"]==564

    assert pair["eligible_islands"]==51 and pair["all_pairs"]==1275
    near(pair["spearman_host_new"],.29980576483377025)
    near(pair["spearman_host_earlier"],.4178255833828514)
    near(pair["partial_rank_host_new"],.2153535842700987)

    near(hist["primary"]["regression_year_change_per_ln1p_effort_unit"],-8.720482664116073)
    near(hist["primary"]["permutation"]["two_sided_centered_p"],.0004)
    assert wing["observed"]["observed_index_events"]==13
    assert gate["coverage"]["all_proxy_fields_and_prebutterfly_plant_index"]==12
    assert gate["coverage"]["all_proxy_distinct_species"]==2
    assert gate["coverage"]["all_proxy_distinct_butterfly_WGSRPD3_cells"]==8

    for item in (claims["port_total_introduced_flora"],
                 claims["port_distinct_introduced_host_taxa"]):
        assert item["n"]==98
    near(claims["host_not_enriched_above_all_introduced_flora"]["beta"],-.041526655753407464)

    assert "not a causal" in draft.lower() or "not a known causal" in draft.lower() or "not a causal" in figures.lower()
    for term in ("118.9%","40.9%","+0.244","+0.164","51 island","1,275"):
        assert term in draft,("Manuscript statement unexpectedly absent",term)
    assert "not" in draft.lower() and "colonization" in draft.lower()
    assert "not" in figures.lower() and "port" in figures.lower()

    abstract=draft.split("## Abstract",1)[1].split("## Introduction",1)[0]
    count=len(re.findall(r"\b[\w-]+\b",abstract))
    if count>310:
        raise AssertionError(f"Abstract exceeds provisional ~300-word JBI target: {count}")
    result={
      "status":"PASS_SOURCE_NUMERIC_AND_CLAIM_LABELS_ONLY",
      "butterfly_panel":239,
      "island_mainland_regions":[126,237],
      "effort_complete_port_flora_regions":98,
      "island_butterfly_network":[51,1275],
      "island_vs_mainland_resource_increase_pct":[118.88589770138745,40.942729201037004],
      "port_effects":[.24447852550694815,.16390595789001258,-.041526655753407464],
      "abstract_word_tokens":count,
      "not_validated":[
         "Full citations/bibliography", "Original data licensing",
         "Causal shipping or plant mediation", "Independent island-pair inference",
         "Actual larval feeding or butterfly colonisation", "Rendered figures",
         "Coauthor or editor publication-overlap approval"
      ],
    }
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    result=validate()
    print("ISLAND_PORT_MANUSCRIPT_V01_SOURCE_AUDIT",json.dumps(result,sort_keys=True))
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(result,indent=2)+"\n")

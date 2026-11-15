#!/usr/bin/env python3
"""Verify v0.3 ONE integrated GEB manuscript, sign conventions and five figures.

This source-claim audit reads original committed observations and source hashes.
It does not verify causality, refit models, license restricted inputs or claim
the review copy has been visually typeset/accepted by collaborators.
"""
from __future__ import annotations
from pathlib import Path
import csv,hashlib,json,math,re,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"manuscript/integrated_butterfly_resource_access_20261109"
M=DIR/"ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
CAP=DIR/"FIGURE_CAPTIONS_V03.md"
F2=ROOT/"figures/integrated_island_port_butterfly_v02"
F3=ROOT/"figures/integrated_island_port_butterfly_v03"

def data(path):return json.loads((ROOT/path).read_text(encoding="utf8"))
def rows(path):
    with (ROOT/path).open(newline="",encoding="utf-8-sig") as fp:return list(csv.DictReader(fp))
def close(observed,expected,tol=1e-7):
    if not math.isclose(float(observed),float(expected),rel_tol=tol,abs_tol=tol):
        raise ValueError(f"Wrong original research number: observed {observed}, expected {expected}")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    t=M.read_text(encoding="utf8");captions=CAP.read_text(encoding="utf8")
    assert t.startswith("# Human-redistributed host plants")
    assert "## Structured abstract" in t and "## 1. Introduction" in t
    abstract=t.split("## Structured abstract",1)[1].split("## 1. Introduction",1)[0]
    wc=len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",abstract))
    if wc>300:raise ValueError(f"GEB abstract exceeds 300-word bound: {wc}")
    # Kew archived source is WCVP v13 (2024-05-21), NOT the v10 DOI in the
    # frozen rWCVPdata v13 package's historical CITATION file.
    for term in ("10.34885/0yex-xv26", "WCVP v13", "21 May 2024",
                 "65bed76bae9d644ccb6ad200c05f9f5071d89e05"):
        if term not in t:raise AssertionError("Correct WCVP v13 source citation missing: "+term)
    if "Govaerts, R. (ed.) (2024)" not in t:
        raise AssertionError("Kew v13 archived dataset citation absent from references")
    for source in ("**Running title:**", "**Keywords:**"):
        if source not in t:raise AssertionError("Missing GEB research article metadata: "+source)
    running=t.split("**Running title:**",1)[1].split("\n",1)[0].strip()
    if len(running)>=40:raise AssertionError(f"GEB running title is >=40 characters: {len(running)}")
    keyword_line=t.split("**Keywords:**",1)[1].split("\n",1)[0]
    keys=[v.strip() for v in keyword_line.split(",") if v.strip()]
    if len(keys)<6 or len(keys)>10 or keys!=sorted(keys,key=str.lower):
        raise AssertionError("GEB keywords must be 6–10 and alphabetized: "+str(keys))
    if not (DIR/"REFERENCE_DATA_ATTRIBUTION_AUDIT_V03.md").exists():
        raise AssertionError("No versioned GBIF/WCVP/GEB pre-submission audit")
    body=t.split("## 1. Introduction",1)[1].split("## Data and code availability",1)[0]
    body_words=len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",body))
    for x in ("239 butterflies","126 insular","237 mainland","1,596",
              "98","155","225 butterflies","1,275","1.731%",
              "6.488%","6.48%","4.699%","5.082%","100/199","4.662%","0.300","0.2445",
              "0.1639","0.0415","17-island","2001–2010"):
        if x not in t:raise AssertionError("v0.3 essential source population or number missing: "+x)
    if "MSE" not in t or "not" not in t.lower() or "unseen" not in t.lower():
        raise AssertionError("Predictive loss direction or historical scope boundary not clearly reported")
    assert t.count("## 3. Results")==1 and t.count("## 4. Discussion")==1
    assert captions.count("## Figure ")==2
    rsrc=data("benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json")["full_resource_eligible_panel"]
    assert rsrc["species"]==239 and rsrc["total_introduced_added_species_units"]==14553
    # Separate original GBIF recording target: source-verified butterfly
    # predictive information, not observed plant arrivals or colonization.
    prediction = rows("provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_validation_all_metrics.csv")
    pairs = rows("provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_paired_comparisons.csv")
    def fixed_predictive_row(mode, model):
        q = [r for r in prediction if r["mode"] == mode and r["model"] == model]
        if len(q) != 1 or int(q[0]["species"]) != 155:
            raise AssertionError("155-species locked predictor missing: "+mode+" "+model)
        return q[0]
    for mode, expected_pct, expected_n in (
        ("temporal", 4.120172, 208482),
        ("spatial", 2.995, 208482),
        ("spatial_buffer100", 2.736, 102728),
        ("species", -0.2085, 208482),
    ):
        m3 = fixed_predictive_row(mode, "M3_both")
        m5 = fixed_predictive_row(mode, "M5_host_composition")
        assert int(m3["rows"]) == int(m5["rows"]) == expected_n
        pct = 100*(float(m3["species_equal_mse"]) - float(m5["species_equal_mse"]))/float(m3["species_equal_mse"])
        if abs(pct - expected_pct) > .008:
            raise AssertionError("Forecast performance source drift: "+mode+" "+str(pct))
        diff = [r for r in pairs if r["mode"] == mode and r["added"] == "M5_host_composition"
                and r["reference"] == "M3_both" and r["metric"] == "mse"]
        assert len(diff) == 1 and int(diff[0]["species"]) == 155
        close(float(diff[0]["difference"]),
              float(m5["species_equal_mse"]) - float(m3["species_equal_mse"]))
        if mode in ("temporal", "spatial", "spatial_buffer100"):
            assert float(diff[0]["ci_high"]) < 0
        else:
            assert float(diff[0]["ci_low"]) <= 0 <= float(diff[0]["ci_high"])
    initial = rows("provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_initially_unrecorded_comparison.csv")
    for mode, expected_pct in (("temporal", 5.231694), ("spatial", 5.284209), ("species", -.761421)):
        q = [r for r in initial if r["mode"] == mode and r["metric"] == "mse"]
        assert len(q) == 1 and (int(q[0]["species"]), int(q[0]["rows"])) == (154, 133557)
        close(float(q[0]["relative_improvement_percent"]), expected_pct, tol=.0001)
    m4_spatial = fixed_predictive_row("spatial", "M4_host_resources")
    m5_spatial = fixed_predictive_row("spatial", "M5_host_composition")
    delta_family = 100 * (float(m4_spatial["species_equal_mse"]) -
                           float(m5_spatial["species_equal_mse"])) / float(m4_spatial["species_equal_mse"])
    close(delta_family, 1.4357, tol=.001)
    assert "### 3.3. Independent butterfly-record prediction" in t
    for value in ("4.12%", "3.00%", "2.74%", "5.23%", "5.28%", "−0.21%", "1.44%"):
        assert value in t, "Missing critical original-155 validation or transfer boundary: "+value

    # Source-complete 239-resource vs 155-forecast selection boundary, NOT
    # unarchived per-species M3/M5 loss estimates.
    original_traits=rows("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
    forecast_eligibility=rows("provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_all_225_final_coverage.csv")
    gap_by_species=rows("results/island_resource_gap_filling/by_species.csv")
    elig={x["species"]:x for x in forecast_eligibility}
    family={x["species"]:x for x in original_traits}
    gaps={}
    for x in gap_by_species:
        gaps.setdefault(x["species"],{})[x["group"]]=x
    assert (len(family),len(elig),len(gaps))==(239,225,239)
    modeled={k for k,v in elig.items() if v["final_disposition"]=="MODELED_ALL_PLANNED_VALIDATIONS"}
    assert len(modeled)==155 and len(set(elig)-modeled)==70 and len(set(family)-set(elig))==14
    delta=lambda name:100*(float(gaps[name]["island"]["fill"])-float(gaps[name]["mainland"]["fill"]))
    modeled_d=sum(delta(z) for z in modeled)/155
    notmodeled_d=sum(delta(z) for z in set(elig)-modeled)/70
    close(modeled_d,-7.910170358059511)
    close(notmodeled_d,.3370650437506869)
    close(sum(float(family[z]["resolved_host_species"]) for z in modeled)/155,13.038709677419355)
    close(sum(float(family[z]["resolved_host_species"]) for z in set(elig)-modeled)/70,7.1)
    assert "13.04" in t and "7.10" in t and "−8.25" in t
    assert "155" in t and "70" in t and "239" in t
    assert (ROOT/"results/resource_forecast_cohort_selection_20261115/README.md").exists()

    # Independent source-limited GEB archival DOI test: geographic host overlap
    # cannot be interpreted as actual butterfly absence or colonization.
    q3=data("results/q3-observed-butterfly-island-opportunity-v02/q3_archival_gbif_overlap_audit.json")
    qp=q3["validation_panel"]
    assert qp["gbif_validation_species"]==32 and qp["gbif_species_with_original_records"]==31
    assert qp["sampled_species_added_island_opportunity_cells"]==920
    assert qp["island_butterfly_species_cells_with_any_gbif_record"]==159
    assert qp["recorded_butterfly_in_new_opportunity_cells"]==13
    obs=rows("results/q3-observed-butterfly-island-opportunity-v02/observed_species_island_cells.csv")
    assert len(obs)==159 and len({(r["species"],r["area_code_l3"]) for r in obs})==159
    assert all(r["observed_butterfly"]=="1" for r in obs)
    assert sum(int(r["introduced_host_resource_added"]) for r in obs)==13
    assert "13" in t and "12.87" in t and "p=0.632" in t
    audit_path=ROOT/"results/island_observed_host_degree_null_20261115/README.md"
    assert audit_path.exists() and "13" in audit_path.read_text(encoding="utf8")

    im=data("results/island-mainland-resource-q1/island_mainland_summary.json")["groups"]
    assert (im["island"]["native_butterfly_region_units"],im["island"]["contemporary_butterfly_region_units"])==(4829,10570)
    assert (im["mainland"]["native_butterfly_region_units"],im["mainland"]["contemporary_butterfly_region_units"])==(21215,29901)
    for record in rows("results/island_resource_gap_filling/pooled.csv"):
        group=im[record["geography"]]
        native=group["native_butterfly_region_units"]
        added=group["contemporary_butterfly_region_units"]-native
        total=239*group["all_regions"]
        close(record["empty"],total-native)
        close(record["coverage"],native/total)
        close(record["fill"],added/(total-native))
        close(record["gain"],added/native)
    joint=rows("results/island51_host_circuit_joint_20261106/comparison_summary.csv")
    jc={z["comparison"]:z for z in joint if z["target"]=="new_butterfly_jaccard"}
    assert set(["G_to_H","G_to_C","C_to_HC","H_to_HC","HC_to_I"])<=set(jc)
    def check_j(k,pct,improved):
        z=jc[k];ref=float(z["MSE_reference"]);new=float(z["MSE_augmented"])
        assert ref>new, f"Model {k} should be a positive MSE reduction"
        close((ref-new)/ref*100,pct,tol=1e-6)
        close(float(z["relative_MSE_reduction_percent"]),pct,tol=1e-6)
        assert int(z["islands_improved"])==improved
    geom=rows("results/island51_observed_ocean_circuit_20261106/model_comparison_summary.csv")
    geom=[x for x in geom if x["comparison"]=="base_vs_uniform" and x["target"]=="new_butterfly_jaccard"]
    assert len(geom)==1
    g=geom[0]
    close(float(g["MSE_reference"]),0.006797760771290343)
    close(float(g["MSE_augmented"]),0.0064523062608844396)
    close(float(g["relative_MSE_reduction_percent"]),5.0818868452225594)
    assert int(g["islands_improved"])==33
    check_j("C_to_HC",1.731162479353141,31)
    check_j("H_to_HC",6.4880250544198574,39)
    total=(float(jc["G_to_H"]["MSE_reference"])-float(jc["H_to_HC"]["MSE_augmented"]))/float(jc["G_to_H"]["MSE_reference"])*100
    close(total,8.6164,tol=.0001)
    cx=rows("results/island51_cross_community_marine_20261107/comparison_scope_summary.csv")
    def cross(panel,contrast):
        q=[z for z in cx if z["panel"]==panel and z["contrast"]==contrast and z["pair_class"]=="all_both_plants_indexed"]
        if len(q)!=1 or int(q[0]["unique_source_pairs"])!=990:raise ValueError("Missing distinct two-endpoint cohort")
        return float(q[0]["MSE_gain_percent"])
    close(cross("later_butterflies_retro_hosts","BASE_to_AIS"),.67755574674502794)
    close(cross("later_butterflies_retro_hosts","BASE_to_BOTH"),-4.6985954047498986)
    close(cross("earlier_butterflies","BASE_to_AIS"),-8.5406192351514054)
    close(cross("early_general_flora","BASE_to_AIS"),-1.2943024214677576)
    lat=rows("results/island51_latitude_controls_20261107/contrast_scope_summary.csv")
    w=[x for x in lat if x["target"]=="new_butterfly_jaccard" and
       x["contrast"]=="HRE_to_HRE_WIND_LAT" and x["pair_class"]=="all_pairs"]
    assert len(w)==1
    close(w[0]["MSE_reduction_percent"],-4.661759477075)
    local=rows("results/island51_AIS_local750_node_null_20261107/observed_vs_local750_node_null.csv")
    a=[x for x in local if x["target"]=="new_butterfly_jaccard" and
       x["contrast"]=="AIS_without_wind" and x["pair_class"]=="all_pairs"]
    assert len(a)==1 and int(a[0]["MSE_null_ge_observed_count"])==100
    assert int(a[0]["null_replicates"])==199
    idx=data("data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json")
    assert idx["cells"]==11475 and idx["totals"]["first_indexed_in_2011_2017"]==221
    assert idx["totals"]["first_indexed_in_2018_2022"]==96
    f17=data("provenance/reviewer_defenses/results/paired_island_port_nonport_coast_20261109/audit.json")
    assert f17["cohort"]["matched_islands"]==17 and f17["cohort"]["matched_port_nonport_pairs"]==24
    close(f17["associations"]["relative_host_enrichment_vs_new_2011_2017"]["rho"],-.08513459437657447)
    # Match the active summary figure labels to the source tables.
    fig1 = F3 / "Fig1_resource_geography_scales_v03.svg"
    labels = " ".join(ET.parse(fig1).getroot().itertext())
    for row in rows("results/island_resource_gap_filling/pooled.csv"):
        for field in ("coverage", "fill", "gain"):
            assert f"{100*float(row[field]):.1f}%" in labels
    fig2 = F3 / "Fig2_global_contributors_at_ports.svg"
    port_labels = {text.strip() for text in ET.parse(fig2).getroot().itertext()}
    port_hosts = data("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json")
    for plant in port_hosts["source_contributor_plants"]:
        assert plant["plant_name"] in port_labels
        assert str(plant["n_port_buffers_with_pre2011_record"]) in port_labels
        assert f"{100*plant['GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant']:.2f}" in port_labels
    for figure in (fig1, fig2):
        assert figure.with_suffix(".pdf").read_bytes().startswith(b"%PDF")
        assert figure.with_suffix(".png").read_bytes().startswith(b"\x89PNG")
        assert figure.with_suffix(".png").name in body
    # Retain the legacy v0.2 source-file checks.
    for i in range(1,5):
        fs=list(F2.glob(f"Fig{i}_*.png"))
        if len(fs)!=1 or not fs[0].read_bytes().startswith(b"\x89PNG\r\n\x1a\n"):
            raise ValueError(f"Missing/invalid retained Fig {i}")
        if not fs[0].with_suffix(".svg").exists():raise ValueError(f"Missing vector retained Fig {i}")
    receipt=data("figures/integrated_island_port_butterfly_v03/FIGURE5_V03_SOURCE_RECEIPT.json")
    assert receipt["status"]=="SOURCE_LOCKED_FIGURE_RENDERED"
    for src,digests in receipt["digests"].items():
        path=ROOT/receipt["sources"][src]
        assert sha(path)==digests["sha256"],"changed scientific figure source "+src
    for name,desc in receipt["figures"].items():
        p=F3/name
        if sha(p)!=desc["sha256"] or p.stat().st_size!=desc["bytes"]:
            raise ValueError("Figure 5 hash mismatch "+name)
        if p.suffix==".svg":ET.parse(p)
    close(receipt["values"]["geometry_only_beyond_distance"],5.0818868452225594)
    assert "uniform" in receipt["sources"]
    close(receipt["values"]["host_given_marine"],float(jc["C_to_HC"]["relative_MSE_reduction_percent"]))
    close(receipt["values"]["marine_given_host"],float(jc["H_to_HC"]["relative_MSE_reduction_percent"]))
    assert "FIGURE_CAPTIONS_V03.md" in t and "Figure 2" in captions
    report={"status":"V03_SINGLE_PAPER_SCIENTIFIC_SOURCE_VERIFICATION_PASS",
            "abstract_word_count":wc,"body_word_count_approx":body_words,
            "independently_verified_model_improvements_pct":{"geometry_only_ocean_beyond_geo_and_history":5.0818868452225594,
              "host_beyond_marine":1.731162479353141,
              "marine_beyond_host":6.4880250544198574,
              "geo_baseline_to_joint":total,
              "double_unseen_AIS_and_wind_beyond_host":-4.6985954047498986},
            "main_figures_verified":2,"new_Figure5_input_and_output_hashes_verified":True,
            "warnings":["The 51-island host+marine positive increments are retrospective and not causal",
               "Historic host list may contain 2001–2024 records after target butterflies",
               "The strict 990-pair geographic holdout is a different population/estimand",
               "Manuscript body may exceed final target word budget depending on journal definitions",
               "No independent raw GBIF DOI, third-party license, coauthor signoff or journal submission check"]}
    print(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2))
    return report
if __name__=="__main__":run()

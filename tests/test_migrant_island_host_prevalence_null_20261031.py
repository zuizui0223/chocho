"""Replay 51-island migration-roster × pooled port flora similarity checks.

A positive ecological association is documented, but does not support a
migrant-specific effect beyond this exploratory prevalence-matched label null.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031"
SRC=D/"analyze_migrant_island_network.py"
SPEC=importlib.util.spec_from_file_location("butterfly_migrant_island_network",SRC)
assert SPEC and SPEC.loader
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def close(left,right):
    assert math.isclose(left,right,rel_tol=0,abs_tol=2e-9),(left,right)


def test_positive_migrant_similarity_but_no_prevalence_null_excess():
    saved=json.loads((D/"migrant_island_host_prevalence_null_results.json").read_text(encoding="utf8"))
    got=MODULE.analyze(draws=999)
    assert saved["status"]==got["status"]=="POSITIVE_ISLAND_SIMILARITY_NOT_MIGRANT_SPECIFIC_ABOVE_PREVALENCE_NULL"
    for k,v in [("islands",51),("island_pairs",1275),
                ("butterfly_first_index_2011_2017_island_events",302),
                ("butterfly_species_first_indexed_2011_2017",106),
                ("author_roster_migrant_selected_species_with_new_island_records",63),
                ("not_listed_unknown_species_with_new_records",43)]:
        assert got[k]==saved[k]==v
    assert got["source_git_blobs"]==saved["source_git_blobs"]
    assert got["source_prevalence_bins"]==saved["source_prevalence_bins"]
    for k in ("first_2011_2017_all_species",
              "first_2011_2017_author_migrant_roster",
              "first_2011_2017_not_listed_UNKNOWN_not_nonmigrants"):
        for metric,value in saved[k].items():
            if metric.endswith("butterfly"):
                close(got[k][metric],value)
            else:
                assert got[k][metric]==value
    migrant=got["first_2011_2017_author_migrant_roster"]
    assert migrant["n_pairs"]==1272
    close(migrant["spearman_host_vs_butterfly"],0.2204520274075129)
    close(migrant["partial_rank_distance_prior_butterfly"],0.17021539879216102)
    for metric in ("spearman_host_vs_butterfly",
                   "partial_rank_distance_prior_butterfly","n_pairs"):
        observed=got["permutation"]["statistics"][metric]
        archived=saved["permutation"]["statistics"][metric]
        assert observed["draws"]==archived["draws"]==999
        assert observed["draws_ge_observed"]==archived["draws_ge_observed"]
        close(observed["mean"],archived["mean"])
        for name in ("reference95","range"):
            for x,y in zip(observed[name],archived[name]):
                close(x,y)
    assert saved["permutation"]["statistics"]["spearman_host_vs_butterfly"]["draws_ge_observed"]==399
    assert saved["permutation"]["statistics"]["partial_rank_distance_prior_butterfly"]["draws_ge_observed"]==340
    for name in ("simple_rho_range","partial_rank_range"):
        for x,y in zip(got["migrant_leave_one_island_sensitivity"][name],
                       saved["migrant_leave_one_island_sensitivity"][name]):
            close(x,y)
    assert saved["permutation"]["posthoc_descriptive_only"] is True
    assert "not sedentary" in saved["explicit_limits"][1].lower()
    assert any("postdating" in s for s in saved["explicit_limits"])

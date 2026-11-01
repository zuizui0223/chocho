"""Unit and source-replay tests for VOС cargo-specific first-indexed island butterflies."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/analyze_VOC_cargo_vs_recorded_island_butterfly_species.py"


def load():
    spec=importlib.util.spec_from_file_location("voc_island_real_butterflies",SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_route_specific_genus_link_requires_independently_indexed_butterfly():
    m=load()
    butterflies={"12":{"eurema hecabe","papilio demoleus"},
                 "19":{"eurema hecabe","vanessa cardui"},
                 "31":{"mycalesis perseus","vanessa cardui"}}
    compatible={"indigo":{"eurema hecabe"},
                "rice":{"mycalesis perseus","vanessa cardui"},
                "sugar":{"eurema hecabe"}}
    r=m.source_route_stats(("12","19"),{"indigo"},butterflies,compatible)
    assert r["has_either_endpoint_with_compatible_record"]==1
    assert r["has_both_endpoints_each_with_compatible_record"]==1
    assert r["has_same_compatible_recorded_butterfly_in_both_islands"]==1
    assert r["same_recorded_butterfly_species_both"]==["eurema hecabe"]
    missing=m.source_route_stats(("12","31"),{"indigo"},butterflies,compatible)
    assert missing["has_both_endpoints_each_with_compatible_record"]==0
    assert missing["has_same_compatible_recorded_butterfly_in_both_islands"]==0
    no_data=m.source_route_stats(("12","31"),set(),butterflies,compatible)
    assert no_data["has_either_endpoint_with_compatible_record"]==0


def test_all_5040_seven_cargo_packet_assignments_exactly_enumerated():
    m=load()
    # Seven positive archive links, no invented negative voyage links.
    links=[("12","19"),("12","31"),("12","44"),
           ("17","19"),("19","31"),("19","43"),("43","44")]
    packets=[{"indigo"},{"rice"},{"none"},{"rice","indigo"},
             {"rice"},{"indigo","none"},{"indigo"}]
    edges={pair:{"historic_cargo_indexed_rows":1,
                 "globally_host_compatible_categories":packet}
           for pair,packet in zip(links,packets)}
    first={"12":{"a"},"17":{"b"},"19":{"a","b"},
           "31":{"b"},"43":{"a"},"44":{"c"}}
    by_cat={"indigo":{"a"},"rice":{"b"},"none":set()}
    _,summary,n=m.exact_permutation_diagnostic(edges,first,by_cat)
    assert n==5040
    assert summary["has_either_endpoint_with_compatible_record"]["observed"]<=7
    assert summary["has_same_compatible_recorded_butterfly_in_both_islands"]["observed"]<=7
    for r in summary.values():
        assert sum(r["null_histogram"].values())==5040
        assert 0<=r["conditional_finite_permutation_right_tail_fraction"]<=1
        assert r["null_min"]<=r["observed"]<=r["null_max"]


def test_original_SIX_islands_49_first_indexed_events_and_seven_cargo_links():
    m=load()
    edges,b,lock=m.load_source_receipts()
    indexed,qa=m.independently_indexed_butterflies()
    assert len(edges)==9
    assert sum(v["historic_cargo_indexed_rows"]>0 for v in edges.values())==7
    assert qa["global_original_new_events"]==5665
    assert sum(len(v) for v in indexed.values())==49
    assert indexed["12"]&indexed["19"]&{"eurema hecabe"}=={"eurema hecabe"}
    assert indexed["19"]&indexed["31"]&{"mycalesis perseus"}=={"mycalesis perseus"}
    assert lock["null"]["n_exact_assignments"]==5040

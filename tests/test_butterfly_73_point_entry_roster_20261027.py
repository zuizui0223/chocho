"""Source-independent controls for exact PR72 positive roster extraction."""
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from derive_butterfly_73_point_entry_roster_20261027 import select_roster, summarize
from analyze_butterfly_niche_entry_country_20261026 import PRIMARY


def row(species, country, level1, kind, native=False):
    return {
        "species": species, "country": country, "level1": level1,
        "native": native, "added": .2,
        "analog": .19 if kind == "analogue" else .01,
        "novel": .01 if kind == "analogue" else .19,
        "country_km2": 100,
    }


def source_fixture():
    rows = [
        row("Test butterfly", "CU", "8", "novel"),
        row("Test butterfly", "DO", "8", "analogue"),
        row("Test butterfly", "HT", "8", "novel"),
        row("Test butterfly", "DE", "1", "novel"),
        row("Test butterfly", "SE", "1", "novel"),
        row("Test butterfly", "FR", "1", "analogue", native=True),
    ]
    count = Counter({("Test butterfly", "CU"): 4,
                     ("Test butterfly", "DO"): 1,
                     ("Test butterfly", "DE"): 3})
    for country in ("CU", "DO", "HT", "DE", "SE"):
        for s in ("Effort one", "Effort two", "Effort three"):
            count[s, country] = 1
    return rows, {("all_accepted_bases", "pooled"): count}, {
        "all_accepted_bases": {"Test butterfly", "Effort one",
                               "Effort two", "Effort three"}
    }


def test_reproduces_only_positive_paired_continent_cases():
    rows, counts, complete = source_fixture()
    hits, strata, audit = select_roster(
        rows, counts, complete, {"Test butterfly": 987654}
    )
    assert len(hits) == 2
    assert [(x["country"], x["resource_class"]) for x in hits] == [
        ("CU", "novel"), ("DO", "analogue")
    ]
    assert {x["gbif_accepted_species_key"] for x in hits} == {987654}
    assert len(strata) == 1
    assert (strata[0]["N"], strata[0]["A"], strata[0]["K"],
            strata[0]["O"]) == (3, 2, 2, 1)
    assert audit["candidate_source_cells"] == 5
    result = summarize(hits, strata)
    assert result["positive_country_species_cases"] == 2
    assert result["distinct_butterfly_species"] == 1
    assert result["source_positive_cases_above_cap_5000"] == 0


def test_no_own_observations_as_effort():
    rows, counts, complete = source_fixture()
    counts["all_accepted_bases", "pooled"] = Counter({
        ("Test butterfly", "CU"): 1, ("Test butterfly", "DO"): 1
    })
    hits, strata, audit = select_roster(rows, counts, complete, {
        "Test butterfly": 987654
    })
    assert hits == [] and strata == []


def test_exact_taxon_key_must_exist_and_be_integer():
    rows, counts, complete = source_fixture()
    with pytest.raises(ValueError, match="missing exact source"):
        select_roster(rows, counts, complete, {})
    with pytest.raises(ValueError, match="missing exact source"):
        select_roster(rows, counts, complete, {"Test butterfly": "not-a-key"})


def test_excludes_source_native_guard_even_when_observed():
    rows, counts, complete = source_fixture()
    rows = [x for x in rows if x["country"] in ("CU", "DO", "FR")]
    counts["all_accepted_bases", "pooled"][("Test butterfly", "FR")] = 8
    hits, _, _ = select_roster(rows, counts, complete, {
        "Test butterfly": 987654
    })
    assert all(c["country"] != "FR" for c in hits)


def test_positives_are_species_country_units_not_index_record_total():
    rows, counts, complete = source_fixture()
    counts["all_accepted_bases", "pooled"][("Test butterfly", "CU")] = 5001
    hits, strata, _ = select_roster(rows, counts, complete, {
        "Test butterfly": 987654
    })
    report = summarize(hits, strata)
    assert report["positive_country_species_cases"] == 2
    assert report["source_positive_cases_above_cap_5000"] == 1
    assert report["source_positive_records_above_cap_5000"] == 5001
    assert report["original_GBIF_indexed_record_total_in_positive_cases"] == 5002

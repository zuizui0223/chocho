"""Cross-check archived independent island checklist denominator and HOSTS matches."""
import csv
import json
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
BASE = (ROOT / "provenance" / "reviewer_defenses" / "results"
        / "independent_island_butterfly_host_coverage_20261117")


def load():
    with (BASE / "species_island_coverage.csv").open(encoding="utf-8", newline="") as f:
        pair = list(csv.DictReader(f))
    summary = json.loads((BASE / "summary.json").read_text(encoding="utf-8"))
    return pair, summary


def test_audit_source_compiles():
    path = BASE / "audit_independent_checklist_hosts.py"
    compile(path.read_text(encoding="utf-8"), str(path), "exec")


def test_independent_denominator_and_host_match():
    pair, s = load()
    assert len(pair) == 227
    assert len({(r["island"], r["scientific_name"]) for r in pair}) == len(pair)
    assert len({r["scientific_name"] for r in pair}) == 42
    assert len({r["island"] for r in pair}) == 18
    assert sum(int(r["hosts_named_host_genus"]) for r in pair) == 168
    assert sum(int(r["hosts_named_host_species"]) for r in pair) == 162
    by_species = {r["scientific_name"]: r for r in pair}
    assert sum(int(r["hosts_named_host_genus"]) for r in by_species.values()) == 19
    assert sum(int(r["hosts_named_host_species"]) for r in by_species.values()) == 18
    assert s["species_island"]["n"] == len(pair)
    assert s["independent_unique_species"]["n"] == len(by_species)
    assert s["independent_unique_species"]["hosts_named_genus"] == 19


def test_archipelagos_and_madeira_status():
    pair, s = load()
    assert Counter(r["archipelago"] for r in pair) == {
        "Azores": 61, "Canary Islands": 141, "Madeira": 25,
    }
    madeira = [r for r in pair if r["island"] == "Madeira, Madeira"]
    native = [r for r in madeira if r["madeira_origin"] == "native_nonendemic"]
    endemic = [r for r in madeira if r["madeira_origin"] == "island_endemic_species"]
    assert len(native) == 7
    assert all(r["hosts_named_host_genus"] == "1" for r in native)
    assert {r["scientific_name"] for r in endemic} == {
        "Hipparchia maderensis", "Pararge xiphia", "Gonepteryx maderensis"
    }
    assert all(r["hosts_named_host_genus"] == "0" for r in endemic)
    old_extinct = [r for r in madeira
                   if r["madeira_origin"] == "historic_extinct_endemic_taxon"]
    assert len(old_extinct) == 1
    assert old_extinct[0]["scientific_name"] == "Pieris wollastoni"
    assert s["madeira_island_origin_strata"]["island_endemic_species"]["n"] == 3


def test_no_auto_synonym_or_unobserved_absence():
    pair, s = load()
    rows = {r["scientific_name"]: r for r in pair}
    assert rows["Cyclyrius webbianus"]["unverified_alias_candidate"] == "Leptotes webbianus"
    assert rows["Cyclyrius webbianus"]["hosts_named_host_genus"] == "0"
    assert rows["Pontia daplidice"]["unverified_alias_candidate"] == "Pieris daplidice"
    assert rows["Pontia daplidice"]["hosts_named_host_genus"] == "0"
    assert s["duplicates_present_and_irregular"] == [
        {"island": "Madeira, Madeira", "taxon": "Issoria lathonia"}
    ]


def test_independent_species_vs_original_frozen239_selection():
    from importlib.util import module_from_spec, spec_from_file_location
    script = BASE / "audit_frozen239_representativeness.py"
    compile(script.read_text(encoding="utf-8"), str(script), "exec")
    data = json.loads((BASE / "frozen239_selection_summary.json").read_text(encoding="utf-8"))
    with (BASE / "independent_42_species_vs_frozen239.csv").open(
        encoding="utf-8", newline=""
    ) as f:
        species = list(csv.DictReader(f))
    assert len(species) == 42
    assert len({r["scientific_name"] for r in species}) == 42
    assert sum(int(r["hosts_named_genus"]) for r in species) == 19
    assert sum(int(r["in_frozen239"]) for r in species) == 12
    assert sum(int(r["hosts_named_genus"]) and not int(r["in_frozen239"])
               for r in species) == 7
    assert data["source_denominator"] == 42
    assert data["in_frozen239"] == 12
    assert data["nonselected_named_host"] == 7
    assert data["nonselected_no_named_host"] == 23
    madeiran_endemics = [
        r for r in species if r["madeira_origin"] == "island_endemic_species"
    ]
    assert len(madeiran_endemics) == 3
    assert all(r["in_frozen239"] == "0" for r in madeiran_endemics)
    native = [r for r in species if r["madeira_origin"] == "native_nonendemic"]
    assert len(native) == 7
    assert sum(int(r["in_frozen239"]) for r in native) == 6
    # Executes second gate entirely from the already archived, frozen CSVs.
    # No web or extra packages are needed.
    import runpy
    result = runpy.run_path(str(script), run_name="audit_test_import")
    recomputed = result["audit"]()
    assert recomputed["in_frozen239"] == 12
    assert recomputed["nonselected_named_host"] == 7

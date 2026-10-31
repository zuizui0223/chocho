"""Fixed tracked-source checks for species-level host/remote-gap eligibility."""
import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030"
SRC = BASE / "audit_species_host_remote_gap_overlap.py"
spec = importlib.util.spec_from_file_location("host_remote_gap_audit", SRC)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_real_source_species_overlap_snapshot():
    computed = mod.audit(ROOT)
    recorded = json.loads((BASE / "species_host_remote_gap_overlap.json").read_text(encoding="utf-8"))
    assert computed == recorded
    assert computed["species_with_any_remote_gap"] == 3
    assert computed["remote_gap_events_in_these_selected_species_any_island"] == 12
    assert computed["remote_gap_cohort"]["species_with_at_least_one_remote_gap"] == 86
    assert computed["joint_island_species_year_host_and_gap_event_status"].startswith("NOT_OBSERVABLE")
    by_name = {r["butterfly"]: r["remote_gap_500km_both_period_events_across_all_islands"]
               for r in computed["selected_early_host_species"]}
    assert by_name == {
        "Colias croceus": 0, "Colias hyale": 1, "Nymphalis antiopa": 0,
        "Pieris rapae": 3, "Vanessa cardui": 8,
    }

def test_bad_source_sha_does_not_silently_pass(tmp_path):
    fake = tmp_path / mod.LINKS
    fake.parent.mkdir(parents=True)
    fake.write_text("species,island_a,island_b,host_name\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Git blob SHA changed"):
        mod.locked_rows(tmp_path, mod.LINKS)

"""Public CI: preserve 230 original host temporal comparisons without re-publishing edges.

Do not treat source non-detection as host absence; do not call re-recording
actual persistence/colonization. CI cannot refetch the private author source.
"""
from __future__ import annotations
from pathlib import Path
import ast,csv,hashlib,json,math
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/original230_dated_host_join_20261031"
EVENTS=ROOT/"provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv"
CANDIDATES=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv"
PLANTS=ROOT/"provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
OUT=D/"original230_prior_host_recurrence_audit.json"
CODE=D/"run_dated_host_recurrence.py"
def blob(f):
 b=f.read_bytes()
 return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def test_private_source_temporal_host_cohort_denominators():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("This exploratory original230 host and port source is withheld from anonymous GEB review bundle")
 x=json.loads(OUT.read_text(encoding="utf8"))
 ast.parse(CODE.read_text(encoding="utf8"))
 pins={Path(p).name:v for p,v in x["public_repo_git_blobs"].items()}
 assert blob(EVENTS)==pins[EVENTS.name]
 assert blob(CANDIDATES)==pins[CANDIDATES.name]
 assert blob(PLANTS)==pins[PLANTS.name]
 source=json.loads(PLANTS.read_text(encoding="utf8"))
 assert x["new_15_port_gzip_sha256"]==source["new_15port_gzip_sha256"]
 assert x["original_23_port_gzip_sha256"]==source["private_23port_gzip_sha256"]
 assert x["author_original_host_edges_sha256"]=="e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22"
 assert len(x["author_plant_candidates_sha256"])==64
 t=x["taxonomic_match"]
 assert t["source_original_edges"]==2524 and t["source_original_butterflies"]==230
 assert t["eligible_plant_names"]==1642
 assert sum(t["resolution_status_counts"].values())==1642
 assert t["resolved_butterfly_host_edges"]==2396
 assert t["resolved_original_butterfly_species"]==228
 e=x["butterfly_event_cohort"]
 events=list(csv.DictReader(EVENTS.open(encoding="utf8",newline="")))
 assert len(events)==e["event_rows_original_all_islands"]==5665
 ids={p["island_id"] for p in source["islands"]}
 subset=[r for r in events if r["landmass_id"] in ids and 2011<=int(r["first_record_year"])<=2017]
 assert len(subset)==e["events_in_17_islands_2011_2017"]==121
 assert len({r["species"] for r in subset})==e["species_with_2011_2017_new_index"]==64
 assert e["target_with_resolved_source_catalogue"]==113
 assert e["source_lost_or_unresolved"]==8
 assert e["retrospective_host_group_column_NOT_USED_TO_ASSIGN_HISTORIC_HOST"] is True
def test_original230_positive_structural_overlap_vs_negative_re_record():
 if "_anonymous_stage" in ROOT.parts:pytest.skip("Private author source provenance excluded from manuscript export")
 x=json.loads(OUT.read_text(encoding="utf8"))
 s=x["rates"]["all_source_eligible"]
 assert s["eligible_first_record_events"]==113 and s["species"]==60 and s["islands"]==17
 assert s["historically_recorded_compatible_host_events"]==34
 assert s["no_indexed_compatible_host_events_NOT_TRUE_ABSENCE"]==79
 a=s["any_later_re_record"]
 assert (a["matched_yes"],a["matched_n"],a["no_index_yes"],a["no_index_n"])==(26,34,70,79)
 assert math.isclose(a["matched_rate"],26/34)
 assert math.isclose(a["no_index_rate"],70/79)
 assert math.isclose(a["difference"],26/34-70/79)
 strict=x["rates"]["botanical_50keys_3years_exploratory"]
 assert strict["eligible_first_record_events"]==54
 assert strict["historically_recorded_compatible_host_events"]==33
 null=x["host_network_structural_null"]
 assert null["butterfly_degree_host_degree_and_species_x_family_breadth_preserved"] is True
 assert null["n_draws"]==299
 assert null["null_draws_ge_observed_fraction"]==7
 assert null["null_draws_ge_observed_delta"]==169
 assert math.isclose(null["observed_prior_host_index_event_fraction"],34/113)
 assert null["null_fraction_mean"]<null["observed_prior_host_index_event_fraction"]
 assert null["null_fraction_95pct"][0]<null["null_fraction_95pct"][1]
 assert x["no_unobserved_plant_treated_as_absent"]
 assert x["raw_original_host_edges_or_individual_plant_point_occurrences_published"] is False
 assert x["does_not_identify_true_colonization_or_island_source"]

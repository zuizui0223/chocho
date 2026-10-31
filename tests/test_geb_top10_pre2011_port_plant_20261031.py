"""GEB top-10 regional-resource contributor, historical pre-butterfly port plants.

CI verifies source identities and frozen aggregate outputs, without
distributing historical original GBIF plant records or source host edges.
"""
from pathlib import Path
import json,hashlib,math,ast
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/geb_port_island_realization_20261031"

def blob(path):
 data=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def test_GEB_resource_top10_in_actual_historic_port_sources():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Historical port host sources are exploratory; not in the independent GEB review pack")
 a=json.loads((D/"GEB_top10_pre2011_port_plant_overlap.json").read_text(encoding="utf8"))
 source=ROOT/"provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json"
 assert blob(source)==a["source_GEB_contributor_git_blob"]
 original=json.loads(source.read_text(encoding="utf8"))
 assert (a["study_islands"],a["study_port_proxies"],a["GEB_top10_species"])==(17,38,10)
 assert a["top10_species_with_plant_index_at_study_island_ports"]==8
 assert math.isclose(a["GEB_top10_fraction_of_added_units"],original["plant_species"]["top_k_share"]["10"],abs_tol=1e-12)
 plants=a["source_contributor_plants"]
 assert len(plants)==10
 assert sum(x["n_study_islands_with_pre2011_record"]>0 for x in plants)==8
 assert all(x["source_strict_accepted_taxon_match"] for x in plants)
 assert math.isclose(sum(x["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"] for x in plants
  if x["n_study_islands_with_pre2011_record"]),a["GEB_top10_species_WITH_any_port_index_total_fraction_of_GLOBAL_added_units"],abs_tol=1e-12)
 assert math.isclose(a["GEB_top10_species_WITH_any_port_index_total_fraction_of_GLOBAL_added_units"],.20701372105789737,abs_tol=1e-12)
 index={x["plant_name"]:x for x in plants}
 assert index["Medicago sativa"]["n_study_islands_with_pre2011_record"]==3
 assert index["Medicago sativa"]["n_port_buffers_with_pre2011_record"]==7
 assert index["Dactylis glomerata"]["n_study_islands_with_pre2011_record"]==5
 assert index["Dactylis glomerata"]["n_port_buffers_with_pre2011_record"]==9
 assert index["Cynodon dactylon"]["first_indexed_plant_year"] is None
 assert index["Oryza sativa"]["first_indexed_plant_year"] is None
 assert a["no_2001_2010_site_specific_alien_status"]
 assert a["top10_contribution_share_not_fraction_of_global_55pct_realized_on_islands"]
 ast.parse((D/"audit_GEB_top10_on_early_island_ports.py").read_text(encoding="utf8"))

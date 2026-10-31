"""Historical GBIF accepted larval-host names precede original butterfly index years.

Proof of indexed chronology ONLY, not dispersal, actual host use or colonization.
"""
from pathlib import Path
import hashlib,json
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/dated_focal_host_taxonomy_20261031"
BF=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/full_original_island_chain_years_observed_20261030.json"
LINKS=ROOT/"provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_links.csv"
def blob(f):
    b=f.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()
def test_gbif_taxa_original_source_annual_host_before_butterfly():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory maritime/port data excluded from independent anonymous GEB manuscript")
    taxa_path=D/"focal_historical_host_taxon_year_audit.json"
    p=json.loads(taxa_path.read_text(encoding="utf8"))
    b=json.loads(BF.read_text(encoding="utf8"))
    o=json.loads((D/"focal_host_prior_to_butterfly_20261031.json").read_text(encoding="utf8"))
    assert blob(LINKS)=="88db773c74dea9cd2d58bde0b67a394d93226869"
    assert blob(taxa_path)==o["plants_source_git_blob"]
    assert blob(BF)==o["butterflies_source_git_blob"]
    assert (o["butterfly_links_checked"],o["butterfly_species_island_endpoints"],o["island_taxon_endpoint_types"])==(3,6,4)
    assert p["source_all_23port_specieskey_year_gzip_sha256"]=="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
    assert {x["query_name"]:x["gbif_species_key"] for x in p["names"]}=={
        "Medicago sativa":9151957,"Alnus incana":9148577}
    assert all(x["rank"]=="SPECIES" and x["status_at_query"]=="ACCEPTED" and x["exact_name_match"] and x["confidence"]==100 for x in p["names"])
    host={(x["host_name"],x["island_id"]):x for x in p["island_taxon_observations"]}
    assert {key:x["first_indexed_year"] for key,x in host.items()}=={
        ("Medicago sativa","169"):2002,("Medicago sativa","276"):2001,
        ("Alnus incana","276"):2001,("Alnus incana","715"):2002}
    expected=[x for x in b["rows"] if x["island_port_representative_distance_km"]<500]
    assert len(expected)==3
    for x in o["focal_links"]:
        orig=next(t for t in expected if t["butterfly"]==x["butterfly"])
        assert x["source_catalogue_compatible_host"]==orig["known_host_before_window"]
        for end in x["island_record_comparisons"]:
            src=next(z for z in orig["endpoints"] if z["island"]==end["island"])
            assert src["first_indexed_year_from_original_year_table"]==end["butterfly_first_indexed_year"]
            assert host[x["source_catalogue_compatible_host"],end["island"]]["first_indexed_year"]==end["host_first_indexed_plant_year"]
            assert end["first_documented_host_precedes_first_documented_butterfly"]
            assert 14<=end["lead_years_of_host_index_before_butterfly_index"]<=19
    assert o["no_historical_island_introduced_status_known"]
    assert o["no_historical_actual_vessel_or_butterfly_island_to_island_route"]
    assert o["no_nonhost_negative_control_observed"]

"""Public consistency checks for 23-port taxonomic-sampling audit."""
from pathlib import Path
import json,hashlib,math
import pytest

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_host_controls_20261031"

def git_blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def test_port_rank_sampling_bias_from_gbif_facets():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory plant source excluded from GEB anonymous review bundle")
    source=json.loads((D/"daily_botanical_all23ports_20261031.json").read_text())
    x=json.loads((D/"port_rank_taxonomic_coverage_20261031.json").read_text())
    assert git_blob(D/"daily_botanical_all23ports_20261031.json")==x["source_all23_manifest_git_blob_sha"]
    assert source["specieskey_year_artifact"]["sha256"]==x["raw_all23_specieskey_year_gzip_sha256"]
    assert source["specieskey_year_artifact"]["rows"]==9788
    assert x["all_historical_index_ports"]==len(x["per_port_aggregated"])==23
    assert sum(p["observed_taxon_keys"]>0 for p in x["per_port_aggregated"])==22
    port_map={str(z["port_id"]):z for z in source["by_port"]}
    per_port=x["per_port_aggregated"]
    assert set(z["port_id"] for z in per_port)==set(port_map)
    assert all(z["indexed_taxon_identified_record_count"]<=port_map[z["port_id"]]["new_annual_query_index_sum"] for z in per_port)
    islands=x["island_coverage_comparison"]
    assert len(islands)==9
    for i in islands:
        assert 0<i["largest_plant_record_count_port_n_taxa"]<=i["across_port_species_key_union"]
        assert math.isclose(i["species_key_fraction_recovered_by_max_record_port"],
            i["largest_plant_record_count_port_n_taxa"]/i["across_port_species_key_union"],abs_tol=1e-6)
        assert i["species_keys_not_indexed_in_max_record_port"]==i["across_port_species_key_union"]-i["largest_plant_record_count_port_n_taxa"]
    imap={i["island_id"]:i for i in islands}
    assert (imap["81"]["largest_plant_record_count_port_id"],imap["81"]["largest_plant_record_count_port_n_taxa"],imap["81"]["across_port_species_key_union"])==("1518",4,142)
    assert imap["81"]["richest_single_port_n_taxa"]==52
    assert (imap["276"]["largest_plant_record_count_port_id"],imap["276"]["richest_species_key_port_id"])==("1436","1433")
    assert x["no_actual_host_non_detection_inference"]

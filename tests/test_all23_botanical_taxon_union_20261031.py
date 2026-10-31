"""Guard the 23-port dated botanical source audit without claiming host absence."""
import json,hashlib
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_host_controls_20261031"
def test_all23_port_taxon_union_snapshot():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Botanical source provenance deliberately excluded from blinded GEB submission")
    a=json.loads((D/"daily_botanical_all23ports_20261031.json").read_text())
    b=json.loads((D/"all23_botanical_taxon_union_20261031.json").read_text())
    assert a["source_ports_all"]==23 and a["source_islands_all"]==9
    assert a["annual_port_year_rows"]==230
    assert a["specieskey_year_artifact"]["rows"]==9788
    assert b["source_private_gzip_sha256"]==a["specieskey_year_artifact"]["sha256"]
    assert b["source_allport_manifest_git_blob_sha"]==(
        lambda f:hashlib.sha1(b"blob "+str(len(f)).encode()+b"\0"+f).hexdigest())(
            (D/"daily_botanical_all23ports_20261031.json").read_bytes())
    assert sum(x["indexed_species_key_union_across_all_ports"] for x in a["by_island_species_key_unions"])==sum(
        x["unique_indexed_species_keys_across_all_ports"] for x in b["island_union_of_indexed_plant_species_keys"])
    ids={x["island_id"]:x["unique_indexed_species_keys_across_all_ports"]
        for x in b["island_union_of_indexed_plant_species_keys"]}
    assert (ids["81"],ids["214"],ids["276"],ids["951"])==(142,219,962,196)
    pairs={(x["island_a"],x["island_b"]):x for x in b["prespecified_priority_and_focal_pair_overlap"]}
    assert len(pairs)==7
    assert pairs["81","214"]["shared_indexed_species_keys_pre2011"]==47
    assert pairs["276","321"]["species_first_indexed_earlier_on_a"]==83
    assert b["historical_introduced_status_unknown"]
    assert b["counts_from_overlapping_port_buffers_are_not_independent"]
    assert b["no_confirmed_nonhost_control"]

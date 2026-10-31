"""Check near-island record-order baseline against the source-locked public 1275-pair table.

The original private year parquet is not in CI. This test validates support
counts, pair identities and arithmetic only; it cannot prove the raw chronology.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
B=ROOT/"provenance/reviewer_defenses/results"
PAIRS=B/"island_host_network_20261030/island_host_network_pairs.csv"
RESULT=B/"maritime_host_steppingstone_20261030/near_island_record_order_baseline_20261030.json"

def blob(path):
    data=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def test_near_pairs_and_shared_species_counts_match_public_source():
    j=json.loads(RESULT.read_text(encoding="utf-8"))
    assert blob(PAIRS)==j["source_island_pair_csv_git_blob"]=="5663d2ffe42419f5c4496b92c9f3bbf76175907b"
    assert j["source_private_parquet_sha256"]=="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
    pairrows=list(csv.DictReader(PAIRS.open(encoding="utf-8",newline="")))
    assert len(pairrows)==j["source_total_island_pairs"]==1275
    published={(r["island_a"],r["island_b"]):r for r in pairrows if float(r["port_distance_km"])<200}
    assert len(published)==j["near_island_pairs"]==14
    near=j["pair_summaries"]
    assert len(near)==len(published)
    for r in near:
        p=published[(r["island_a"],r["island_b"])]
        assert abs(float(p["port_distance_km"])-r["port_representative_distance_km"])<1e-9
        assert int(p["shared_new_species"])==r["shared_new_butterfly_species_pair_count"]
        n=r["shared_new_butterfly_species_pair_count"]
        distinct=r["distinct_first_record_year_species_pairs"]
        early=r["upstream_re_record_before_downstream_first_year"]
        late=r["downstream_re_record_within_five_years"]
        both=r["upstream_and_downstream_re_record"]
        assert 0<=both<=min(early,late) and early<=distinct<=n and late<=n
    def total(key):
        return sum(r[key] for r in near)
    assert total("shared_new_butterfly_species_pair_count")==j["all_shared_new_butterfly_species_pair_records"]==41
    assert total("distinct_first_record_year_species_pairs")==j["all_distinct_first_record_years"]==36
    assert total("upstream_re_record_before_downstream_first_year")==j["all_upstream_re_record_before_downstream"]==21
    assert total("downstream_re_record_within_five_years")==j["all_downstream_re_record"]==26
    assert total("previously_prior_host_matched_focal_links")==3
    assert j["other_new_pair_species_records_host_status_unclassified"]["n"]==38
    assert j["other_new_pair_species_records_host_status_unclassified"]["upstream_re_record_before_downstream"]==18
    assert j["no_statistical_test_due_to_postselection_clustered_pairs_host_unclassified"]
    assert j["raw_rows_published"] is False

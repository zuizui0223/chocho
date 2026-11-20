#!/usr/bin/env python3
"""Offline original-WPI1596, NEW GBIF datasetKey snapshot verifier.

No GBIF API calls, no original private source files, no DOI registration.
The current public snapshot does NOT retroactively attribute the older 2026
original speciesKey/WKT port botanical dataset.
"""
import csv,hashlib,json,uuid
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/geb_v05_gbif_source_attribution_20261120"
FILE=BASE/"NEW_GBIF_ORIGINAL1596_PORTS_DATASETKEY_INDEX_20261120.json"
PORTS=ROOT/"provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"
SOURCE_SHA="73aa16ce190aee301c790786acbee01beda3ca9b53b02b8226a7b4a762b628e7"
PORTS_BLOB="fdf73b11ddcb1644f46730123350ae52cc61cae6"

def main():
    raw=PORTS.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    assert blob==PORTS_BLOB, "Original port roster drift"
    src=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    selected={int(x["port_id"]) for x in src if x["query_complete"]=="True"}
    assert len(src)==3630 and len(selected)==1596
    data=FILE.read_bytes()
    assert hashlib.sha256(data).hexdigest()==SOURCE_SHA, "New snapshot bytes changed"
    q=json.loads(data)
    assert q["schema"]=="chocho_full1596_originalWPI_ports_NEW_2026_GBIF_botanical_datasetKey_index_v1"
    assert q["original_source_complete_1596_ports"]==1596
    assert q["original_source_all_WPI_candidates"]==3630
    assert q["original_WPI2019_coordinates_match_to_6dp"]
    assert q["original_pooled_botanical_speciesKey_cache_not_reproduced"]
    assert q["original_cached_query_used_WKT_buffer_nonidentical"]
    assert q["original_port_plant_source_datasetKey_citation_NOT_recovered"]
    assert q["read_only_NO_GBIF_DOI_minted"]
    assert q["new_GBIF_query_period"]=="2001,2024"
    rows=q["source_rows"]
    assert len(rows)==1596
    ids=[int(x["original_WPI2019_OBJECTID"]) for x in rows]
    assert set(ids)==selected and len(set(ids))==1596
    assert ids==sorted(ids)
    assert sum(x["new_2026_GBIF_index_count"]>0 for x in rows)==1461
    assert sum(x["new_2026_GBIF_index_count"]==0 for x in rows)==135
    tot=sum(x["new_2026_GBIF_index_count"] for x in rows)
    assert tot==3138675
    provider_counts=defaultdict(int)
    for r in rows:
        portcounts=r["current_GBIF_record_counts_by_datasetKey"]
        assert sum(x["count"] for x in portcounts)==r["new_2026_GBIF_index_count"]
        assert len(portcounts)==len({x["datasetKey"] for x in portcounts})
        for item in portcounts:
            uuid.UUID(item["datasetKey"])
            assert item["count"]>=0
            provider_counts[item["datasetKey"]]+=item["count"]
    assert len(provider_counts)==1482
    assert sum(provider_counts.values())==tot
    assert len(q["new_indexed_record_counts_by_datasetKey_WITH_OVERLAP"])==1482
    saved={z["datasetKey"]:z["sum_port_buffer_record_hits"]
           for z in q["new_indexed_record_counts_by_datasetKey_WITH_OVERLAP"]}
    assert saved==dict(provider_counts)
    receipt=q["new_snapshot_shards_acquired_at_distinct_times"]
    assert sorted(z["part"] for z in receipt)==[0,1,2,3]
    assert [z["records_per_port_buffer_sum"] for z in sorted(receipt,key=lambda z:z["part"])] == [745423,787905,820841,784506]
    for x in receipt:
        assert x["new_source_start_utc"].startswith("2026-11-20T")
        assert x["new_source_end_utc"].startswith("2026-11-20T")
    byrank=sorted(provider_counts.items(),key=lambda z:-z[1])
    assert byrank[0][1]==851148
    print("GBIF_NEW1596_WPI_PORT_SNAPSHOT_OFFLINE_SOURCE_SHA_AND_PROVIDER_SUMS_PASS",
      json.dumps({"ports":len(rows),"positive":1461,"zero":135,
      "unique_current_datasetKeys":len(provider_counts),
      "summed_overlapping_record_index":tot,
      "largest_provider_port_buffer_hits":byrank[0][1],
      "top10_source_port_buffer_hit_share":sum(z[1] for z in byrank[:10])/tot}))

if __name__=="__main__":main()

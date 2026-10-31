#!/usr/bin/env python3
"""Replay source-locked all23 port taxonomic coverage from downloaded original ZIP.

Observation richness is not ecological flora richness, host absence, actual
introductions, shipping pathways or butterfly establishment.
"""
from pathlib import Path
from collections import defaultdict
import argparse,csv,gzip,hashlib,io,json,zipfile

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
SOURCE_SHA="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
INNER="local_all23ports_specieskey_year.csv.gz"

def audit(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        if z.namelist()!=[INNER]:
            raise ValueError("Unexpected original artifact structure")
        compressed=z.read(INNER)
    if hashlib.sha256(compressed).hexdigest()!=SOURCE_SHA:
        raise ValueError("Original GBIF all23 speciesKey×port×year source SHA changed")
    rows=list(csv.DictReader(io.StringIO(gzip.decompress(compressed).decode("utf8"))))
    if len(rows)!=9788:raise ValueError("Source 9788 taxon×port×year rows changed")
    portsets=defaultdict(set)
    portcounts=defaultdict(int)
    islandsets=defaultdict(set)
    for row in rows:
        iid,pid,sp=row["island_id"],row["port_id"],row["gbif_speciesKey"]
        y,n=int(row["year"]),int(row["indexed_record_count"])
        if not 2001<=y<=2010 or n<=0:raise ValueError("Invalid source date or record count")
        portsets[iid,pid].add(sp)
        islandsets[iid].add(sp)
        portcounts[iid,pid]+=n
    saved=json.loads((HERE/"port_rank_taxonomic_coverage_20261031.json").read_text(encoding="utf8"))
    stats=saved["per_port_aggregated"]
    assert len(stats)==23
    for p in stats:
        key=p["island_id"],p["port_id"]
        assert len(portsets[key])==p["observed_taxon_keys"]
        assert portcounts[key]==p["indexed_taxon_identified_record_count"]
    for island in saved["island_coverage_comparison"]:
        assert len(islandsets[island["island_id"]])==island["across_port_species_key_union"]
        assert sum(x["island_id"]==island["island_id"] for x in stats)==island["source_wpi_port_count"]
    assert len(islandsets)==9
    print(json.dumps({"raw_sha_verified":True,"islands":len(islandsets),"source_rows":len(rows)}))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-zip",type=Path,required=True)
    a=p.parse_args()
    audit(a.source_zip)

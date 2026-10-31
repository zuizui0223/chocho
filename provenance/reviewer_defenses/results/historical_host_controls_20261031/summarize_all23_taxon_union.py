#!/usr/bin/env python3
"""Replay 23-port public botanical species-key union and year-first summaries.

Requires the downloaded Actions ZIP containing the GBIF speciesKey×port×year
facet results. A missing source is NOT a biological absence. No host edges.
"""
from pathlib import Path
from collections import defaultdict
import argparse,csv,gzip,hashlib,io,json,zipfile

HERE=Path(__file__).resolve().parent
SHA="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
MEMBER="local_all23ports_specieskey_year.csv.gz"
PAIRS=[("81","214"),("821","951"),("276","321"),("169","715"),
       ("715","1577"),("169","276"),("276","715")]

def compute(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        if z.namelist()!=[MEMBER]:raise ValueError("Unexpected private source ZIP members")
        raw=z.read(MEMBER)
    if hashlib.sha256(raw).hexdigest()!=SHA:
        raise ValueError("Original GBIF all-port facet gzip SHA mismatch")
    rows=list(csv.DictReader(io.StringIO(gzip.decompress(raw).decode("utf8"))))
    if len(rows)!=9788:raise ValueError("Wrong source speciesKey×port×year row count")
    early=defaultdict(dict)
    presence_years=defaultdict(set)
    ports=defaultdict(set)
    for r in rows:
        i=r["island_id"];k=r["gbif_speciesKey"];y=int(r["year"])
        if not 2001<=y<=2010 or int(r["indexed_record_count"])<=0:
            raise ValueError("Invalid source year or facet count")
        early[i][k]=min(y,early[i].get(k,y))
        presence_years[i].add(y);ports[i].add(r["port_id"])
    pairs=[]
    for a,b in PAIRS:
        x,y=early[a],early[b]
        shared=set(x)&set(y);union=set(x)|set(y)
        pairs.append(dict(island_a=a,island_b=b,
            shared_indexed_species_keys_pre2011=len(shared),
            union_indexed_species_keys_pre2011=len(union),
            jaccard=round(len(shared)/len(union),6),
            species_first_indexed_earlier_on_a=sum(x[k]<y[k] for k in shared),
            species_first_indexed_earlier_on_b=sum(y[k]<x[k] for k in shared),
            species_first_indexed_same_year=sum(x[k]==y[k] for k in shared)))
    return early,pairs

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-zip",required=True,type=Path)
    args=p.parse_args()
    early,pairs=compute(args.source_zip)
    saved=json.loads((HERE/"all23_botanical_taxon_union_20261031.json").read_text())
    assert saved["prespecified_priority_and_focal_pair_overlap"]==pairs
    assert {r["island_id"]:r["unique_indexed_species_keys_across_all_ports"]
            for r in saved["island_union_of_indexed_plant_species_keys"]}=={
        k:len(v) for k,v in early.items()}
    print(json.dumps({"source_verified":True,"islands":len(early),"pairs":len(pairs)}))
if __name__=="__main__":main()

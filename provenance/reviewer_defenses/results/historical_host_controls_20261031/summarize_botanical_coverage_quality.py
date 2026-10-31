#!/usr/bin/env python3
"""Verify private GBIF speciesKey×port×year artifact and derive public quality summary.

The ZIP has publicly indexed botanical species keys but excludes raw occurrence
locations and all restricted butterfly-host catalogue edges. No species non-
detection, historical alien status, shipping mechanism or establishment claim.
"""
import argparse,csv,gzip,io,zipfile,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
HERE=Path(__file__).resolve().parent
GZIP_SHA="26ad181ee2d021e3ed8d7ac5568c451562c1114e89c3543f07cee71cfd520add"

def audit(zip_path):
    with zipfile.ZipFile(zip_path) as archive:
        names=archive.namelist()
        if names!=["local_private_gbif_specieskey_year.csv.gz"]:
            raise ValueError(f"Unexpected source ZIP members {names}")
        raw=archive.read(names[0])
    if hashlib.sha256(raw).hexdigest()!=GZIP_SHA:
        raise ValueError("Source GBIF species-year gzip SHA256 differs from first acquisition")
    rows=list(csv.DictReader(io.StringIO(gzip.decompress(raw).decode("utf8"))))
    if len(rows)!=4843:raise ValueError("Expected fixed dated plant facet row count 4843")
    plants=defaultdict(Counter);years=defaultdict(Counter)
    for r in rows:
        island=r["island_id"];key=r["gbif_speciesKey"]
        year=int(r["year"]);amount=int(r["indexed_records"])
        if not 2001<=year<=2010 or amount<1:
            raise ValueError("Bad original botanical year/facet count")
        plants[island][key]+=amount
        years[island][year]+=amount
    if len(plants)!=9:raise ValueError("Different source island count")
    return plants,years

def compute(plants,years):
    site=[]
    for i,counts in sorted(plants.items(),key=lambda x:int(x[0])):
        total=sum(counts.values());ordered=sorted(counts.values(),reverse=True)
        p=[x/total for x in ordered]
        site.append(dict(
            island_id=i,unique_species_keys=len(counts),records_with_species_key=total,
            top_species_records=ordered[0],top_species_fraction=round(ordered[0]/total,6),
            inverse_simpson_effective_species=round(1/sum(x*x for x in p),3),
            exp_shannon_effective_species=round(math.exp(-sum(x*math.log(x) for x in p)),3),
            distinct_recorded_years=len(years[i]),max_record_year=max(years[i],key=years[i].get)))
    selected=[("81","214"),("821","951"),("276","321"),
        ("169","715"),("715","1577"),("169","276"),("276","715")]
    pairs=[]
    for a,b in selected:
        A=set(plants[a]);B=set(plants[b]);shared=len(A&B);union=len(A|B)
        pairs.append(dict(island_a=a,island_b=b,
            shared_early_botanical_species_keys=shared,jaccard=round(shared/union,6),
            union_species_keys=union,
            observations_are_historical_port_search_not_verified_larval_hosts=True))
    return site,pairs

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-zip",type=Path,required=True)
    args=p.parse_args()
    plants,years=audit(args.source_zip)
    sites,pairs=compute(plants,years)
    saved=json.loads((HERE/"botanical_coverage_quality_20261031.json").read_text(encoding="utf8"))
    assert saved["sites"]==sites and saved["pairs"]==pairs
    assert saved["source_species_year_facet_rows"]==4843
    assert saved["unique_species_keys_across_nine_ports"]==len(set().union(*[set(x) for x in plants.values()]))
    assert sum(x["records_with_species_key"] for x in sites)==14081
    print(json.dumps({"validated":True,"sites":len(sites),"speciesYearRows":4843,"uniqueSpecies":1447}))
if __name__=="__main__":main()

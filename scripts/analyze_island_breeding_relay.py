"""Verify breeding-stage witnesses against the frozen index and count ordered links."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from acquire_island_breeding_relay import OUT,ROOT,SOURCE


def verify(workspace):
    base=workspace/"work/barriers"
    hits=pd.read_csv(OUT/"stage_search_hits.csv.gz")
    pairs=pd.read_csv(OUT/"candidate_links.csv")
    geography=pd.read_parquet(base/"original_coordinate_geography.parquet")
    ids=set(pairs.island)|set(pairs.destination)
    geography=geography[geography.landmass_id.isin({str(i) for i in ids})].copy()
    geography["island"]=geography.landmass_id.astype(int)
    xy=["decimalLatitude","decimalLongitude"]
    needed_keys=set(pairs.key)
    dest=pairs[["key","destination","first_year"]].drop_duplicates().rename(columns={"destination":"island","first_year":"year"})
    witness=[];dest_records=[]
    for path in sorted((workspace/"work/temporal/occurrence_parts").glob("*.parquet")):
        chunk=pd.read_parquet(path,columns=["gbifID","speciesKey","year","eventDate",*xy],
                              filters=[("speciesKey","in",list(needed_keys)),("year",">=",2011)])
        chunk=chunk.rename(columns={"gbifID":"gbif_id","speciesKey":"key"}).merge(geography[xy+["island"]],on=xy,validate="many_to_one")
        witness.append(chunk[chunk.gbif_id.isin(hits.gbif_id)])
        dest_records.append(chunk.merge(dest,on=["key","island","year"]))
    original=pd.concat(witness).drop_duplicates(["gbif_id","key","year","island"])
    verified=hits.merge(original[["gbif_id","key","year","island"]],on=["gbif_id","key","year","island"],validate="many_to_one")
    verified["larva"]=verified["query"].eq("Larva") & verified.stage.fillna("").str.lower().isin(["larva","larvae","caterpillar"])
    verified["egg"]=verified["query"].eq("Egg") & verified.stage.fillna("").str.lower().isin(["egg","eggs"])
    verified=verified.groupby(["gbif_id","key","year","island"],as_index=False).agg(larva=("larva","max"),egg=("egg","max"),laying=("laying","max"),dataset_key=("dataset_key","first"),uncertainty_m=("uncertainty_m","first"))
    verified.to_csv(OUT/"verified_stage_records.csv",index=False)
    dest_records=pd.concat(dest_records).sort_values(["eventDate","gbif_id"]).drop_duplicates(["key","island","year"]).rename(columns={"island":"destination","year":"first_year","gbif_id":"destination_gbif_id"})
    assert len(dest_records)==len(dest)
    dest_records[["key","destination","first_year","destination_gbif_id","eventDate"]].to_csv(OUT/"destination_witnesses.csv",index=False)
    flora=pd.read_csv(SOURCE/"port_host_status.csv").query("status in ['native','introduced']")
    edges=pd.read_parquet(base/"all225_named_host_edges_local.parquet").drop_duplicates()
    plants=pairs[["island","species"]].drop_duplicates().merge(edges).merge(flora,on=["island","host_name"])
    plants.to_csv(OUT/"compatible_port_hosts.csv",index=False)


def count(g):
    return dict(middle_islands=g.island.nunique(),destination_islands=g.destination.nunique(),species=g.key.nunique(),
                links=len(g.drop_duplicates(["island","destination","key"])),destination_firsts=len(g.drop_duplicates(["destination","key"])),
                middle_species=len(g.drop_duplicates(["island","key"])))


def analyze():
    pairs=pd.read_csv(OUT/"candidate_links.csv")
    stage=pd.read_csv(OUT/"verified_stage_records.csv").rename(columns={"year":"stage_year","gbif_id":"stage_gbif_id"})
    dest=pd.read_csv(OUT/"destination_witnesses.csv")
    z=pairs.merge(stage,on=["key","island"]).merge(dest,on=["key","destination","first_year"],validate="many_to_one")
    z=z[z.stage_year.lt(z.first_year)].copy()
    z["direct_stage"]=z.larva|z.laying
    assert z.stage_year.gt(2010).all() and (z.island!=z.destination).all()
    z.to_csv(OUT/"ordered_witnesses.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    links=z[z.direct_stage].sort_values(["km","species","stage_year","stage_gbif_id"]).drop_duplicates(["key","island","destination"])
    links.to_csv(OUT/"ordered_links.csv",index=False)
    checks=[]
    for km in [100,250,500]:
        for kind,g in [("candidates",pairs),("larva or explicit laying",z[z.direct_stage]),("larva",z[z.larva]),("explicit laying",z[z.laying]),
                       ("including egg records",z[z.direct_stage|z.egg]),("nearest destination only",z[z.direct_stage&z.neighbor_rank.eq(1)]),
                       ("introduced host confirmed",z[z.direct_stage&z.introduced.gt(0)]),
                       ("stage uncertainty <=1 km",z[z.direct_stage&z.uncertainty_m.gt(0)&z.uncertainty_m.le(1000)])]:
            checks.append(dict(scope=kind,km=km,**count(g[g.km.le(km)])))
    counts=pd.DataFrame(checks);counts.to_csv(OUT/"counts.csv",index=False)
    print(counts.to_string(index=False));print(links[links.km.le(250)][["species","island","destination","km","stage_year","first_year","stage_gbif_id","destination_gbif_id"]].to_string(index=False))
    return links,counts


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--workspace",type=Path)
    w=ap.parse_args().workspace
    if w:verify(w)
    analyze()

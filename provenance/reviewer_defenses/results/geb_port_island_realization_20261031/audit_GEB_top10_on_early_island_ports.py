#!/usr/bin/env python3
"""Join GEB's actual top larval-resource contributors to PRE-2011 GBIF port plants.

Source-locked 2019 World Port Index ports, independently queried 2001–2010
GBIF taxon-by-port-by-year facets. No plant's WCVP native/introduced status
at port in 2001–2010 is demonstrated; no shipping mode, original plant
introduction year, feeding, successful butterfly arrival or colonization.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import argparse, csv,gzip,hashlib,io,json
ROOT=Path(__file__).resolve().parents[4]
SOURCE=Path("provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json")
CONTRIB_SHA="201a762de94b9ecf42a766c8cb5d8d1739a2ca3c"
TAX_SHA="d8d21e6686e3f069d4374d4be1c4e47744788f9244c46a1bbbbb9015e5cd7fb9"
PLANT23_SHA="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
PLANT15_SHA="44dc11ad3a0b73aef50e92ca28f49db3927f217de9496398f1b9ba9687a2bf15"

def sha(data):return hashlib.sha256(data).hexdigest()

def content(root):
 b=(root/SOURCE).read_bytes()
 name=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()
 if name!=CONTRIB_SHA:raise ValueError("GEB original contributor source Git blob changed")
 d=json.loads(b)
 assert len(d["plant_species"]["top10"])==10
 assert abs(sum(float(x["share"]) for x in d["plant_species"]["top10"])-d["plant_species"]["top_k_share"]["10"])<1e-12
 return d["plant_species"]["top10"]

def original_csv(path,expected,name):
 b=path.read_bytes()
 if sha(b)!=expected:raise ValueError("Author original "+name+" bytes differ from source SHA256")
 return b

def audit(root,taxonomy,plant23,plant15):
 top=content(root)
 tax=json.loads(gzip.decompress(original_csv(taxonomy,TAX_SHA,"taxonomic crosswalk")))
 assert tax["schema"]=="chocho_accepted_original_host_species_gbif_crosswalk_v1"
 d=tax["name_resolution"]
 rows=[]
 for src,expected,old in [(plant23,PLANT23_SHA,True),(plant15,PLANT15_SHA,False)]:
  data=gzip.decompress(original_csv(src,expected,"dated plant facets")).decode("utf8")
  source=list(csv.DictReader(io.StringIO(data)))
  if len(source)!=(9788 if old else 327):raise ValueError("GBIF port species-year source cohort changed")
  for x in source:
   rows.append((x["island_id"],x["port_id"],int(x["year"]),
      x["gbif_speciesKey"] if old else x["speciesKey"],
      int(x["indexed_record_count"] if old else x["indexed_records"])))
 out=[]
 for p in top:
  name=p["accepted_name"];share=float(p["share"])
  meta=d.get(name)
  if not meta or meta.get("gate")!="EXACT_ACCEPTED" or not meta.get("gbif_speciesKey"):
   raise ValueError("GEB primary host requires ambiguous/unresolved GBIF taxonomy: "+name)
  key=meta["gbif_speciesKey"]
  records=[x for x in rows if x[3]==key]
  out.append({"plant_name":name,"GBIF_accepted_speciesKey":str(key),
    "GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant":share,
    "source_strict_accepted_taxon_match":True,
    "observed_island_ids":sorted(set(r[0] for r in records),key=int),
    "n_study_islands_with_pre2011_record":len(set(r[0] for r in records)),
    "n_port_buffers_with_pre2011_record":len(set(r[1] for r in records)),
    "first_indexed_plant_year":min((r[2] for r in records),default=None),
    "sum_indexed_records_in_potentially_overlapping_WPI_5km_port_buffers":sum(r[4] for r in records),
    "zero_means_not_indexed_in_this_sample_NOT_absence":not records})
 total=sum(x["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"] for x in out)
 observed=sum(bool(x["n_study_islands_with_pre2011_record"]) for x in out)
 observed_share=sum(x["GEB_fraction_of_added_species_region_opportunity_attributed_to_this_plant"] for x in out if x["n_study_islands_with_pre2011_record"])
 assert observed==8 and abs(total-.2510105136113853)<1e-12
 return {
  "schema":"chocho_GEB_top10_host_contributors_on_pre2011_island_ports_v1",
  "status":"SOURCE_VERIFIED_SHARED_PLANT_TAXA_NOT_PORT_REALIZATION_OF_GLOBAL_GAIN",
  "source_GEB_contributor_git_blob":CONTRIB_SHA,
  "taxon_accepted_name_crosswalk_gzip_sha256":TAX_SHA,
  "historical_23port_plant_gzip_sha256":PLANT23_SHA,
  "historical_15port_plant_gzip_sha256":PLANT15_SHA,
  "original_GEB_total_added_species_WGSRPD3_units":14553,
  "study_islands":17,"study_port_proxies":38,"plant_record_years":"2001–2010",
  "GEB_top10_species":len(top),"GEB_top10_fraction_of_added_units":total,
  "top10_species_with_plant_index_at_study_island_ports":observed,
  "GEB_top10_species_WITH_any_port_index_total_fraction_of_GLOBAL_added_units":observed_share,
  "source_contributor_plants":out,
  "no_2001_2010_site_specific_alien_status":True,
  "no_host_growth_density_or_local_butterfly_feeding_evidence":True,
  "no_shipping_route_or_causal_introduction":True,
  "top10_contribution_share_not_fraction_of_global_55pct_realized_on_islands":True,
  "source_GBIF_index_reporting_heterogeneous_and_WPI_retroactive":True
 }
if __name__=="__main__":
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--root",type=Path,default=ROOT)
 p.add_argument("--taxonomy",type=Path,required=True)
 p.add_argument("--plant23",type=Path,required=True)
 p.add_argument("--plant15",type=Path,required=True)
 p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 r=audit(a.root,a.taxonomy,a.plant23,a.plant15)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
 print("GEB_TOP10_HISTORIC_PORTS",json.dumps({"species":len(r["source_contributor_plants"]),
   "indexed":r["top10_species_with_plant_index_at_study_island_ports"],
   "share_of_GEB_global_added_units":r["GEB_top10_species_WITH_any_port_index_total_fraction_of_GLOBAL_added_units"]}))

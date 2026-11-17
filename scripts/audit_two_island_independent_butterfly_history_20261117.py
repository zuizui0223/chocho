#!/usr/bin/env python3
"""Extend independent island-history audit to Mauritius (1993) without changing GBIF.

The historical record-list is a heterogeneous-source *known presence* baseline,
not unbiased ecological occupancy. Butterfly origin/endemism stay unassessed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[1]
LEGACY=ROOT/"scripts/audit_independent_island_historical_butterfly_presence_20261117.py"
NEW=ROOT/"results/butterfly_GBIF_effort_origin_gate_20261117/MAURITIUS_1993_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv"
SHA="e06fb0032b8dd2a9bd6344e2e88db0f201f3a4ee"
BHL="https://www.biodiversitylibrary.org/part/332210"
DOI="10.5962/p.332210"
SOURCE_BEFORE=1993
EXPECTED={"Zizula hylax":1922244,"Lampides boeticus":1924713,"Danaus chrysippus":7642610}
UNCERTAIN={"Zizina otis":"taxonomy_concept_of_1993_Zizina_antanossa_vs_2008_Zizina_otis_requires_exact_GBIF_match",
           "Zerynthia polyxena":"out_of_historical_Mauritius_checklist_and_European_range; inspect_original_GBIF_occurrence_id_taxon_coordinates"}

def git_blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()

def load_legacy():
    spec=importlib.util.spec_from_file_location("island_hist_2007",LEGACY)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def new_hist():
    if git_blob(NEW)!=SHA:
        raise ValueError("Mauritius source table SHA drift")
    with NEW.open(newline="",encoding="utf-8-sig") as f:
        r=list(csv.DictReader(f))
    if len(r)!=3:
        raise ValueError("Mauritius exact species evidence must be 3 rows")
    return r

def analyze():
    m=load_legacy()
    rows=m.read_pinned("gbif")
    geo=m.read_pinned("islands")
    existing=m.read_pinned("history")
    base,annotated=m.audit(rows,geo,existing)
    if len(annotated)!=317 or base["Phu_Quoc"]["documented_2006_or_2007_butterflies"]!=12:
        raise ValueError("Legacy historical checklist source drift")
    islands={int(z["landmass_id"]):z for z in geo}
    island=islands[218]
    if island["example_port_name"]!="PORT LOUIS" or island["example_port_country_code"]!="MU":
        raise ValueError("Mauritius 218 not Port Louis original WPI anchor")
    if not 57<float(island["representative_longitude"])<58 or not -21<float(island["representative_latitude"])<-20:
        raise ValueError("Mauritius geographic original representative mismatch")
    by_key={(int(x["speciesKey"]),int(x["landmass_id"])):x for x in annotated}
    # Source types must not be conflated with native/establishment/endemism.
    audit=[]
    seen=set()
    for record in new_hist():
        name=record["species"].strip()
        key=(int(record["speciesKey"]),int(record["landmass_id"]))
        if name not in EXPECTED or EXPECTED[name]!=key[0] or key[1]!=218:
            raise ValueError("Unapproved taxon or island in Mauritius source table")
        if key in seen or key not in by_key:
            raise ValueError("Historical evidence key duplicate or not originally first-indexed")
        seen.add(key)
        x=by_key[key]
        if x["species"]!=name or x["earliest_GBIF_index_period"]!="2011-2017":
            raise ValueError("Mauritius source does not predate 2011-2017 GBIF first index")
        if (record["island"],int(record["publication_year_upper_bound"]))!=("Mauritius",1993):
            raise ValueError("Mauritius publication year or locality drift")
        if int(record["publication_year_upper_bound"])>=2011 or record["name_match_strength"]!="source_exact_species":
            raise ValueError("Wrong evidence age or not exact name match")
        if record["source_url"]!=BHL or record["source_doi"]!=DOI:
            raise ValueError("Source-URL or DOI provenance drift")
        if (record["origin_status"],record["establishment_status"],record["endemicity_status"])!=("unknown","unknown","unassessed"):
            raise ValueError("Historical checklist is not evidence of native origin, breeding, or endemicity")
        if int(record["printed_page"]) not in (147,154):
            raise ValueError("Primary published page source drift")
        x["historical_field_presence"]="source_verified_by_1993"
        x["last_possible_historical_presence_year"]=1993
        x["source_url"]=BHL
        audit.append(dict(species=name,speciesKey=key[0],landmass_id=218,
                          earliest_independent_publication_upper_bound=1993,
                          printed_page=int(record["printed_page"]),
                          GBIF_first_index_period=x["earliest_GBIF_index_period"],
                          minimum_index_window_gap_years=2011-1993,
                          independent_source_doi=DOI,source_url=BHL,
                          native_origin_not_assessed=True,establishment_not_assessed=True,
                          endemicity_not_assessed=True))
    if set(EXPECTED)!=set(x["species"] for x in audit):
        raise ValueError("Missing 1993 exact-species evidence")
    allMaur=[x for x in annotated if x["landmass_id"]==218]
    if len(allMaur)!=5:
        raise ValueError("Original five Mauritius first-index taxon cohort drift")
    if {x["species"] for x in allMaur}!=set(EXPECTED)|set(UNCERTAIN):
        raise ValueError("Mauritius original taxon names drift")
    for x in allMaur:
        if x["species"] in UNCERTAIN:
            x["historical_field_presence"]="unresolved_taxonomic_or_geospatial_quality"
            x["source_url"]=BHL
            x["audit_note"]=UNCERTAIN[x["species"]]
    # Crosswalk sets use only exact historically corroborated species.
    early={i:set() for i in islands}
    late={i:set() for i in islands}
    names=set()
    for r in rows:
        sk=int(r["speciesKey"]);i=int(r["landmass_id"])
        names.add(r["species"])
        a,b,c=[int(r[z]) for z in (
            "pre2001_2010_indexed_records",
            "future2011_2017_indexed_records","future2018_2022_indexed_records")]
        if a: early[i].add(sk)
        if a+b+c: late[i].add(sk)
    phu_aug={i:set(z) for i,z in early.items()}
    for h in existing:
        phu_aug[int(h["landmass_id"])].add(int(h["speciesKey"]))
    both_aug={i:set(z) for i,z in phu_aug.items()}
    for key in seen:
        both_aug[key[1]].add(key[0])
    ids=sorted(islands)
    j0,p0=m.mean_pair_jaccard(ids,early)
    jp,pp=m.mean_pair_jaccard(ids,phu_aug)
    jm,pm=m.mean_pair_jaccard(ids,both_aug)
    jl,pl=m.mean_pair_jaccard(ids,late)
    if not all(p==1275 for p in (p0,pp,pm,pl)):
        raise ValueError("Changed island pair sample")
    confirmed=sum(x["historical_field_presence"].startswith("source_verified_by_") for x in annotated)
    if confirmed!=15:
        raise ValueError("Source-supported historical confirmations not 15")
    if sum(x["butterfly_origin"]!="unknown" or x["butterfly_establishment"]!="unknown" or x["butterfly_endemicity"]!="unassessed" for x in annotated):
        raise ValueError("Independent historical presence must not be treated as status determination")
    native_examples=[
        {"historical_name":"Euploea euphon","type":"1993 explicitly calls endemic Mauritius species"},
        {"historical_name":"Cyclyrius mandersi","type":"1993 source says confined to Mauritius; taxonomy needs current reconciliation"}
    ]
    if any(x["historical_name"] in names for x in native_examples):
        raise ValueError("Named historical source-native island species already represented in selected 225")
    report={
      "schema":"chocho_two_island_source_verified_historical_225x51_consumer_GBIF_record_lag_v2",
      "status":"HISTORICAL_OBSERVATION_COVERAGE_AUDIT_NO_TRUE_COLONIZATION_OR_NATIVE_RANGE_INFERENCE",
      "frozen_source_git_blobs":{**base["source_git_blob_shas"],"Mauritius_1993_crosswalk":SHA},
      "cohort":{"butterfly_species":len(names),"island_landmasses":len(islands),
                "species_island_cells":len(rows),"first_GBIF_index_cells_2011_2022":len(annotated),
                "historically_confirmed_before_2011":confirmed,
                "not_confirmed_by_these_historical_sources":len(annotated)-confirmed,
                "Phu_Quoc_confirmed":12,"Mauritius_confirmed":3,
                "Mauritius_unresolved_identity_or_taxon_cases":2,
                "independent_islands_with_historical_matches":2},
      "Mauritius":{
          "id":218,"WPI_port":"PORT LOUIS","source_year_conservative":SOURCE_BEFORE,
          "source_title":"The Butterflies of Mauritius",
          "source_authors":"P. M. H. Davis and M. J. C. Barnes",
          "source_bibliographic_date":"1993-12-01 per BHL metadata (paper imprint may show 1991)",
          "DOI":DOI,"source_url":BHL,
          "GBIF_first_indexes_2011_2022":len(allMaur),
          "confirmed_historical_same_name":len(audit),
          "source_exact_matches":audit,
          "other_two_taxa_outcome":UNCERTAIN,
          "source_reports_historical_insular_taxa_missing_from_225":native_examples,
          "note":"The 1993 source has a species-level island checklist; the new source-verified 3 are exact binomials, without asserting native origin or breeding."
      },
      "multi_source_same_51_island_pairs":{
          "pairs":1275,"early_GBIF_indexed_cells":sum(map(len,early.values())),
          "early_augmented_with_PhuQuoc_2006_2007":sum(map(len,phu_aug.values())),
          "early_augmented_with_PhuQuoc_and_Mauritius_1993":sum(map(len,both_aug.values())),
          "cumulative_GBIF_2022_cells":sum(map(len,late.values())),
          "Jaccard_early_GBIF_only":j0,
          "Jaccard_early_plus_PhuQuoc":jp,
          "Jaccard_early_plus_both_historical_islands":jm,
          "Jaccard_cumulative_GBIF_2022":jl,
          "delta_early_GBIF_to_late_GBIF":jl-j0,
          "delta_early_PhuQuoc_augmented_to_late":jl-jp,
          "delta_early_two_historical_island_augmented_to_late":jl-jm,
          "change_in_delta_from_fifteen_historical_rows":j0-jm,
          "boundary":"Heterogeneous historical checklist additions; not time-matched ecological survey completeness or true occurrence."
      },
      "limits":["This selected 15/317 is NOT the fraction of all butterfly GBIF first indices reflecting overlooked prior presence.",
                "The two chosen islands and single sources per island are not independent random replicate sites.",
                "Historical absence from a checklist is NOT proven absence; keep both unidentified Mauritius cases unresolved.",
                "Zizina otis may share a taxon concept with historical Zizina antanossa under a 2008 molecular revision; exact GBIF accepted taxon mapping must be checked.",
                "The lone Zerynthia polyxena Mauritius GBIF record is geographically anomalous for a European butterfly; requires raw occurrence verification before deletion.",
                "Species omitted from the 225 HOSTS-based candidate panel cannot establish true endemic capture fractions.",
                "Do not infer established, native, introduced, vagrant, or endemic state from GBIF first records or these historical species lists."]
    }
    return report,annotated

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=Path("_butterfly_two_island_history"))
    a=p.parse_args()
    report,rows=analyze()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"TWO_ISLAND_HISTORICAL_GBIF_GATE_20261117.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    with (a.out_dir/"ALL317_FIRST_GBIF_INDEX_HISTORICAL_2_ISLAND_AUDIT.csv").open("w",newline="",encoding="utf-8") as f:
        cols=sorted({k for x in rows for k in x})
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader()
        w.writerows(sorted(rows,key=lambda x:(x["landmass_id"],x["speciesKey"])))
    print(json.dumps({"cohort":report["cohort"],"Mauritius":{k:v for k,v in report["Mauritius"].items() if k!="source_exact_matches"},
                      "pair_index":report["multi_source_same_51_island_pairs"]},indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()

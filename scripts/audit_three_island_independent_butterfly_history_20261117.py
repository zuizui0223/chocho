#!/usr/bin/env python3
"""Three-island historical source augmentation of 225x51 GBIF butterfly indices.

Adds strictly named pre-2011 Cape Breton taxa from 1948/1984 field literature,
including one independent 1983 wild larval-stage observation on host plants.
Historic specimens are not native origin, endemicity, or stable establishment.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"scripts/audit_two_island_independent_butterfly_history_20261117.py"
CAPE=ROOT/"results/butterfly_GBIF_effort_origin_gate_20261117/CAPE_BRETON_1948_1984_INDEPENDENT_FIELD_MATCHES_20261117.csv"
PIN="7bf77344e0be697ff4a2497a36b91a3fd577d78c"
EXPECTED={"Erynnis icelus":5977250,"Glaucopsyche lygdamus":1924106,
          "Lycaena phlaeas":1929697,"Colias philodice":5137672}

def blob(p):
    raw=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def load_module(path):
    spec=importlib.util.spec_from_file_location("butterfly_three_island_base",path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def caperows():
    if blob(CAPE)!=PIN:raise ValueError("Cape Breton exact source CSV blob drift")
    with CAPE.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    if len(rows)!=4:raise ValueError("Expected four exact island source rows, no extras")
    return rows

def analyze():
    tw=load_module(BASE)
    prior,previous=tw.analyze()
    legacy=tw.load_legacy()
    original=legacy.read_pinned("gbif")
    geography=legacy.read_pinned("islands")
    existing_phu=legacy.read_pinned("history")
    old_maur=tw.new_hist()
    cape=caperows()
    first,all_species=legacy.pre2011_first_indices(original)
    landmasses={int(x["landmass_id"]):x for x in geography}
    area=landmasses.get(79)
    if area is None or area["example_port_name"]!="SYDNEY" or area["example_port_country_code"]!="CA":
        raise ValueError("Original Cape Breton study port geography changed")
    if not -61<float(area["representative_longitude"])<-60:
        raise ValueError("Cape Breton longitude doesn't match original source")
    targets=[(k,z) for k,z in first.items() if k[1]==79]
    if len(targets)!=17:raise ValueError("Original late first-index Cape Breton cohort drift")
    output={(int(x["speciesKey"]),int(x["landmass_id"])):dict(x) for x in previous}
    sourceKeys=set()
    life=[]
    for e in cape:
        name=e["species"].strip()
        key=(int(e["speciesKey"]),int(e["landmass_id"]))
        if name not in EXPECTED or EXPECTED[name]!=key[0] or key[1]!=79:
            raise ValueError("Unapproved Cape Breton species/island id")
        if key in sourceKeys or key not in output:
            raise ValueError("Duplicated or absent original GBIF first-index taxon")
        sourceKeys.add(key)
        x=output[key]
        if x["species"]!=name or x["earliest_GBIF_index_period"]!="2011-2017":
            raise ValueError("Mismatch original later GBIF record")
        yr=int(e["observation_year_upper_bound"])
        pub=int(e["publication_year_upper_bound"])
        if yr>1983 or pub>1984 or pub<yr:
            raise ValueError("Historical year is not source-backed prior to 2011")
        if e["source_url"] not in (
          "https://images.peabody.yale.edu/lepsoc/nls/1980s/1984/1984_v26_n2.pdf",
          "https://images.peabody.yale.edu/lepsoc/jls/1940s/1948/1948-2%28S%29-SEASON_SUMMARY.pdf"):
            raise ValueError("Source URL not from checked historical original")
        if not e["exact_island_place"] or not e["name_match_type"]:
            raise ValueError("No explicit original island locality or name resolution")
        if (e["origin_status"],e["establishment_status"],e["endemicity_status"])!=("unknown","unknown","unassessed"):
            raise ValueError("Historic source cannot automatically decide butterfly origin, stability or endemism")
        life_stage=e["life_stage"]
        if name=="Glaucopsyche lygdamus":
            if life_stage!="wild_larvae" or e["named_larval_host"]!="Vicia cracca; Lathyrus japonicus":
                raise ValueError("The one original larvae+host observation changed")
            reproduction="local_wild_larvae_observed_in_1983"
        else:
            if life_stage!="adult" or e["named_larval_host"]:
                raise ValueError("No unsupported larval use allowed for this taxon")
            reproduction="unknown_no_larval_record_in_source"
        x["historical_field_presence"]="source_verified_by_"+str(pub)
        x["last_possible_historical_presence_year"]=yr
        x["source_url"]=e["source_url"]
        x["independent_life_stage_evidence"]=reproduction
        x["source_confidence_note"]=e["name_match_type"]
        life.append({"species":name,"speciesKey":key[0],"landmass_id":79,
           "historical_observation_year_upper_bound":yr,
           "historical_publication_year_upper_bound":pub,
           "GBIF_first_index_window":x["earliest_GBIF_index_period"],
           "minimum_lag_start":2011-yr,
           "evidence_type":e["evidence_type"],
           "life_stage":life_stage,
           "named_larval_host":e["named_larval_host"],
           "source":e["source_url"],
           "source_section":e["source_section"],
           "matched_name_confidence":e["name_match_type"],
           "butterfly_origin":"unknown",
           "stable_establishment":"unknown",
           "single_island_endemism":"unassessed"})
    if set(EXPECTED)!=set(x["species"] for x in life):
        raise ValueError("Cape Breton source candidates missing")
    known=sum(str(x["historical_field_presence"]).startswith("source_verified_by_") for x in output.values())
    if known!=19:raise ValueError("Expected 12+3+4 exact historical source support")
    if any(x["butterfly_origin"]!="unknown" or x["butterfly_establishment"]!="unknown"
           or x["butterfly_endemicity"]!="unassessed" for x in output.values()):
        raise ValueError("Evidence leaked into native/established/endemic status")

    old={i:set() for i in landmasses}
    late={i:set() for i in landmasses}
    for row in original:
        i=int(row["landmass_id"]);s=int(row["speciesKey"])
        a=int(row["pre2001_2010_indexed_records"])
        b=int(row["future2011_2017_indexed_records"])
        c=int(row["future2018_2022_indexed_records"])
        if a:old[i].add(s)
        if a+b+c:late[i].add(s)
    source_aug={i:set(x) for i,x in old.items()}
    for r in existing_phu:source_aug[int(r["landmass_id"])].add(int(r["speciesKey"]))
    for r in old_maur:source_aug[int(r["landmass_id"])].add(int(r["speciesKey"]))
    old_three={i:set(x) for i,x in source_aug.items()}
    for s,i in sourceKeys:old_three[i].add(s)
    ids=sorted(landmasses)
    j0,p0=legacy.mean_pair_jaccard(ids,old)
    j2,p2=legacy.mean_pair_jaccard(ids,source_aug)
    j3,p3=legacy.mean_pair_jaccard(ids,old_three)
    jl,pL=legacy.mean_pair_jaccard(ids,late)
    if not (p0==p2==p3==pL==1275):raise ValueError("Original 51 x 1275 pair cohort drift")
    if len(output)!=317:raise ValueError("Original 317 late first-index risk universe drift")
    return {
      "schema":"chocho_225x51_three_island_1948_1984_1993_2006_historical_butterfly_source_gate_v1",
      "status":"SOURCE_LOCKED_HISTORICAL_PRESENCE_LARVAE_NOT_TRUE_NATIVE_OR_COLONIZATION",
      "frozen_GBIF_225x51_blob":legacy.PINS["gbif"][1],
      "frozen_Cape_Breton_evidence_git_blob":PIN,
      "previous_two_islands_Mauritius_Phu_Quoc":prior["cohort"],
      "cohort":{"butterflies":225,"landmass_groups":51,"species_island_cells":11475,
                "later_GBIF_first_index_cells":317,
                "historically_confirmed_on_three_selected_islands":19,
                "not_confirmed_by_these_historical_sources":298,
                "source_audited_early_first_indices_by_landmass":{"471":12,"218":3,"79":4},
                "Cape_Breton_first_indices_all":17,
                "Cape_Breton_old_documented_same_butterfly":4,
                "historical_field_larvae_pre2011_species_island_cells":1},
      "Cape_Breton":{"landmass_id":79,"island":"Cape Breton Island","WPI_port":"SYDNEY",
                     "source_record_rows":sorted(life,key=lambda x:x["species"]),
                     "larval_record":"Glaucopsyche lygdamus on Vicia cracca and Lathyrus japonicus at distinct sites in 1983",
                     "larval_evidence_boundary":"Wild larvae at one time support local larval development, not long-term establishment or native origin",
                     "taxon_name_boundary":"Lycaena phlaeas appears as abbreviated genus in 1984; 2010 exact-name local atlas independently corroborates; Erynnis icelus needs l/I PDF text normalization"},
      "same_1275_pairs_heterogeneous_source_augmented_known_presence":{
         "original_early_GBIF_cells":sum(len(x) for x in old.values()),
         "after_prior_two_historical_islands_cells":sum(len(x) for x in source_aug.values()),
         "after_three_historical_islands_cells":sum(len(x) for x in old_three.values()),
         "late_cumulative_GBIF_cells":sum(len(x) for x in late.values()),
         "Jaccard_early_GBIF_only":j0,
         "Jaccard_early_two_islands_augmented":j2,
         "Jaccard_early_three_islands_augmented":j3,
         "Jaccard_2022_late":jl,
         "delta_GBIF_only":jl-j0,
         "delta_two_islands":jl-j2,
         "delta_three_islands":jl-j3,
         "delta_change_from_19_historical_matches":j0-j3,
         "estimand":"historically augmented documentary record-list comparison, not adjusted ecological species turnover"},
      "limits":["Three islands were selected and audited from available historical reports; 19/317 is NOT a representative global initial-index error rate",
                "The 1984 and 1948 newsletters are independent of the GBIF 2001-2022 record window but the same field sites can appear in subsequent compilations",
                "Wild larvae can provide stronger evidence for local life-stage use than adult sightings; no multiyear population persistence is concluded",
                "Do not auto-promote native/introduced/established/endemic labels from historical occurrence.",
                "The original 225 HOSTS-linked consumer species are not an exhaustive butterfly-endemic fauna; other endemic taxa may not be represented."]
    },list(output.values())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,default=Path("_butterfly_three_island_histories"))
    args=ap.parse_args()
    d,r=analyze()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"THREE_ISLAND_SOURCE_AUDIT_20261117.json").write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    columns=sorted({k for x in r for k in x})
    with (args.out_dir/"ALL317_INDEXED_CELLS_HISTORICAL_GATES_3ISLAND.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=columns)
        writer.writeheader()
        writer.writerows(sorted(r,key=lambda x:(x["landmass_id"],x["speciesKey"])))
    print(json.dumps({"cohort":d["cohort"],"Jaccard":d["same_1275_pairs_heterogeneous_source_augmented_known_presence"],
       "Cape_Breton":d["Cape_Breton"]["source_record_rows"]},ensure_ascii=False,indent=2))
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Three-gate butterfly island case/trait feasibility census.

Outputs are a *source provenance screen*. Not an establishment effect test.
The 47 literature entries include repeated physical islands within WGSRPD3,
uncertain years, records only and historical events before NASA POWER 1981.
"""
import argparse,csv,json,re
from collections import Counter,defaultdict
from pathlib import Path

def read(path):
    with path.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f);return list(rd),rd.fieldnames

def native(v):
    return bool(v is not None and str(v).strip() not in ("","NA","NaN","nan","None","NULL","Unknown","unknown"))

def main():
    ap=argparse.ArgumentParser()
    for k in ("literature-seed","leptraits","resource-species","original-island-species","output-dir"):
        ap.add_argument("--"+k,type=Path,required=True)
    a=ap.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    literature,_=read(a.literature_seed)
    traits,fields=read(a.leptraits)
    resource,_=read(a.resource_species)
    orig,_=read(a.original_island_species)
    res={r["species"] for r in resource}
    original={r["species"] for r in orig}
    listed={r["species"] for r in literature}
    if len(res)!=239 or len(original)!=225 or len({(r["species"],r["physical_island"]) for r in literature})!=len(literature):
        raise RuntimeError("Source-panel mismatch")
    trait_columns={
      "wingspan":[c for c in fields if re.search(r"wing.?span|forewing|wing.?length|^(WS|FW)_",c,re.I)],
      "movement_or_migration":[c for c in fields if re.search("migra|dispers|mobil|movement",c,re.I)],
      "diapause_or_overwinter":[c for c in fields if re.search("diapause|overwinter|hibernat",c,re.I)],
      "oviposition":[c for c in fields if re.search("oviposit|egg.?lay",c,re.I)],
      "habitat_disturbance":[c for c in fields if re.search("disturbance|urban|habitat.?edge",c,re.I)],
      "voltinism":[c for c in fields if re.search("voltinism|voltine|brood",c,re.I)],
    }
    candidates=[c for c in fields if re.search(r"^species$|^scientific.?name$|^accepted.?name$|species.?name",c,re.I)]
    species_col=None
    rates=[]
    for col in candidates:
        intersection=len({r[col].strip() for r in traits if native(r.get(col))}&res)
        rates.append((intersection,col))
    if rates:
        _,species_col=max(rates)
    matching={}
    if species_col:
        for r in traits:
            key=r.get(species_col,"").strip()
            if key in res:
                matching.setdefault(key,[]).append(r)
    def coverage(subset, columns):
        return len([sp for sp in subset if any(
            native(r.get(c)) for r in matching.get(sp,[]) for c in columns)])
    trait_counts={
      group:{"columns":cols,
             "all_239_species_with_value":coverage(res,cols),
             "literature_seed_species_with_value":coverage(listed&res,cols),
             "nonmissing_in_each_source_column":{col:coverage(res,[col]) for col in cols}}
      for group,cols in trait_columns.items()
    }
    metrics_nonmissing={
      "wing_size_proxy":sum(native(r.get("wing_size_proxy")) for r in resource),
      "voltinism":sum(native(r.get("voltinism")) for r in resource),
      "host_family_count":sum(native(r.get("host_family_count")) for r in resource),
      "host_contemporary_units":sum(native(r.get("contemporary_resource_units")) for r in resource)
    }
    rows=[]
    for event in literature:
        sp=event["species"]
        start=int(event["first_evidence_earliest_year"]);end=int(event["first_evidence_latest_year"])
        if end<start or start<1750 or end>2026:
            raise RuntimeError("Invalid first-year event in literature seed: "+str(event))
        rows.append({
            **event,"in_GEB239":int(sp in res),"in_original_island_GBIF225":int(sp in original),
            "can_link_to_GEB_and_island_GBIF":int(sp in res and sp in original),
            "first_year_precise":int(start==end),
            "first_record_after_POWER_1981":int(start>=1981),
            "published_breeding_or_establishment_by_year":int(bool(event["established_by_year"])),
            "has_native_species_trait":int(sp in matching)
        })
    group=Counter((r["species"],r["wgsrpd3_code"]) for r in literature if r["wgsrpd3_code"])
    source_overlap=sum(r["can_link_to_GEB_and_island_GBIF"] for r in rows)
    summary={
       "schema":"chocho_three_gates_island_butterfly_case_traits_feasibility_v0.1",
       "status":"LITERATURE_CASE_LEADS_NOT_INDEPENDENT_ESTABLISHMENTS",
       "data":{
           "literature_events":len(rows),
           "unique_butterflies":len(listed),
           "distinct_physical_islands_or_localities":len({r["physical_island"] for r in literature}),
           "distinct_wgsrpd3_codes":len({r["wgsrpd3_code"] for r in literature if r["wgsrpd3_code"]}),
           "unique_butterfly_x_WGSRPD3":len(group),
           "within_WGSRPD3_redundant_events":sum(v-1 for v in group.values()),
           "source_events_missing_WGSRPD3_mapping":sum(not r["wgsrpd3_code"] for r in literature),
           "exact_year_first_record_cases":sum(r["first_year_precise"] for r in rows),
           "first_year_1981plus_cases":sum(r["first_record_after_POWER_1981"] for r in rows),
           "events_with_subsequent_establishment_or_breeding_upper_year":sum(r["published_breeding_or_establishment_by_year"] for r in rows),
           "event_rows_in_both_GEB239_and_GBIF225":source_overlap,
           "distinct_species_in_both_panels":len({r["species"] for r in rows if r["can_link_to_GEB_and_island_GBIF"]})
       },
       "GEB_resource_species_panel_size":len(res),
       "GEB_225_butterfly_port_observation_panel_size":len(original),
       "shared_butterfly_species_count":len(res&original),
       "leptraits":{"source_rows":len(traits),
                    "source_columns":fields,
                    "candidate_taxon_key_cols":rates,
                    "taxon_key_selected":species_col,
                    "frozen_species_matched":len(matching),
                    "by_trait":trait_counts},
       "GEB_precomputed_species_summary_completeness":metrics_nonmissing,
       "literature_unique_species_panel_link":{
          "in_GEB239_and_original225": sorted(listed&res&original),
          "missing_from_original225": sorted(listed-original),
       },
       "further_gates":[
         "Historical first record is not establishment and differs from known arrival.",
         "Independent observations on individual islands are collapsed within WGSRPD3 and cannot be treated as independent WGSRPD3 exposures.",
         "NASA POWER meteorological time series start 1981; pre-1981 butterfly first-record cases cannot receive direct climate-threshold chronology from POWER.",
         "For tested island outcomes, stable pre-arrival botanical and climatic denominators require period-specific effort and source history.",
         "Original 224/225-island butterfly occurrence parquet with 14.4+ million records is private/not mounted in this workflow; without source it cannot be used to infer all island first years.",
         "Natural dispersal, wind transport and human import are latent alternatives even for a butterfly with small wings.",
         "Wing size is a proxy (with limited prediction), not a deterministic flight-vs-ship classifier.",
         "Do not infer that occupation of newly available resource is a colonization event for native butterflies.",
       ]
    }
    for p,vals in [
      (a.output_dir/"case_coverage_by_panels.csv",rows)
    ]:
      with p.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(vals[0]))
        w.writeheader();w.writerows(vals)
    (a.output_dir/"case_trait_feasibility_audit.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k!="leptraits"},indent=2,ensure_ascii=False))
    print("TRAITS",json.dumps({k:v for k,v in summary["leptraits"].items() if k!="source_columns"},ensure_ascii=False))
if __name__=="__main__":
    main()

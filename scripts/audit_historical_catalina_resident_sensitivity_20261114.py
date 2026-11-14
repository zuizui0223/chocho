#!/usr/bin/env python3
"""Literature-dated historic resident counterexample to GBIF "first index" (Catalina).

This is an **externally motivated exclusion**: Powell (2004, rev.2012) lists
Pieris rapae on Santa Catalina by 1927, and Miller (1984/85) categorizes
Catalina Pieris rapae as resident. The GBIF search had none indexed during
2001-2010 and thus declared an index event in 2011-2017. That event is NOT
new island arrival. No other first-index cells are auto-reclassified.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path

from analyze_port_host_ocean_interaction_20261114 import source_join, model_sets
from audit_port_host_later_reindex_20261114 import evaluate

KEY=("1920496","951")
HISTORY={
    "species":"Pieris rapae",
    "GBIF_speciesKey":"1920496",
    "GSHHG_landmass_id":"951",
    "named_island":"Santa Catalina Island, California, USA",
    "historic_pre2001_occurrence_year":1927,
    "historic_residency_status":"R (resident) in Miller 1984/85 Santa Catalina column",
    "historical_sources":[
        {"author":"Powell","document_year":"2004; revised 2012",
         "url":"https://citeseerx.ist.psu.edu/document?doi=4696d7a674ce7b8b2360f7115929cce897c203b0&repid=rep1&type=pdf",
         "observation":"P. rapae in Catalina butterfly checklist: earliest 1927; later III.2011"},
        {"author":"Miller","year":"1984/85","journal":"Journal of Research on the Lepidoptera 23(4):282–296",
         "url":"https://repository.si.edu/server/api/core/bitstreams/3db06b68-e225-4ed1-9144-c994ad93779b/content",
         "observation":"Table 1 designates P. rapae resident (R) on Santa Catalina"},
    ],
    "note":"Old published records, not inferred solely from 2020/2022 larval photos"
}

def csv_rows(root,path):
    b=(root/path).read_bytes()
    return list(csv.DictReader(b.decode("utf-8").splitlines())),hashlib.sha256(b).hexdigest()

def reindex_case_rows(root):
    butterflies,bsha=csv_rows(root,"data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv")
    hosts,hsha=csv_rows(root,"results/route_free_port_hosts/species_island_hosts.csv")
    ports,psha=csv_rows(root,"results/route_free_port_hosts/port_coverage.csv")
    host={(x["key"],x["island"]):x for x in hosts}
    effort=defaultdict(int)
    for p in ports:effort[p["island"]]+=int(p["records"])
    cells=[]
    for b in butterflies:
        if int(b["eligible_2011_2017"])!=1 or int(b["first_indexed_in_2011_2017"])!=1:
            continue
        sp,isl=b["speciesKey"],b["landmass_id"]
        h=host.get((sp,isl))
        if not h or effort[isl]<10:continue
        cells.append({
            "island":isl,"species":sp,
            "introduced":int(int(h["introduced"])>0),
            "native":int(int(h["native"])>0),
            "early_records":int(b["future2011_2017_indexed_records"]),
            "reindexed":int(int(b["future2018_2022_indexed_records"])>0)
        })
    return cells,{"butterflies":bsha,"hosts":hsha,"ports":psha}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument("--out-dir",type=Path,default=Path("results/port_host_historic_resident_20261114"))
    a=ap.parse_args()
    rows,hashes,sea_beta,support=source_join(a.repo_root)
    old=model_sets(rows)
    matched=[r for r in rows["2011_2017"] if (r["sp"],r["isl"])==KEY]
    if len(matched)!=1 or matched[0]["event"]!=1 or matched[0]["host"]!=1:
        raise ValueError(f"Expected one early event with intro-host index, got {matched}")
    filtered={p:[r for r in sub if (r["sp"],r["isl"])!=KEY] for p,sub in rows.items()}
    changed=model_sets(filtered)
    before=old["2011_2017"]["plant_records_ge_10"]
    after=changed["2011_2017"]["plant_records_ge_10"]
    assert before["n"]-after["n"]==1 and before["events"]-after["events"]==1
    assert before["introduced_event_cells"]-after["introduced_event_cells"]==1
    assert old["2018_2022"]==changed["2018_2022"]

    cases,reindex_hashes=reindex_case_rows(a.repo_root)
    before_cases=[r for r in cases if (r["species"],r["island"])==KEY]
    if len(before_cases)!=1 or before_cases[0]["reindexed"]!=1 or before_cases[0]["introduced"]!=1:
        raise ValueError(f"Historic resident missing from observed re-index numerator: {before_cases}")
    cohort=evaluate(cases)
    corrected=evaluate([r for r in cases if (r["species"],r["island"])!=KEY])
    assert cohort["n_total"]==159 and corrected["n_total"]==158
    assert cohort["introduced_host_indexed"]["n"]==15 and corrected["introduced_host_indexed"]["n"]==14
    assert cohort["introduced_host_indexed"]["later_reindexed"]==11 and corrected["introduced_host_indexed"]["later_reindexed"]==10

    report={
        "study":"Historical independent butterfly checklist counterexample; source-index not novel colonization",
        "historical_case":HISTORY,
        "source_sha256":{**hashes,**{f"reindex_{k}":v for k,v in reindex_hashes.items()}},
        "source_status":support,
        "first_index_interaction":{
            "2011_2017_original":before,
            "2011_2017_excluding_historically_resident_951":after,
            "2018_2022_unchanged":True
        },
        "later_reindex":{
            "original_post_index_main_cohort":cohort,
            "excluding_historically_resident_951":corrected
        },
        "inference_limits":[
            "The historic resident counterexample applies to P. rapae on Catalina only",
            "Historical evidence is separate published literature, not original private GBIF parquet",
            "Reindexing is not novel settlement, breeding sustainability, or independent colonization",
            "Other apparent first indexes may likewise represent preexisting residents, unknown",
            "No plant host was directly observed feeding in the two original 2020/2022 GBIF textual records",
            "Original positive 239-butterfly host redistribution and 51-island community Circuitscape results are different endpoints and unchanged"
        ]
    }
    a.out_dir.mkdir(exist_ok=True,parents=True)
    path=a.out_dir/"historical_catalina_pre2001_counterexample.json"
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "output":str(path),
        "historic":HISTORY["historic_pre2001_occurrence_year"],
        "original_early_interaction":before["introduce_x_sea_access_effect_pp_per_SD"],
        "excl_historic_early_interaction":after["introduce_x_sea_access_effect_pp_per_SD"],
        "reindex_before":cohort["risk_difference_pp"],
        "reindex_excl_historic":corrected["risk_difference_pp"],
        "reindex_excl_intro_n":corrected["introduced_host_indexed"]["n"],
        "reindex_excl_intro_reindexed":corrected["introduced_host_indexed"]["later_reindexed"]
    },ensure_ascii=False))
if __name__=="__main__":
    main()

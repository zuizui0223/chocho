#!/usr/bin/env python3
"""Source-pinned crosswalk: historic island butterfly reports versus GBIF first indexes.

Independent field-source presence is NOT native origin, reproduction, or endemism.
A not-listed taxon is NOT a true absence; a first GBIF record is NOT arrival.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "gbif": ("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv",
             "b6a3c98eb76864052d45d520a55b7d830816126c"),
    "islands": ("results/island51_host_circuit_atlas_20261106/island_geography.csv",
                "f0c5b1139e54b4872d28682f0991f89f2368ab43"),
    "history": ("results/butterfly_GBIF_effort_origin_gate_20261117/PHUQUOC_PRE2011_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv",
                "c7dcc85a216ff5f76339f277bdeb6329d522a611"),
}
CHRISTMAS_SOURCE = "https://www.aph.gov.au/binaries/senate/committee/eca_ctte/estimates/add_0809/ewha/christmas_island_report_sept08.pdf"


def git_blob(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def read_pinned(key: str):
    rel, sha = PINS[key]
    p = ROOT / rel
    actual = git_blob(p)
    if actual != sha:
        raise ValueError(f"source drift for {key}: {actual}, expected {sha}")
    with p.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def pre2011_first_indices(rows):
    out = {}
    all_species = set()
    for row in rows:
        species = row["species"].strip()
        sk = int(row["speciesKey"])
        i = int(row["landmass_id"])
        a, b, c = [int(row[k]) for k in (
            "pre2001_2010_indexed_records",
            "future2011_2017_indexed_records",
            "future2018_2022_indexed_records")]
        if min(a, b, c) < 0:
            raise ValueError("negative source counts")
        all_species.add(species)
        if a == 0 and (b + c) > 0:
            k = (sk, i)
            if k in out:
                raise ValueError("duplicate species-island key")
            out[k] = {
                "species": species, "speciesKey": sk, "landmass_id": i,
                "earliest_GBIF_index_period": "2011-2017" if b else "2018-2022",
                "earliest_GBIF_index_year_lower_bound": 2011 if b else 2018,
                "GBIF_records_early": a,
                "GBIF_records_mid": b,
                "GBIF_records_late": c,
                "historical_field_presence": "unknown_not_in_source_scope",
                "last_possible_historical_presence_year": "",
                "source_url": "",
                "butterfly_origin": "unknown",
                "butterfly_establishment": "unknown",
                "butterfly_endemicity": "unassessed",
            }
    return out, all_species


def audit(rows, islands, history):
    if len(rows) != 11475:
        raise ValueError("original 225x51 source extent changed")
    risks, sp = pre2011_first_indices(rows)
    if len(risks) != 317:
        raise ValueError("later first-index count changed")
    ii = {int(x["landmass_id"]): x for x in islands}
    if len(ii) != 51:
        raise ValueError("island source changed")
    phu = ii[471]
    if (phu["example_port_name"], phu["example_port_country_code"]) != ("DUONG DONG", "VN"):
        raise ValueError("original Phu Quoc sampling identity changed")
    if not 103 < float(phu["representative_longitude"]) < 105:
        raise ValueError("Phu Quoc source coordinate mismatch")
    if not 9 < float(phu["representative_latitude"]) < 11:
        raise ValueError("Phu Quoc source coordinate mismatch")

    all_471 = [k for k in risks if k[1] == 471]
    if len(all_471) != 12:
        raise ValueError("changed 12-first-index candidate cohort on Phu Quoc")
    evidence_keys = set()
    detail = []
    for e in history:
        sk = int(e["speciesKey"])
        island = int(e["landmass_id"])
        k = (sk, island)
        if k in evidence_keys or k not in risks:
            raise ValueError("historical evidence taxon-island key not in first-index candidates or duplicated")
        evidence_keys.add(k)
        x = risks[k]
        if island != 471 or x["species"] != e["species"]:
            raise ValueError("historical source mapped to wrong island or wrong name")
        if x["earliest_GBIF_index_period"] != "2018-2022":
            raise ValueError("event must be post 2017")
        yr = int(e["earliest_independent_record_year_upper_bound"])
        publication = int(e["publication_year"])
        if yr not in (2006, 2007) or publication != yr or publication >= 2011:
            raise ValueError("historical source not actually before 2011")
        if e["island_name"] != "Phu Quoc":
            raise ValueError("wrong source locality")
        if not e["source_url"].startswith("https://wildlifeatrisk.org/wp-content/uploads/"):
            raise ValueError("unapproved historical source")
        if e["source_evidence_type"] not in ("field_checklist", "dated_specimen_or_field_image"):
            raise ValueError("unknown source type")
        if yr == 2006 and (e["printed_page"] != "4" or e["source_evidence_type"] != "field_checklist"):
            raise ValueError("2006 checklist not on expected page")
        if yr == 2007 and e["source_evidence_type"] != "dated_specimen_or_field_image":
            raise ValueError("2007 provenance or identity mismatch")
        if (e["butterfly_origin_status"],e["butterfly_establishment_status"],e["butterfly_endemicity_status"]) != ("unknown","unknown","unassessed"):
            raise ValueError("historical presence must not populate true origin/endemic/establishment")
        x["historical_field_presence"] = "source_verified_by_2007"
        x["last_possible_historical_presence_year"] = yr
        x["source_url"] = e["source_url"]
        detail.append({
            "species":e["species"],"speciesKey":sk,"landmass_id":island,
            "historical_year_upper_bound":yr, "source_printed_page":e["printed_page"],
            "independent_evidence_type":e["source_evidence_type"],
            "GBIF_first_index_window":x["earliest_GBIF_index_period"],
            "minimum_calendar_year_index_lag":2018-yr
        })
    if evidence_keys != set(all_471):
        raise ValueError("incomplete Phu Quoc 12 matched historical checklist cohort")
    if Counter(x["historical_year_upper_bound"] for x in detail) != {2006:9, 2007:3}:
        raise ValueError("2006 checklist + 2007 additional evidence count mismatch")

    mel = [x for x in rows if x["species"] == "Melanitis leda" and int(x["landmass_id"]) == 1210]
    if len(mel) != 1 or sum(int(mel[0][f]) for f in (
      "pre2001_2010_indexed_records","future2011_2017_indexed_records","future2018_2022_indexed_records")) != 0:
        raise ValueError("Christmas island independent vagrant/GBIF zero example no longer holds")
    if any(x in sp for x in ("Polyura andrewsii","Polyura andrewsi","Charaxes andrewsi","Charaxes andrewsii")):
        raise ValueError("the documented Christmas-Island endemic appears in 225 cohort; update bias statement")
    return {
        "schema":"chocho_independent_historic_butterfly_checklist_vs_GBIF_225x51_20261117_v1",
        "status":"SOURCE_PINNED_HISTORICAL_PRESENCE_DIAGNOSTIC_NOT_BIOLOGICAL_COLONISATION",
        "source_git_blob_shas":{x:sha for x,(path,sha) in PINS.items()},
        "full_cohort":{"butterflies":len(sp),"landmasses":len(ii),"species_landmass_cells":len(rows),
                       "GBIF_later_first_index_cells":len(risks),
                       "crosschecked_historical_independent_presence_cells":len(detail)},
        "Phu_Quoc":{
            "id":471,"original_WPI_port":"DUONG DONG","country_code":"VN",
            "GBIF_first_indices_2018_2022":12,
            "documented_2006_or_2007_butterflies":12,
            "independently_known_before_2011_fraction":1.0,
            "2006_field_checklist_matches":9,
            "additional_2007_field_specimen_matches":3,
            "all_12_early2001_2010_GBIF_counts":0,
            "minimum_evidence_to_GBIF_index_lag_years":11,
            "conclusion":"all 12 were already recorded on the same island by 2007; not new 2018-2022 arrivals",
            "source_timing_caveat":"source covers field sightings 2006-2007, first GBIF index interval begins 2018; not dated biological first occurrence",
            "evidence_entries":sorted(detail,key=lambda x:x["species"]),
        },
        "Christmas_Island_context":{
            "id":1210,"source":CHRISTMAS_SOURCE,"publication_year":2008,
            "Melanitis_leda":"government source identifies an earlier Christmas-Island vagrant (ssp. bankia) while selected 2001-2022 GBIF panel has 0 recorded specimens",
            "GBIF_counts_for_Melanitis_leda_2001_2022":0,
            "island_endemic_cited":"Polyura andrewsii (Christmas emperor; spelling/taxon concept in 2008 government report)",
            "island_endemic_present_in_225":False,
            "endemic_source_generalization":"At least one confirmed government-listed island endemic is not sampled by original HOSTS-linked 225 butterfly species, without estimating fraction missing"
        },
        "limits":["Evidence from Phu Quoc is selected, not an estimate of 317-cell population-wide historical misclassification rate",
                  "Both Phu Quoc brochures share Wildlife at Risk provenance; independent of GBIF period index, not independent biological sampling replicates",
                  "Earlier field presence is not proof of butterfly native range, breeding or endemicity",
                  "Source-missing butterfly is not assumed absent",
                  "Historical Christmas Island vagrant status should not be projected as current status without revalidation",
                  "195? no: the original denominators are 225 butterflies, 51 landmasses, 317 first-index records; preserve exact source populations"],
    }, list(risks.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,default=Path("_butterfly_historical_record_screen"))
    args = ap.parse_args()
    rows, islands, hist = (read_pinned(i) for i in ("gbif","islands","history"))
    report, annotated = audit(rows,islands,hist)
    # Sanity summary: unique all 225 species from the fixed source, 51 island groups
    args.out_dir.mkdir(parents=True,exist_ok=True)
    with (args.out_dir / "PHUQUOC_PRE2011_PRESENCE_VS_GBIF_20261117.json").open("w",encoding="utf8") as f:
        json.dump(report,f,indent=2,ensure_ascii=False)
        f.write("\n")
    with (args.out_dir / "ALL317_FIRSTINDEX_HISTORICAL_CROSSCHECK.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(annotated[0]))
        writer.writeheader()
        writer.writerows(sorted(annotated,key=lambda x:(x["landmass_id"],x["speciesKey"])))
    print(json.dumps({"cohort":report["full_cohort"],
          "Phu_Quoc":{k:v for k,v in report["Phu_Quoc"].items() if k!="evidence_entries"},
          "Christmas_Island_context":report["Christmas_Island_context"]},ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()

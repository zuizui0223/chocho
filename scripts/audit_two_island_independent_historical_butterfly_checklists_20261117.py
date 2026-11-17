#!/usr/bin/env python3
"""Audit independently published pre-2011 butterfly checklists against 317 later GBIF first indices.

The original 225 x 51 occurrence-index source remains unchanged.
Phu Quoc source has island-specific records, Sri Lanka's 2006 national checklist
has country-island-level scope, NOT within-island point coordinates. Neither is
a source of native origin, breeding, endemism, colonisation, or biological absence.
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
    "gbif": ("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv", "b6a3c98eb76864052d45d520a55b7d830816126c"),
    "islands": ("results/island51_host_circuit_atlas_20261106/island_geography.csv", "f0c5b1139e54b4872d28682f0991f89f2368ab43"),
    "phuquoc": ("results/butterfly_GBIF_effort_origin_gate_20261117/PHUQUOC_PRE2011_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv", "c7dcc85a216ff5f76339f277bdeb6329d522a611"),
    "srilanka": ("results/butterfly_GBIF_effort_origin_gate_20261117/SRILANKA_2006_IUCN_CHECKLIST_15_CANDIDATE_CROSSWALK_20261117.csv", "ba952428113ccb56936c2cf62e3f808ff328c7fc"),
}
SRI_LANKA_SOURCE = "https://portals.iucn.org/library/efiles/documents/2006-030.pdf"
SRI_LANKA_MIRROR = "https://lk.chm-cbd.net/sites/lk/files/2022-06/2006-030.pdf"


def git_blob(p):
    b = p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode("ascii")+b"\0"+b).hexdigest()


def pinned_table(key):
    relative, expected = PINS[key]
    p = ROOT / relative
    got = git_blob(p)
    if got != expected:
        raise ValueError(f"{key} source Git blob changed: expected {expected}, observed {got}")
    with p.open(newline="", encoding="utf-8-sig") as file:
        return list(csv.DictReader(file))


def jaccard(ids, sets):
    total = 0.
    n = 0
    for k,i in enumerate(ids):
        for j in ids[k+1:]:
            a,b = sets[i],sets[j]
            union = len(a | b)
            if not union:
                raise ValueError("Undefined original sampled pair Jaccard")
            total += len(a & b) / union
            n += 1
    return total / n, n


def audit(raw, islands, phuquoc, srilanka):
    if len(raw) != 11475:
        raise ValueError("Original fixed 225x51 species-island records changed")
    original_ids = {int(x["landmass_id"]): x for x in islands}
    if len(original_ids) != 51:
        raise ValueError("Original fixed island set changed")
    if original_ids[471]["example_port_name"] != "DUONG DONG" or original_ids[471]["example_port_country_code"] != "VN":
        raise ValueError("Phu Quoc original island geometry mismatch")
    if original_ids[31]["example_port_name"] != "KANKESANTURAI" or original_ids[31]["example_port_country_code"] != "LK":
        raise ValueError("Sri Lanka main island geometry mismatch")
    # The original WPI group 31 has three sampled ports on main-island Sri Lanka.
    if json.loads(original_ids[31]["source_port_names"]) != ["KANKESANTURAI", "COLOMBO", "GALLE HARBOR"]:
        raise ValueError("Sri Lanka WPI port roster changed")

    species = {x["speciesKey"]: x["species"] for x in raw}
    if len(species) != 225:
        raise ValueError("Butterfly source species changed")
    risk = {}
    early = {i: set() for i in original_ids}
    late = {i: set() for i in original_ids}
    mid_events = last_events = 0
    for x in raw:
        key = (int(x["speciesKey"]), int(x["landmass_id"]))
        name = x["species"]
        a,b,c = (int(x[z]) for z in ("pre2001_2010_indexed_records",
                                    "future2011_2017_indexed_records",
                                    "future2018_2022_indexed_records"))
        if min(a,b,c) < 0 or key[1] not in original_ids:
            raise ValueError("Invalid original cell")
        if a:
            early[key[1]].add(key[0])
        if a+b+c:
            late[key[1]].add(key[0])
        if not a and (b+c):
            if key in risk:
                raise ValueError("Duplicate first-index source identity")
            period = "2011-2017" if b else "2018-2022"
            mid_events += bool(b)
            last_events += not b
            risk[key] = dict(species=name, speciesKey=key[0], landmass_id=key[1],
                             earliest_GBIF_first_index_period=period,
                             GBIF_index_2001_2010=a, GBIF_index_2011_2017=b,
                             GBIF_index_2018_2022=c,
                             historical_status="not_audited",
                             historical_source_year="",
                             historical_geographic_resolution="",
                             original_historical_source_url="",
                             native_status="unknown", establishment_status="unknown",
                             island_endemicity="unassessed")
    if (len(risk),mid_events,last_events)!=(317,221,96):
        raise ValueError("Original later first-index cohort changed")

    checked = set()
    for e in phuquoc:
        k = (int(e["speciesKey"]), int(e["landmass_id"]))
        if k in checked or k not in risk or k[1] != 471 or risk[k]["species"] != e["species"]:
            raise ValueError("Phu Quoc source not matched to an original GBIF first index")
        checked.add(k)
        if int(e["earliest_independent_record_year_upper_bound"]) not in (2006,2007) or risk[k]["earliest_GBIF_first_index_period"]!="2018-2022":
            raise ValueError("Phu Quoc source timing mismatch")
        if (e["butterfly_origin_status"], e["butterfly_establishment_status"],
                e["butterfly_endemicity_status"]) != ("unknown","unknown","unassessed"):
            raise ValueError("Source historical presence is not origin/establishment/endemicity")
        risk[k].update(historical_status="source_confirmed_before_GBIF",
                       historical_source_year=e["earliest_independent_record_year_upper_bound"],
                       historical_geographic_resolution="named_island_field_report",
                       original_historical_source_url=e["source_url"])
    if len(checked)!=12:
        raise ValueError("Incomplete Phu Quoc prior source")

    checked_sri = set()
    sri_unmatched = set()
    for e in srilanka:
        k = (int(e["speciesKey"]),int(e["landmass_id"]))
        if k in checked_sri or k in sri_unmatched or k not in risk or k[1]!=31 or risk[k]["species"]!=e["species"]:
            raise ValueError("Sri Lanka 2006 source key/name/landmass mismatch")
        if int(e["checklist_publication_year"])!=2006 or e["island_name"]!="Sri Lanka" or e["source_url"]!=SRI_LANKA_SOURCE or e["government_mirror_url"]!=SRI_LANKA_MIRROR:
            raise ValueError("Sri Lanka publication or source mismatch")
        if (e["butterfly_origin_status"],e["butterfly_establishment_status"],e["butterfly_endemicity_status"]) != ("unknown","unknown","unassessed"):
            raise ValueError("A 2006 national checklist cannot prove origin, breeding or native endemicity")
        status=e["checklist_exact_species_match"]
        if status=="yes":
            page = int(e["printed_page"])
            if page not in (59,60,61,63,64) or not e["source_title"].startswith("Perera and Bambaradeniya (2006)"):
                raise ValueError("Historical appendix page or citation mismatch")
            checked_sri.add(k)
            risk[k].update(historical_status="source_confirmed_before_GBIF",
                           historical_source_year="2006",
                           historical_geographic_resolution="national_island_checklist_not_exact_point",
                           original_historical_source_url=SRI_LANKA_SOURCE)
        elif status=="not_found_in_2006_checklist":
            if e["species"]!="Ypthima huebneri" or e["printed_page"].strip():
                raise ValueError("Unexpected unresolved species or invalid source page")
            sri_unmatched.add(k)
            risk[k]["historical_status"]="searched_2006_source_not_found_unknown"
        else:
            raise ValueError("Unknown source status, fail closed")
    if (len(checked_sri),len(sri_unmatched))!=(14,1):
        raise ValueError("Expected 14 exact Sri Lanka 2006 source matches and one unresolved")
    if checked & checked_sri:
        raise ValueError("Historical same-cell duplication across island sources")
    if {k for k in risk if k[1]==31} != checked_sri | sri_unmatched:
        raise ValueError("Incomplete Sri Lanka 15 candidate cohort")
    if any(k not in risk for k in checked | checked_sri):
        raise ValueError("Historical source could not be matched")

    early_phu = {i:set(z) for i,z in early.items()}
    for sp,island in checked:
        early_phu[island].add(sp)
    early_2 = {i:set(z) for i,z in early_phu.items()}
    for sp,island in checked_sri:
        early_2[island].add(sp)
    index_ids = sorted(original_ids)
    initial, pi = jaccard(index_ids, early)
    phu, pp = jaccard(index_ids, early_phu)
    augmented, pa = jaccard(index_ids, early_2)
    last, pl = jaccard(index_ids, late)
    if len({pi,pp,pa,pl})!=1 or pi!=1275:
        raise ValueError("Fixed island-pair cohort changed")
    if (sum(map(len,early.values())),sum(map(len,early_phu.values())),
        sum(map(len,early_2.values())),sum(map(len,late.values()))) != (736,748,762,1053):
        raise ValueError("Cumulative inventory source-size mismatch")
    if abs(initial-0.1209632547327149)>1e-10 or abs(phu-0.1243590226479526)>1e-10:
        raise ValueError("Earlier independently audited results do not reproduce")
    if abs(augmented-0.12530899503574813)>1e-10:
        raise ValueError("Historical 2-island augmentation expected Jaccard mismatch")
    report = {
       "schema":"chocho_225x51_two_island_independent_2006_2007_historical_checklists_vs_later_GBIF_index_20261117_v1",
       "status":"SOURCE_PINNED_HISTORICAL_RECORD_LIST_DIAGNOSTIC_NOT_BIOLOGICAL_FIRST_COLONISATION",
       "input_git_blob_shas":{k:sha for k,(_,sha) in PINS.items()},
       "original":{"butterflies":225,"landmasses":51,"species_island_cells":11475,
                   "later_first_GBIF_index_cells":317,"first_index_2011_2017":221,
                   "first_index_2018_2022":96,"pairs":1275},
       "historical_checklist_matches":{
           "Phu_Quoc_landmass471":{"island_specific_earlier_matches":12,
                                   "published_2006":9,"published_2007":3,
                                   "all_late_GBIF_first_index_2018_2022":True,
                                   "source_geographic_resolution":"named_island_field_report"},
           "Sri_Lanka_landmass31":{"national_island_checklist_2006_exact_binomials":14,
                                    "later_GBIF_first_index_2011_2017":sum(risk[k]["earliest_GBIF_first_index_period"]=="2011-2017" for k in checked_sri),
                                    "later_GBIF_first_index_2018_2022":sum(risk[k]["earliest_GBIF_first_index_period"]=="2018-2022" for k in checked_sri),
                                    "unmatched_not_absent":["Ypthima huebneri"],
                                    "source_geographic_resolution":"national_island_checklist_not_exact_point",
                                    "source_year":2006,"source_url":SRI_LANKA_SOURCE},
           "total_source_confirmed_before_later_GBIF_index":26,
           "first_indices_not_audited_or_without_exact_2006_source_name":291,
           "never_infer_population_rate_from_selected_sources":True
       },
       "source_augmented_recorded_similarity":{
           "original_early_GBIF_indexed_cells":736,"Phu_Quoc_augmented_early_cells":748,
           "Phu_Quoc_and_Sri_Lanka_augmented_early_cells":762,
           "cumulative_2022_GBIF_indexed_cells":1053,
           "pairs":1275,
           "original_GBIF_only_early_Jaccard":initial,
           "Phu_Quoc_early_Jaccard":phu,
           "Phu_Quoc_Sri_Lanka_early_Jaccard":augmented,
           "GBIF_late_cumulative_Jaccard":last,
           "original_GBIF_list_delta":last-initial,
           "Phu_Quoc_augmentation_delta":last-phu,
           "two_sources_augmentation_delta":last-augmented,
           "two_sources_delta_change_from_original":initial-augmented,
           "interpretation":"Different observational sources augment early document-supported species lists; not same-protocol effort correction, not proof of native presence or biological range change"
       },
       "limitations":[
          "Sri Lanka national island checklist is island-wide but not a point-level record at WPI ports; the country and main island may differ at offshore islets",
          "Checklist absence for Ypthima huebneri is not true absence; no unverified synonym substitution",
          "Phu Quoc references are same-publisher continuing field work, not independent ecological replications",
          "Native butterfly distribution, endemism, breeding, and established nonnative populations cannot be determined by this crosswalk",
          "26 confirmed source matches reflect selected accessible published sources, not a 26/317 population-wide error rate",
          "Unreviewed remaining 291 first-index cases are neither true new island arrivals nor known old residents",
          "Historical report year bounds earliest evidence, not biological first observation or residency",
          "The original 225 HOSTS-associated butterflies are not a census of native endemic island butterflies"
       ]
    }
    return report,list(risk.values())


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=Path("_butterfly_two_island_historical_audit"))
    args=p.parse_args()
    tables=[pinned_table(k) for k in ("gbif","islands","phuquoc","srilanka")]
    report, rows=audit(*tables)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"TWO_ISLAND_HISTORICAL_GBIF_INDEX_SOURCE_AUDIT_20261117.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    with (args.out_dir/"ALL317_TWO_ISLAND_EARLIER_PUBLISHED_SOURCE_FLAGS_20261117.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows,key=lambda z:(z["landmass_id"],z["speciesKey"])))
    print(json.dumps({"historical":report["historical_checklist_matches"],
        "similarity":report["source_augmented_recorded_similarity"]},indent=2))
if __name__=="__main__":
    main()

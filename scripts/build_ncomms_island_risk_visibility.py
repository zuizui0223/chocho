#!/usr/bin/env python3
"""Build auditable source-record VISIBILITY for the NC butterfly x island risk set.

DO NOT use this output as an event/non-event training set. The 2011–2017
positive species-island event identities exist, while historic pre-2011
species-island identities and 2018–2022 residual positive identities are
NOT in this public source snapshot. Aggregate marginals alone do not
identify case/non-case rows. Never convert unknown into ecological absence.
All archived input files are checked against their original Git blob SHA1.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    "event_2011_2017": (
        "provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv",
        "efbc770469a720e2dea82478db2f6196ba03b353",
    ),
    "island_roster": (
        "provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv",
        "18cecde5bf3066f276380167f451f1fa5cc1caa1",
    ),
    "species_roster": (
        "provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv",
        "39fe84981cab2b0e47063399f0a5c0c9731e800d",
    ),
    "historical_general_flora": (
        "results/island51_pre2011_public_general_flora_20261107/source/GBIF_original51_early_general_flora_islands.csv",
        "048fa6d644f5579879620deb394ce74ace653e26",
    ),
    "original_port_coverage": (
        "provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_islands.csv",
        "143f9956ce4f871c7a97b4777ab69d416963c2bb",
    ),
    "original_17island_risk": (
        "provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_port_island_realization_risk_20261031.json",
        "25103acdd5cd0e8e28624ec92c340544f7a622ae",
    ),
}

MATRIX_FIELDS = [
    "species", "taxon_key", "landmass_id", "port_count",
    "historical_general_plant_speciesKeys_indexed_on_island",
    "original_pre2011_botanical_ports_complete",
    "butterfly_pre2011_index_at_THIS_species_island",
    "butterfly_first_index_2011_2017_at_THIS_species_island",
    "first_index_year_IF_2011_2017_case",
    "followup_years_count_IF_2011_2017_case",
    "modern_host_group_IF_2011_2017_case",
    "historical_2001_2010_HOSTS_record_at_THIS_species_island",
    "first_index_2018_2022_at_THIS_species_island",
    "full_index_2011_2022_at_THIS_species_island",
    "source_scope",
]
ISLAND_FIELDS = [
    "landmass_id", "ports", "early_butterfly_species_recorded_2001_2010",
    "new_butterfly_species_indexed_2011_2022",
    "identified_new_species_2011_2017",
    "unassigned_new_species_2018_2022",
    "historical_general_flora_speciesKeys", "historical_general_flora_ports_indexed",
    "historical_botanical_port_queries_complete",
]
SPECIES_FIELDS = [
    "species", "taxon_key", "early_indexed_islands_2001_2010",
    "new_indexed_islands_2011_2022", "identified_new_islands_2011_2017",
    "unassigned_new_islands_2018_2022",
]


def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()


def load_original(root: Path = ROOT) -> tuple[dict, dict]:
    inputs, provenance = {}, {}
    for key, (rel, expected) in FILES.items():
        raw = (root/rel).read_bytes()
        digest = git_blob_sha(raw)
        if digest != expected:
            raise ValueError("Source SHA mismatch: " + rel)
        inputs[key] = (json.loads(raw.decode("utf8")) if rel.endswith(".json")
                       else list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))))
        provenance[key] = {"path": rel, "git_blob_sha1": digest,
                           "source_sha256": hashlib.sha256(raw).hexdigest()}
    return inputs, provenance


def unique_rows(rows: list[dict], field: str, expected: int | None = None) -> dict[str, dict]:
    out = {}
    for row in rows:
        key = str(row[field])
        if key in out:
            raise ValueError("Duplicate "+field+": "+key)
        out[key] = row
    if expected is not None and len(out) != expected:
        raise ValueError("Unexpected "+field+" count "+str(len(out)))
    return out


def get_positive_int(value, name: str, *, allow_zero: bool = True) -> int:
    try:
        i = int(value)
    except (ValueError, TypeError) as error:
        raise ValueError(name+" not an integer") from error
    if str(value) != str(i) or i < 0 or (i == 0 and not allow_zero):
        raise ValueError(name+" invalid nonnegative integer")
    return i


def check_17island_archive(risk: dict) -> dict:
    if risk.get("schema") != "chocho_geb_port_island_species_at_risk_2011_2022_v1":
        raise ValueError("Original 17-island source schema changed")
    components = risk["geographic_components"]
    assert len(components) == 6
    geo = [island for block in components for island in block]
    if len(geo) != 17 or len(set(geo)) != 17:
        raise ValueError("Original 17-island source geography is inconsistent")
    periods = risk["periods"]
    result = {}
    for key, n, events, hosts, matched, nullmean, tail in [
        ("2011_2017_first_index_among_no_pre2011_butterfly_index",
         905, 113, 207, 34, 33.91983967935872, 1047),
        ("2018_2022_first_index_among_still_unindexed_after_2017",
         792, 33, 173, 20, 19.50367401469606, 770),
    ]:
        r = periods[key]["recording_risk"]
        null = periods[key]["geographic_margin_null"]
        a = r["host_early_recorded"]
        b = r["host_early_NOT_indexed"]
        if (int(r["cells"]) != n or int(r["indexed_first_events"]) != events
            or int(a["risk_cells"]) != hosts or int(a["new_indexed_events"]) != matched
            or int(b["risk_cells"])+hosts != n
            or int(b["new_indexed_events"])+matched != events):
            raise ValueError("Original 17-island period margins were modified")
        comp = r["by_geographic_component"]
        if len(comp) != 6 or {tuple(z["islands"]) for z in comp} != {tuple(z) for z in components}:
            raise ValueError("Original regional component identities were modified")
        for field, expected in [
            ("risk_cells", n),
            ("first_index_events", events),
            ("early_compatible_host_indexed_risk_cells", hosts),
            ("first_index_events_at_preindexed_compatible_host", matched),
        ]:
            if sum(int(z[field]) for z in comp) != expected:
                raise ValueError("Regional component margins not reconciled: "+field)
        if (int(null["total_draws"]) != 1497
            or int(null["observed_matched_events"]) != matched
            or int(null["reference_tail_draws"]) != tail
            or abs(float(null["null_mean_matched_events"])-nullmean)>1e-10):
            raise ValueError("Original conditional null mismatch")
        raw_host_rate = matched/hosts
        raw_without_rate = (events-matched)/(n-hosts)
        result[key] = {
            "risk_cells": n, "indexed_first_events": events,
            "historical_2001_2010_host_indexed_risk_cells": hosts,
            "historical_host_indexed_first_events": matched,
            "historical_host_NOT_indexed_risk_cells": n-hosts,
            "historical_host_NOT_indexed_first_events": events-matched,
            "raw_index_probability_given_recorded_prior_host": raw_host_rate,
            "raw_index_probability_given_NO_documented_prior_host": raw_without_rate,
            "raw_relative_record_index_rate": raw_host_rate/raw_without_rate,
            "record_difference_percent_points": 100*(raw_host_rate-raw_without_rate),
            "conditional_margin_null_mean_matched_events": nullmean,
            "observed_minus_conditional_null_expected_matched_events": matched-nullmean,
            "conditional_margin_null_draws_at_least_observed": tail,
            "conditional_margin_null_draws_total": 1497,
            "NOT_CAUSAL": True,
            "NONINDEXED_PLANT_IS_NOT_ABSENT": True,
            "selected_60_species_have_prior_2011_2017_index_event": True,
            "row_level_input_publicly_available": False,
        }
    if (int(risk["selected_name_join"]["resolved_larval_host_species"]) != 60
        or set(geo) - {str(x) for x in geo}):
        raise ValueError("Original 17-island source selection changed")
    return {"population_islands": len(geo), "population_butterfly_species": 60,
            "original_17_island_ids": geo, "periods": result,
            "source_observational_not_causal": True}


def build_visibility(source: dict, provenance: dict) -> tuple[list[dict], list[dict], list[dict], dict]:
    island_full = unique_rows(source["island_roster"],"landmass_id",103)
    islands = {k:z for k,z in island_full.items() if z["eligible"] == "True"}
    if len(islands) != 51:
        raise ValueError("The original fixed 51 eligible islands changed")
    species = unique_rows(source["species_roster"],"species",225)
    keys = unique_rows(source["species_roster"],"taxon_key",225)
    if any(z["species"] not in species or keys[z["taxon_key"]] != z for z in species.values()):
        raise ValueError("Species name and taxon key are inconsistent")
    early_general = unique_rows(source["historical_general_flora"],"island_id",51)
    early_ports = unique_rows(source["original_port_coverage"],"landmass_id",51)
    if set(early_general) != set(islands) or set(early_ports) != set(islands):
        raise ValueError("Current 51 island groups do not match both historical inventories")
    events = {}
    by_island = Counter()
    by_species = Counter()
    by_year = Counter()
    modern_groups = Counter()
    for row in source["event_2011_2017"]:
        name, isl, year = row["species"], row["landmass_id"], int(row["first_record_year"])
        if not 2011 <= year <= 2017:
            raise ValueError("Non-2011–2017 first record in the archived event file")
        if not (name in species and isl in islands):
            continue
        key=(name,isl)
        if key in events:
            raise ValueError("Duplicated archived positive butterfly×island event")
        if row["speciesKey"] != species[name]["taxon_key"]:
            raise ValueError("Historical event speciesKey does not match original 225 species roster")
        if row["host_group"] == "NOT_EVALUABLE":
            raise ValueError("Original 51 island groups must have evaluable pooled host context")
        n_follow=get_positive_int(row["followup_record_years"],"followup record years")
        if n_follow > 5 or n_follow != (len(row["followup_year_list"].split(";")) if row["followup_year_list"] else 0):
            raise ValueError("Archived later re-record years inconsistent")
        events[key]=row
        by_island[isl]+=1
        by_species[name]+=1
        by_year[year]+=1
        modern_groups[row["host_group"]]+=1
    if len(events) != 302 or sum(by_island.values()) != 302:
        raise ValueError("Unexpected originally published 2011–2017 first-index count")

    species_rows = []
    for name, row in sorted(species.items(), key=lambda x:x[0]):
        old = get_positive_int(row["earlier_records_on_eligible_islands"],"earlier species margin")
        later = get_positive_int(row["new_records_on_eligible_islands"],"later species margin")
        first = by_species[name]
        if old + later > 51 or first > later:
            raise ValueError("Inconsistent species early/later marginal records: "+name)
        species_rows.append({
            "species":name, "taxon_key":row["taxon_key"],
            "early_indexed_islands_2001_2010":old,
            "new_indexed_islands_2011_2022":later,
            "identified_new_islands_2011_2017":first,
            "unassigned_new_islands_2018_2022":later-first,
        })
    island_rows=[]
    for isl, row in sorted(islands.items(),key=lambda x:int(x[0])):
        general=early_general[isl]
        ports=early_ports[isl]
        nports=get_positive_int(row["ports"],"original ports",allow_zero=False)
        if (nports != int(general["original_port_count"])
            or nports != int(ports["planned_ports"])
            or nports != int(ports["complete_ports"])
            or ports["all_ports_complete"] != "True"):
            raise ValueError("Original island port counts or query completeness changed")
        old=get_positive_int(row["earlier_butterfly_species"],"island earlier margin")
        later=get_positive_int(row["new_butterfly_species"],"island later margin")
        nfirst=by_island[isl]
        if old + later > 225 or nfirst > later:
            raise ValueError("Inconsistent earlier/later island marginals for "+isl)
        plant=get_positive_int(general["unique_early_indexed_general_plant_speciesKeys"],"indexed general plants")
        island_rows.append({
            "landmass_id":isl, "ports":nports,
            "early_butterfly_species_recorded_2001_2010":old,
            "new_butterfly_species_indexed_2011_2022":later,
            "identified_new_species_2011_2017":nfirst,
            "unassigned_new_species_2018_2022":later-nfirst,
            "historical_general_flora_speciesKeys":plant,
            "historical_general_flora_ports_indexed":int(general["ports_with_plant_record"]),
            "historical_botanical_port_queries_complete":"TRUE",
        })
    fields = {}
    source_unknown = "NOT_RECOVERED_IN_PUBLIC_COMMITTED_ROW_LEVEL_SOURCE"
    modern_unknown = "NOT_AN_INDEXED_2011_2017_EVENT_ROW"
    for row in island_rows:
        fields[row["landmass_id"]]=row
    matrix=[]
    for s in species_rows:
        name=s["species"]
        for island in island_rows:
            lid=island["landmass_id"]
            event=events.get((name,lid))
            matrix.append({
                "species":name, "taxon_key":s["taxon_key"],"landmass_id":lid,
                "port_count":island["ports"],
                "historical_general_plant_speciesKeys_indexed_on_island":island["historical_general_flora_speciesKeys"],
                "original_pre2011_botanical_ports_complete":"TRUE",
                "butterfly_pre2011_index_at_THIS_species_island":(
                    "KNOWN_NOT_PRE2011" if event else source_unknown),
                "butterfly_first_index_2011_2017_at_THIS_species_island":(
                    "CONFIRMED_FIRST_INDEX" if event else source_unknown),
                "first_index_year_IF_2011_2017_case":event["first_record_year"] if event else "",
                "followup_years_count_IF_2011_2017_case":event["followup_record_years"] if event else "",
                "modern_host_group_IF_2011_2017_case":event["host_group"] if event else modern_unknown,
                "historical_2001_2010_HOSTS_record_at_THIS_species_island":source_unknown,
                "first_index_2018_2022_at_THIS_species_island":(
                    "KNOWN_NOT_FIRST_2018_2022" if event else source_unknown),
                "full_index_2011_2022_at_THIS_species_island":(
                    "CONFIRMED_FIRST_INDEX_2011_2017" if event else source_unknown),
                "source_scope":"ORIGINAL_225_SPECIES_51_ISLANDS_INDEX_VISIBILITY_NOT_A_RISK_PREDICTION_TABLE",
            })
    if len(matrix) != 225*51:
        raise ValueError("Full source-pair grid is not rectangular")
    prior_island=sum(int(z["early_butterfly_species_recorded_2001_2010"]) for z in island_rows)
    prior_species=sum(int(z["early_indexed_islands_2001_2010"]) for z in species_rows)
    new_island=sum(int(z["new_butterfly_species_indexed_2011_2022"]) for z in island_rows)
    new_species=sum(int(z["new_indexed_islands_2011_2022"]) for z in species_rows)
    for lhs,rhs,expected,name in [(prior_island,prior_species,524,"prior"),
                                  (new_island,new_species,430,"later")]:
        if lhs!=rhs or lhs!=expected:
            raise ValueError("Original island/species "+name+" margins disagree")
    late_island=sum(z["unassigned_new_species_2018_2022"] for z in island_rows)
    late_species=sum(z["unassigned_new_islands_2018_2022"] for z in species_rows)
    if late_island != 128 or late_species != 128:
        raise ValueError("Unknown later-event identities do not reconcile")
    allpossible=225*51
    missing=allpossible-prior_island-new_island
    if missing!=10521:
        raise ValueError("Unindexed-through-2022 aggregate mismatches")
    original_17=check_17island_archive(source["original_17island_risk"])
    if not set(original_17["original_17_island_ids"]).issubset(set(islands)):
        raise ValueError("Previously published 17 islands are not a subset of original 51")

    audit = {
        "schema":"chocho_NC_island225x51_temporal_risk_visibility_v1",
        "status":"SOURCE_VERIFIED_VISIBILITY_PANEL_NOT_READY_FOR_EVENT_PREDICTION",
        "original_source_provenance":provenance,
        "original_candidates":225,"original_island_groups":51,
        "all_potential_species_by_island_cells":allpossible,
        "published_pre2011_butterfly_indexed_species_island_cells":prior_island,
        "published_first2011_2022_indexed_species_island_cells":new_island,
        "archived_row_identified_first2011_2017_species_island_cells":len(events),
        "first2011_2017_species":len(by_species),
        "first2011_2017_islands":len(by_island),
        "first2011_2017_by_year":dict(sorted(by_year.items())),
        "first2011_2017_event_modern_pooled_HOST_group":dict(sorted(modern_groups.items())),
        "unassigned_2018_2022_positive_identities_but_known_marginal_total":late_species,
        "all_cells_not_indexed_in_original_butterfly_2001_2010_by_count":allpossible-prior_island,
        "all_cells_without_butterfly_index_in_original_2001_2022_by_count":missing,
        "these_nonindexed_cells_are_NOT_established_absences":True,
        "positive_first2011_2017_event_rows_can_be_joined_by_speciesKey_and_island":True,
        "original_complete_species_by_island_pre2011_index_identity_publicly_available":False,
        "original_full_positive_2018_2022_species_island_event_identities_publicly_available":False,
        "original_species_by_island_prior_HOSTS_identity_publicly_available":False,
        "all_51_island_plant_query_port_completeness_status":"51_OF_51",
        "early_general_plant_islands_with_any_indexed_speciesKey":sum(int(x["historical_general_flora_speciesKeys"])>0 for x in island_rows),
        "modern_HOST_group_is_2001_2024_pooled_and_NOT_time_prior":True,
        "original_17island_prehost_risk_subset_aggregate":original_17,
        "cannot_fit_valid_host_times_access_risk_model":True,
        "prohibited_recode":"Never make 11,173 unassigned matrix rows into true negative/at-risk cases; the source only gives marginal pre-2011 and later counts.",
        "next_source_recovery":[
            "Original 83,573-row 2001–2024 species×island×year Parquet with frozen SHA256 8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f",
            "Original 2001–2010 historical HOSTS×port species identity, or an independently dated host accepted-name key panel",
            "Observational opportunities / detection effort for both recorded and unrecorded species×island windows",
        ],
        "no_evidence_of_causal_shipping_or_establishment":True,
    }
    return matrix,island_rows,species_rows,audit


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w",encoding="utf8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    source,provenance=load_original()
    grid,islands,species,report=build_visibility(source,provenance)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    files={
        "visibility":args.output_dir/"species_island_225x51_visibility_NOT_EVENT_RISK.csv",
        "island":args.output_dir/"original51_source_index_margins.csv",
        "species":args.output_dir/"original225_source_index_margins.csv",
        "audit":args.output_dir/"risk_visibility_audit.json",
    }
    write_csv(files["visibility"],grid,MATRIX_FIELDS)
    write_csv(files["island"],islands,ISLAND_FIELDS)
    write_csv(files["species"],species,SPECIES_FIELDS)
    # The audit JSON cannot contain its own SHA-256: doing so would make
    # the receipt self-referential and the published hash immediately stale.
    report["output_csv_sha256"]={
        name:hashlib.sha256(path.read_bytes()).hexdigest()
        for name,path in files.items() if name!="audit"
    }
    files["audit"].write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
        "status":report["status"],
        "visible_2011_2017_events":report["archived_row_identified_first2011_2017_species_island_cells"],
        "unknown_2018_2022_event_IDs":report["unassigned_2018_2022_positive_identities_but_known_marginal_total"],
        "known_prior_species_island_marginal":report["published_pre2011_butterfly_indexed_species_island_cells"],
        "17island_2011_2017_case_count":report["original_17island_prehost_risk_subset_aggregate"]["periods"][
            "2011_2017_first_index_among_no_pre2011_butterfly_index"]["indexed_first_events"],
        "17island_null_nonexceptional":True,
        "out_dir":str(args.output_dir),
    },indent=2))


if __name__=="__main__":
    main()

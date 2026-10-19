#!/usr/bin/env python3
"""Lossless merge of the 135 originally complete GBIF taxa and 96 recovery attempts.

Never convert non-fetched taxa into ecological absence. Raw GBIF coordinates
are not retained. Output only public sampled species×WGSRPD3×period×basis
record counts and complete/partial source receipts. The original source
ordinal sampling budget remains exactly two 300-index windows per period.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

FIELDS=("species","window","wgsrpd3","basis","sampled_record_count")
EXPECTED_ORIGINAL_STATUS={"COMPLETE":135,"PARTIAL_FETCH_NOT_ELIGIBLE":96,"TAXON_UNRESOLVED":8}


def source_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as fh:
        reader=csv.DictReader(fh)
        if tuple(reader.fieldnames or ())!=FIELDS:
            raise ValueError(f"Original/recovery mapped occurrence schema drift: {path}")
        result={}
        for row in reader:
            key=tuple(row[k] for k in FIELDS[:-1])
            n=int(row["sampled_record_count"])
            if key in result or n<=0:
                raise ValueError("Nonpositive or duplicated GBIF source count row")
            result[key]=n
        return result


def assemble(original_manifest,original_counts,recovery_manifest,recovered_counts):
    raw=original_manifest
    rev=recovery_manifest
    if raw.get("schema")!="chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01":
        raise ValueError("Original schema is not the versioned source")
    if raw.get("species_requested")!=239 or raw.get("species_statuses")!=EXPECTED_ORIGINAL_STATUS:
        raise ValueError("Original GBIF 135/96/8 count audit changed")
    if rev.get("schema")!="chocho_GBIF_96_partial_species_transportreread_20261019_v01":
        raise ValueError("Recovery source schema mismatch")
    if rev.get("species_attempted")!=96 or rev.get("original_ordinal_size")!=300 or rev.get("transport_only_subrequest_size")!=75:
        raise ValueError("Recovery has a different sampling scope")
    old_species={r["species"]:r for r in raw["species_receipts"]}
    recovered={r["species"]:r for r in rev["species_receipts"]}
    assert len(old_species)==239 and len(recovered)==96
    partial={s for s,row in old_species.items() if row["status"]=="PARTIAL_FETCH_NOT_ELIGIBLE"}
    if set(recovered)!=partial:
        raise ValueError("Recovery not exactly original incomplete 96")
    source_keys={}
    accepted=[]
    for name,row in recovered.items():
        prior=old_species[name]
        key=int(prior["taxon"]["key"])
        if key!=int(row["accepted_species_key"]):
            raise ValueError("GBIF taxon key changed in recovery, nonreproducible taxonomy")
        if row["status"]=="COMPLETE":
            if set(row["windows"])!={"2010_2017","2018_2025"} or not all(x.get("complete") for x in row["windows"].values()):
                raise ValueError("Incomplete window called complete")
            accepted.append(name)
    old_rows=source_csv(original_counts)
    added_rows=source_csv(recovered_counts)
    accepted=set(accepted)
    if any(name not in accepted for name,_,_,_ in added_rows):
        raise ValueError("Recovery partial species records illegally counted as presence")
    if any(name not in {s for s,row in old_species.items() if row["status"]=="COMPLETE"} for name,_,_,_ in old_rows):
        raise ValueError("Source original contains uncensored partial species records")
    if set(old_rows)&set(added_rows):
        raise ValueError("Original and recovery duplicated aggregated occurrence")
    combined={**old_rows,**added_rows}
    new_receipts=[]
    for name,prior in sorted(old_species.items()):
        now=prior.copy()
        if name in accepted:
            r=recovered[name]
            now["status"]="COMPLETE"
            now["windows"]=r["windows"]
            now["errors"]=[]
            now["originally_incomplete_source_recovered"]=True
            now["transport_microchunk_ordinal_budget_preserved"]=True
            now["raw_unmapped_after_recovery"]=r.get("source_validity_counters",{})
            now["original_curl_error_count"]=len(prior.get("errors",[]))
        elif name in recovered:
            now["microchunk_recovery_status"]=recovered[name]["status"]
            now["microchunk_recovery_errors"]=recovered[name]["errors"]
        new_receipts.append(now)
    status=Counter(r["status"] for r in new_receipts)
    n_complete=status["COMPLETE"]
    if n_complete!=135+len(accepted) or status["TAXON_UNRESOLVED"]!=8:
        raise ValueError("Taxonomic completeness arithmetic mismatch")
    expected_n=sum(combined.values())
    result={
        **raw,
        "status":"PARTIAL_GBIF_8_UNRESOLVED_OR_TRANSPORT_NOT_BIOLOGICAL_ABSENCE" if n_complete<239 else "SOURCE_COMPLETED",
        "merged_from_original_GitHub_run_id":39233336804,
        "recovered_from_original_partial96_GitHub_run_id":"see linked recovered artifact",
        "merged_utc":datetime.now(timezone.utc).isoformat(),
        "species_taxon_resolved_and_complete":n_complete,
        "species_statuses":dict(sorted(status.items())),
        "species_receipts":new_receipts,
        "outputs":{"aggregated_region_basis_rows":len(combined),"source_records_mapped":expected_n},
        "sampling_is_systematic_and_not_random_or_exhaustive":True,
        "merge_metadata":{
            "original_complete":135,
            "recovered_new_complete":len(accepted),
            "still_transport_incomplete":status["PARTIAL_FETCH_NOT_ELIGIBLE"],
            "taxon_unresolved":8,
            "no_new_species_outside_initial_239":True,
            "no_ordinal_sample_budget_change":True,
            "recovery_count_drift_is_exposed_in_receipts":True,
            "no_absence_inferred_from_missing_recovered_taxa":True,
        }
    }
    return combined,result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("original-manifest","original-counts","recovery-manifest","recovery-counts",
                "output-manifest","output-counts"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    for f in (a.original_manifest,a.recovery_manifest):
        if not f.is_file():raise FileNotFoundError(f)
    source=json.loads(a.original_manifest.read_text(encoding="utf-8"))
    recovery=json.loads(a.recovery_manifest.read_text(encoding="utf-8"))
    out,metadata=assemble(source,a.original_counts,recovery,a.recovery_counts)
    a.output_counts.parent.mkdir(parents=True,exist_ok=True)
    with a.output_counts.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.writer(fh)
        writer.writerow(FIELDS)
        for key,n in sorted(out.items()):writer.writerow([*key,n])
    a.output_manifest.parent.mkdir(parents=True,exist_ok=True)
    a.output_manifest.write_text(json.dumps(metadata,indent=2,sort_keys=True)+"\n",
                                 encoding="utf-8")
    print(json.dumps({
        "newly_completed_species":metadata["merge_metadata"]["recovered_new_complete"],
        "all_complete_species":metadata["species_taxon_resolved_and_complete"],
        "remaining_missing":metadata["merge_metadata"]["still_transport_incomplete"],
        "unresolved_taxa":metadata["merge_metadata"]["taxon_unresolved"],
        "unique_aggregates":metadata["outputs"]["aggregated_region_basis_rows"],
        "all_mapped_completed_records":metadata["outputs"]["source_records_mapped"],
        "status":metadata["status"],
    },sort_keys=True))

#!/usr/bin/env python3
"""Auditable location and evidence screen for the pinned NHM HOSTS source.

Source audit only. It does not treat any named location as a butterfly
occurrence site, larval wild feeding, or a WGSRPD3 mapping until classified.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re
from collections import Counter, defaultdict
from pathlib import Path

EXPECTED_COLUMNS = ("_id","Hostplant Genus","Hostplant Species",
                    "Insect Genus","Insect Species","Location","Lab Rearing")
EXPECTED_HOSTS_COMMIT = "808e0b869f9ec1adf8efff87cf6a395adda103e0"


def norm(x: str) -> str:
    return re.sub(r"\s+", " ", x.strip())


def main() -> None:
    a=argparse.ArgumentParser()
    a.add_argument("--hosts",required=True,type=Path)
    a.add_argument("--crop",required=True,type=Path)
    a.add_argument("--ports",required=True,type=Path)
    a.add_argument("--out-dir",required=True,type=Path)
    x=a.parse_args()
    locations=Counter()
    raw_locations=Counter()
    rearing=Counter()
    damage=Counter()
    samples=defaultdict(list)
    locale_rearing=defaultdict(Counter)
    host_status=Counter()
    rows=0
    ids=set()
    with x.hosts.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f)
        fieldnames=reader.fieldnames or []
        missing=set(EXPECTED_COLUMNS)-set(fieldnames)
        if missing:raise RuntimeError(f"Original HOSTS headers missing: {sorted(missing)}")
        for r in reader:
            rows+=1
            ident=r.get("_id","")
            if ident and ident in ids:
                host_status["duplicate_id"]+=1
            if ident:ids.add(ident)
            loc=norm(r.get("Location") or "")
            lab=norm(r.get("Lab Rearing") or "")
            dmg=norm(r.get("Damage") or "")
            rearing[lab]+=1
            damage[dmg]+=1
            host=norm(r.get("Hostplant Genus") or "")+" "+norm(r.get("Hostplant Species") or "")
            insect=norm(r.get("Insect Genus") or "")+" "+norm(r.get("Insect Species") or "")
            if loc:
                locations[loc.casefold()]+=1
                raw_locations[loc]+=1
                locale_rearing[loc][lab]+=1
                if len(samples[loc])<2 and host.strip() and insect.strip():
                    samples[loc].append({"host":host.strip(),"insect":insect.strip(),"lab":lab,"damage":dmg})
            else:host_status["no_location"]+=1
            if not all(r.get(k) for k in ("Hostplant Genus","Hostplant Species",
                                          "Insect Genus","Insect Species")):
                host_status["incomplete_species"]+=1
    crop=json.loads(x.crop.read_text(encoding="utf8"))
    if crop.get("schema")!="chocho_fao_crop_taxa_v0.1":raise ValueError("FAO taxonomy schema changed")
    with x.ports.open(newline="",encoding="utf8") as f:
        ports=list(csv.DictReader(f))
    if len(ports)<300 or len({r["area_code_l3"] for r in ports})!=len(ports):
        raise ValueError("WGSRPD3 WPI port source incomplete")
    x.out_dir.mkdir(parents=True,exist_ok=True)
    top=sorted(raw_locations.items(),key=lambda z:(-z[1],z[0]))[:600]
    with (x.out_dir/"top600_HOSTS_raw_localities.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.writer(f)
        w.writerow(["HOSTS_Location_exact","rows","lab_rearing_value_count","example_insect_host_record_1","example_insect_host_record_2"])
        for loc,n in top:
            w.writerow([loc,n,json.dumps(locale_rearing[loc],ensure_ascii=False),
                        json.dumps(samples[loc][0],ensure_ascii=False) if samples[loc] else "",
                        json.dumps(samples[loc][1],ensure_ascii=False) if len(samples[loc])>1 else ""])
    audit={
        "status":"PINNED_HOSTS_SOURCE_LOCATION_INVENTORY_NOT_LOCAL_USE_VERIFICATION",
        "hosts_commit":EXPECTED_HOSTS_COMMIT,
        "hosts_csv_sha256":hashlib.sha256(x.hosts.read_bytes()).hexdigest(),
        "hosts_source_rows":rows,"source_column_names":fieldnames,
        "record_id_count":len(ids),
        "distinct_exact_nonempty_location_strings":len(raw_locations),
        "distinct_normalized_location_strings":len(locations),
        "rows_with_location":sum(raw_locations.values()),
        "rows_without_location":host_status["no_location"],
        "rows_incomplete_species":host_status["incomplete_species"],
        "duplicate_ids":host_status["duplicate_id"],
        "lab_rearing_values":dict(rearing),
        "damage_values":dict(damage),
        "FAO_exact_binomials":len(crop["exact_binomials"]),
        "FAO_genus_wildcards":len(crop["genus_wildcards"]),
        "FAO_html_sha256":crop["source_html_sha256"],
        "WGSRPD3_port_regions":len(ports),
        "top60_original_localities":[{"location":a,"n":n} for a,n in top[:60]],
        "interpretation_gates":[
            "HOSTS Location is a bibliographic free-text field, not a verified geographic coordinate.",
            "HOSTS Lab Rearing must be audited; records may be laboratory/experimental or cited from secondary literature.",
            "Converting an insect-host literature row to a field-confirmed native/introduced plant feeding event requires source verification.",
            "WCVP introduced applies to WGSRPD3 distributions, not individual field sites.",
            "FAO crop exact binomials, plus genus wildcard stress test, do not establish crop provenance of every HOSTS row.",
            "No location-name inference should be imputed for missing or ambiguous raw labels.",
        ],
    }
    (x.out_dir/"HOSTS_location_source_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("HOSTS_ALL_RAW_LOCALITY_COUNTS",json.dumps(
        [{"location": a, "n": n} for a,n in top],ensure_ascii=False),flush=True)
    print("HOSTS_ORIGINAL_SOURCE_LOCATION_AUDIT",json.dumps({
        k:audit[k] for k in ("hosts_source_rows","rows_with_location","distinct_exact_nonempty_location_strings",
                             "lab_rearing_values","FAO_exact_binomials","FAO_genus_wildcards","WGSRPD3_port_regions",
                             "top60_original_localities")},ensure_ascii=False),flush=True)


if __name__=="__main__":main()

#!/usr/bin/env python3
"""Source-lock NHM HOSTS host plant SPECIES for the independently classified Azores butterflies.

No deduction from taxa selection. Two Azores endemics are included here
because the independent 2023 University of Azores source confirms their
origin; source scientific name H. azorinus is linked to INBO H. azorina by
the 2018 European butterfly checklist original combination Satyrus azorinus.
"""
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path

ORIG=Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117/independent_42_species_vs_frozen239.csv")
HOSTS_SHA="808e0b869f9ec1adf8efff87cf6a395adda103e0"
NAMES={"Colias crocea","Danaus plexippus","Hipparchia azorina",
       "Hipparchia miguelensis","Lampides boeticus",
       "Pieris brassicae","Vanessa atalanta","Vanessa cardui"}
AZORES_END={"Hipparchia azorina","Hipparchia miguelensis"}
HOST_FIELDS=("Insect Genus","Insect Species","Insect Family","Hostplant Genus","Hostplant Species")

def src_rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))

def validate(hosts_file,source_json,output):
    independent=src_rows(ORIG)
    if len(independent)!=42:raise ValueError("Source independent 42 changed")
    azores={x["scientific_name"]:x for x in independent if "Azores" in x["archipelagos"]}
    if set(azores)!=NAMES:raise ValueError("Original independent eight Azores butterflies changed")
    official=json.loads(source_json.read_text(encoding="utf-8"))
    assert official["counts"]["taxa"]==4614 and official["counts"]["distribution"]==8630
    source_status={}
    for r in official["exact_taxon_results"]:
        name=r["independent_INBO_species"]
        means=r["official_establishment_means"] or r.get("alias_official_origin_status",[])
        if len(means)!=1 or means[0] not in ("Native","Endemic"):
            raise ValueError("Unresolved original 2023 source status: "+name)
        source_status[name]=means[0]
    if set(source_status)!=NAMES or set(x for x in source_status if source_status[x]=="Endemic")!=AZORES_END:
        raise ValueError("Official Azores endemic identity changed")
    index=defaultdict(list)
    total=0
    with hosts_file.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        cols=set(reader.fieldnames or ())
        if not set(HOST_FIELDS)<=cols:
            raise ValueError("NHM HOSTS original raw fields changed")
        reference_fields=[x for x in (reader.fieldnames or ()) if any(
            token in x.lower() for token in ("reference","source","citation","locality","country","date","note","page","method"))]
        for x in reader:
            total+=1
            if (x.get("Insect Family") or "").lower() not in {
                "nymphalidae","pieridae","lycaenidae"}:continue
            name=((x.get("Insect Genus") or "").strip()+" "+(x.get("Insect Species") or "").strip()).strip()
            if name in NAMES:
                index[name].append({
                    "name":name,
                    "family":x.get("Insect Family",""),
                    "host_genus":x.get("Hostplant Genus",""),
                    "host_epithet":x.get("Hostplant Species",""),
                    "host_binomial":" ".join((x.get("Hostplant Genus") or "").strip().split()[:1]+(x.get("Hostplant Species") or "").strip().split()[:1]).strip(),
                    "original_reference_fields":{key:x.get(key) for key in reference_fields},
                })
    if total!=140485:raise ValueError("Pinned HOSTS mirror row count drift")
    endemic_records=[]
    for name in sorted(NAMES):
        r=azores[name]
        records=index[name]
        named=[e for e in records if e["host_genus"] and e["host_epithet"] and
          e["host_epithet"].lower() not in ("sp","spp","unknown","indet")]
        if r["hosts_named_species"]!="1" or not named:
            raise ValueError("Frozen independent eight exact named-species HOSTS status drift: "+name)
        for e in records:
            e["independent_azores_butterfly_origin"]=source_status[name]
            e["independent_azores_species"]=name
            if name in AZORES_END:
                endemic_records.append(e)
    if len(endemic_records)<2:raise ValueError("At least one official Azores endemic original host missing")
    output.mkdir(parents=True,exist_ok=True)
    fieldnames=["independent_azores_species","independent_azores_butterfly_origin",
        "host_binomial","host_genus","host_epithet","family","name",
        "original_reference_fields"]
    with (output/"Azores_8_host_species_NHM_original.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames);writer.writeheader()
        for name in sorted(NAMES):
            for obj in index[name]:
                rec={k:obj.get(k,"") for k in fieldnames}
                rec["original_reference_fields"]=json.dumps(rec["original_reference_fields"],ensure_ascii=False)
                writer.writerow(rec)
    summary={
      "schema":"chocho_azores2023_official_NHM_HOSTS_species_source_v1",
      "source_original_HOSTS_commit":HOSTS_SHA,
      "source_original_HOSTS_record_count":total,
      "source_Azores_2023_original_archive_sha":official["original_archive_sha256"],
      "original_n_Azores_butterfly_species":8,
      "source_official_native_nonendemic":6,
      "source_official_endemic":2,
      "original_HOSTS_references_headers":reference_fields,
      "source_records":{
          name:{
            "Azores_origin":source_status[name],
            "raw_insect_host_rows":len(index[name]),
            "named_plant_species":sorted({x["host_binomial"] for x in index[name]
                                         if x["host_genus"] and x["host_epithet"]}),
            "raw_rows":index[name]
          } for name in sorted(NAMES)},
      "warning":[
        "A host plant species recorded for a butterfly in NHM HOSTS does NOT guarantee this plant occurs on or is used on each Azores island.",
        "Original HOSTS raw evidence does not establish confidence, experimental larval performance or local native-versus-alien plant origin without independently checking the source metadata.",
        "The two island-endemic butterflies are absent from the original GEB 339-source species manifest before final 239; reason can be the major LepTraits environmental-affinity requirements, not necessarily a missing host."
      ]
    }
    (output/"Azores_original_NHM_HOSTS_eight_source_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print("AZORES_NATIVE_ENDEMIC_ORIGINAL_HOSTS_REPLAY="+json.dumps({
        "host_total_raw":total,
        "endemic_host_records":{name:summary["source_records"][name]["named_plant_species"] for name in sorted(AZORES_END)},
        "source_original_host_reference_fields":reference_fields,
        "source_official_status":source_status,
    },ensure_ascii=False),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts",required=True,type=Path)
    p.add_argument("--azores-source-summary",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    args=p.parse_args()
    validate(args.hosts,args.azores_source_summary,args.out_dir)

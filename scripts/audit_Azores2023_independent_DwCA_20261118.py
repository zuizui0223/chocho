#!/usr/bin/env python3
"""Source-lock 2023 Azores independent arthropod checklist and inspect butterfly names.

Official University of Azores GBIF IPT Darwin Core Archive v1.7:
https://doi.org/10.15468/vzpxhx
Do not derive in-situ status from HOSTS or GBIF butterfly occurrences.
Core Taxon and Distribution extension are parsed through META.XML field indices.
Only exact independent INBO binomials are counted in the initial join;
unmatched island endemics and source taxonomy conflicts remain explicit.
"""
import argparse,csv,hashlib,io,json,re,time,urllib.request,urllib.error,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

SOURCE="https://ipt.gbif.pt/ipt/archive.do?r=checklist_arthropoda_azores&v=1.7"
DOI="10.15468/vzpxhx"
INBO=Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117/independent_42_species_vs_frozen239.csv")
DWC="http://rs.tdwg.org/dwc/terms/"
REQUIRED_FIELDS={"Taxon":{"taxonID","scientificName","family"},
                 "Distribution":{"taxonID","locality","establishmentMeans"}}

def src_fetch(out,offline=False):
    dest=out/"azores_arthropoda_checklist_v1_7_dwca.zip"
    if dest.exists() and dest.stat().st_size>1000:
        return dest
    if offline:raise FileNotFoundError("Missing original source archive for offline replay")
    req=urllib.request.Request(SOURCE,headers={"User-Agent":"chocho-source-independent-island-checklist-scientific-audit/2026-11-18",
      "Accept":"application/zip,application/octet-stream,*/*"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=50) as f:
                payload=f.read(8_000_000)
                if f.status!=200:raise RuntimeError("Archive source HTTP status "+str(f.status))
            if not payload.startswith(b"PK") or len(payload)>7_500_000:
                raise ValueError("Original IPТ source is not appropriately sized ZIP")
            dest.write_bytes(payload)
            return dest
        except (OSError,urllib.error.HTTPError,urllib.error.URLError,ValueError):
            if attempt==2:raise
            time.sleep(2*(attempt+1))

def inspect_meta(z):
    names=z.namelist()
    opts=[x for x in names if x.endswith("meta.xml")]
    if len(opts)!=1:raise ValueError("Expected one Darwin Core Archive meta.xml")
    root=ET.fromstring(z.read(opts[0]))
    structs={}
    for child in root:
        kind=child.tag.split("}")[-1]
        if kind not in ("core","extension"):continue
        txt=child.find("./{*}files/{*}location")
        if txt is None or not txt.text:raise ValueError("Meta file location missing")
        uri=child.attrib.get("rowType","")
        label="Taxon" if "Taxon" in uri else "Distribution" if "Distribution" in uri else kind
        cols={}
        for c in child:
            if c.tag.split("}")[-1]=="field":
                term=c.attrib.get("term","")
                pos=c.attrib.get("index","")
                if pos!="":
                    cols[term.rsplit("/",1)[-1]]=int(pos)
        ident=child.find("./{*}id")
        if ident is not None and ident.attrib.get("index") is not None:
            cols["taxonID"]=int(ident.attrib["index"])
        # The independent source may use core taxon ID and extension FK coreid
        if kind=="extension":
            coreid=child.find("./{*}coreid")
            if coreid is not None and coreid.attrib.get("index") is not None:
                cols["_coreid"]=int(coreid.attrib["index"])
        structs[label]={"file":txt.text.strip(),
            "delimiter":child.attrib.get("fieldsTerminatedBy","\t"),
            "quote":child.attrib.get("fieldsEnclosedBy",'"'),
            "ignore_header":int(child.attrib.get("ignoreHeaderLines","0")),
            "columns":cols,"rowType":uri}
    if "Taxon" not in structs or "Distribution" not in structs:
        raise ValueError("Missing original authoritative core/extension metadata")
    return structs

def safe_token(v):
    # Darwin Core 2017 text names store species authorities; only species binomial
    tokens=re.findall(r"[A-Za-z][A-Za-z.-]*",str(v or ""))
    if len(tokens)<2 or not re.fullmatch(r"[A-Z][a-zA-Z-]+",tokens[0]) or not re.fullmatch(r"[a-z][a-z-]+",tokens[1]):
        return ""
    return tokens[0]+" "+tokens[1]

def table(z,spec):
    filename=spec["file"]
    files=z.namelist()
    if filename not in files:
        matching=[x for x in files if x.endswith("/"+filename)]
        if len(matching)!=1:raise ValueError("Source meta file missing: "+filename)
        filename=matching[0]
    text=z.read(filename).decode("utf-8-sig")
    delimiter=spec["delimiter"]
    if delimiter=="\\t":delimiter="\t"
    if len(delimiter)!=1:raise ValueError("Unsupported original delimiter: "+repr(delimiter))
    f=io.StringIO(text,newline="")
    reader=csv.reader(f,delimiter=delimiter,quotechar=spec["quote"] or '"')
    for _ in range(spec["ignore_header"]):next(reader,None)
    rows=[]
    keys=spec["columns"]
    for values in reader:
        item={field:(values[pos] if pos<len(values) else "") for field,pos in keys.items()}
        rows.append(item)
    return rows

def main(dest,offline=False):
    dest.mkdir(parents=True,exist_ok=True)
    path=src_fetch(dest,offline)
    with zipfile.ZipFile(path) as z:
        structs=inspect_meta(z)
        taxa=table(z,structs["Taxon"])
        distribution=table(z,structs["Distribution"])
    source_names={r["scientific_name"]:r for r in list(csv.DictReader(INBO.open(encoding="utf-8-sig",newline=""))) if "Azores" in r["archipelagos"]}
    if len(source_names)!=8:
        raise ValueError("Original exact INBO Azores 8 species changed")
    if len(taxa)<4000 or len(distribution)<8000:
        raise ValueError("Unexpectedly incomplete official original 2023 taxon/distribution dataset")
    if not all(s in structs["Taxon"]["columns"] for s in ("taxonID","scientificName","family")):
        raise ValueError("Missing core names or family")
    if not all(s in structs["Distribution"]["columns"] for s in ("locality","establishmentMeans")):
        raise ValueError("Missing original independently classified establishment status")
    # Separate official taxon identity key; archive may refer to core-id only.
    idmap={r["taxonID"]:r for r in taxa if r.get("taxonID")}
    # Use actual source taxon ID, NEVER infer as offset/index without verified core ID
    source_by_name={}
    for row in taxa:
        nm=safe_token(row.get("scientificName",""))
        # Species-rank only: never infer endemism of the whole Pieris brassicae
        # species from an endemic named subspecies row in the official source.
        if (nm in source_names and row.get("family","").lower()==source_names[nm]["checklist_family"].lower()
            and row.get("taxonRank","").strip().lower()=="species"):
            source_by_name.setdefault(nm,[]).append(row)
    exact_join=[]
    for name, source in sorted(source_names.items()):
        matched=source_by_name.get(name,[])
        ids={r["taxonID"] for r in matched}
        dist=[r for r in distribution if r.get("taxonID") in ids or r.get("_coreid") in ids]
        status=sorted({r.get("establishmentMeans","").strip() for r in dist if r.get("establishmentMeans","").strip()})
        exact_join.append({"independent_INBO_species":name,"family":source["checklist_family"],
               "source_taxon_rows":len(matched),"source_taxonIDs":sorted(ids),
               "source_distribution_rows":len(dist),"official_establishment_means":status,
               "source_localities":sorted({r.get("locality","") for r in dist if r.get("locality")}),
               "original_HOSTS_named_genus_exact":int(source["hosts_named_genus"]),
               "in_original_GEB239_binomial_exact":int(source["in_frozen239"])})
    taxon_concept_diagnostics=[
      {"taxonID":r.get("taxonID"),"scientificName":r.get("scientificName"),
       "acceptedNameUsage":r.get("acceptedNameUsage"),
       "taxonRank":r.get("taxonRank"),"taxonomicStatus":r.get("taxonomicStatus"),
       "family":r.get("family"),"infraspecificEpithet":r.get("infraspecificEpithet")}
      for r in taxa if r.get("family","").lower() in ("pieridae","nymphalidae") and (
          r.get("scientificName","").startswith("Colias croce") or
          r.get("scientificName","").startswith("Hipparchia azor") or
          r.get("scientificName","").startswith("Hipparchia miguel") or
          r.get("scientificName","").startswith("Pieris brassicae")
      )]
    # Only one homotypic species equivalence was independently verified
    # against the original GEB239 source; never auto-collapse azorina/azorinus.
    guarded_verified_alias={"Colias crocea":"Colias croceus"}
    species_rank_by_name={}
    for r in taxa:
        if r.get("taxonRank","").strip().lower()!="species":continue
        nm=safe_token(r.get("scientificName",""))
        species_rank_by_name.setdefault(nm,[]).append(r)
    for row in exact_join:
        name=row["independent_INBO_species"]
        source_alias=guarded_verified_alias.get(name,"")
        if row["source_taxon_rows"]==0 and source_alias:
            alt=[x for x in species_rank_by_name.get(source_alias,[])
                  if x.get("family","").lower()==row["family"].lower()]
            row["independently_verified_alias"]=source_alias if len(alt)==1 else ""
            row["alias_taxon_rank_only"]=len(alt)
            row["alias_official_origin_status"]=[]
            if len(alt)==1:
                rr=[x for x in distribution if x.get("_coreid")==alt[0]["taxonID"]]
                row["alias_official_origin_status"]=sorted({x["establishmentMeans"] for x in rr if x.get("establishmentMeans")})
                row["alias_official_named_islands"]=sorted({x.get("locality","") for x in rr if x.get("locality")})
        if name=="Hipparchia azorina":
            row["taxonomy_conflict_note"]="Official 2023 source may use H. azorinus and infraspecific names; not automatically counted as synonym of H. azorina."
    report={
       "schema":"chocho_azores_official2023_independent_origin_initial_crosswalk_v1",
       "source":SOURCE,"source_doi":DOI,"source_description":"Borges et al 2023 Updated Checklist of Arthropods from Azores v1.7 University Azores",
       "original_archive_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
       "source_archive_bytes":path.stat().st_size,
       "meta":structs,
       "counts":{"taxa":len(taxa),"distribution":len(distribution),"original_INBO_Azores_unique":len(source_names)},
       "source_taxon_sample":taxa[:4],
       "source_distribution_sample":distribution[:4],
       "original_source_taxon_concept_diagnostics":taxon_concept_diagnostics,
       "exact_taxon_results":exact_join,
       "known_limits":[
         "Only source-exact original butterfly binomials are joined; accepted synonyms and subspecies are a separate, explicit subsequent gate.",
         "Source distribution establishmentMeans must be read at actual named island or archipelago scope, not silently assigned to unrelated butterfly islands.",
         "An exact source match for butterfly names does not prove its larval host or success."
       ]
    }
    (dest/"azores_official_2023_SOURCE_AUDIT.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("AZORES_OFFICIAL2023_DWCA_AUDIT="+json.dumps({
       "archive_sha":report["original_archive_sha256"],
       "counts":report["counts"],
       "taxon_headers":structs["Taxon"]["columns"],
       "dist_headers":structs["Distribution"]["columns"],
       "independent_eight":exact_join,
       "taxonomy_source_conflicts":taxon_concept_diagnostics
    },ensure_ascii=False)[:12000],flush=True)
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",required=True,type=Path)
    p.add_argument("--offline",action="store_true")
    args=p.parse_args()
    main(args.out_dir,args.offline)

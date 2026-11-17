#!/usr/bin/env python3
"""2006 Sri Lanka explicit-E endemic butterflies versus fixed HOSTS-selected 225 names.

The 2006 source says 20 historical endemics but 18 exact, explicitly E-marked
appendix binomials were transcribed. This is NOT a 2026 taxonomic-endemism
classification or modern accepted-name reconciliation. Two additional endemic
species in source summary remain unassigned pending original tabular audit.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE="https://portals.iucn.org/library/efiles/documents/2006-030.pdf"
COHORT="data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
FROZEN_BLOB="b6a3c98eb76864052d45d520a55b7d830816126c"
E_LIST="results/butterfly_GBIF_effort_origin_gate_20261117/SRILANKA_2006_IUCN_18_EXPLICIT_E_ENDEMICS_20261117.csv"
E_BLOB="a1ed7d50be385fd9fcb99ef9660bda0ec9919a91"
HISTORICAL_REPORTED_TOTAL_SPECIES=243
HISTORICAL_REPORTED_ENDEMICS=20
FAMILY_REPORT={"Papilionidae":2,"Nymphalidae":8,"Lycaenidae":8,"Hesperiidae":2}

def blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode("ascii")+b"\0"+b).hexdigest()

def source(path,expected):
    if blob(path)!=expected:
        raise ValueError("Source Git blob mismatch: "+str(path))
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def analyze():
    rows=source(ROOT/COHORT,FROZEN_BLOB)
    labels=source(ROOT/E_LIST,E_BLOB)
    names={r["species"].strip() for r in rows}
    ids={r["speciesKey"] for r in rows}
    if (len(rows),len(names),len(ids))!=(11475,225,225):
        raise ValueError("Original 225 consumer names/species keys have changed")
    if len(labels)!=18 or len({r["historical_binomial_2006"] for r in labels})!=18:
        raise ValueError("Expected unique 18 explicit E-marked names from published appendix")
    family=Counter(r["historical_family"] for r in labels)
    if family != {"Papilionidae":2,"Nymphalidae":6,"Lycaenidae":8,"Hesperiidae":2}:
        raise ValueError("Explicit-E transcription family counts have changed")
    output=[]
    for row in labels:
        n=row["historical_binomial_2006"].strip()
        year=int(row["source_year"])
        page=int(row["butterfly_appendix_printed_page"])
        if not 59<=page<=64 or year!=2006 or row["source_url"]!=SOURCE:
            raise ValueError("Original source page/year or URL lost")
        if row["endemic_code_as_extracted"] not in ("separate E","attached E"):
            raise ValueError("Label does not contain explicitly transcribed E marking")
        if row["current_taxonomic_status"]!="NOT_RECONCILED":
            raise ValueError("Historical E label must not become modern native or endemic status")
        output.append(dict(historical_binomial=n,family=row["historical_family"],
                           printed_page=page,source_year=year,
                           match_original_225_exact_binomial=n in names,
                           modern_synonymy_reconciled=False,
                           current_native_endemic_classification="UNASSESSED"))
    matches=[z["historical_binomial"] for z in output if z["match_original_225_exact_binomial"]]
    if matches: raise ValueError("Historical exact E-marked binomial unexpectedly represented: "+str(matches))
    return {
      "schema":"chocho_2006_IUCN_SriLanka_18_explicit_historical_E_marks_vs_original225_v1",
      "status":"HISTORICAL_BINOMIAL_SOURCE_SELECTION_BIAS_NOT_MODERN_TAXONOMY_OR_2026_ENDENISM",
      "evidence":{"source_pdf":SOURCE,"printed_appendix_pages":"59-64",
                  "2006_source_reported_fauna_species":HISTORICAL_REPORTED_TOTAL_SPECIES,
                  "2006_source_prose_reported_endemics":HISTORICAL_REPORTED_ENDEMICS,
                  "appendix_explicit_E_marks_transcribed":len(labels),
                  "reported_total_vs_explicit_E_marks_difference":HISTORICAL_REPORTED_ENDEMICS-len(labels),
                  "difference_unresolved":"Two historic endemic species under the prose definition not recovered as explicit E-marked binomials from the text transcription; cannot safely assign their names or reconcile 2026 taxonomy.",
                  "prose_family_distribution":FAMILY_REPORT,
                  "explicit_E_family_distribution":dict(sorted(family.items())),
                  "E_legend":"Appendix p64: E - Endemic species"},
      "original_selected_consumer_panel":{"species":len(names),"species_island_cells":len(rows),"source_git_blob":FROZEN_BLOB},
      "exact_historical_name_overlap":{"2006_explicitly_E_marked_names":len(labels),
                                        "original225_exact_binomial_matches":len(matches),
                                        "not_exactly_named_in_original225":len(labels)-len(matches),
                                        "present":matches},
      "historical_endemic_name_table":output,
      "interpretation":["Complete absence of exact matching historical E-marked names in selected 225 is direct sampling-frame/taxonomic-name coverage evidence for Sri Lanka, NOT that butterflies lack endemics on islands.",
                        "This is not a valid estimate of percentage of 2026 biologically endemic butterflies missed: remaining 2 historic taxa, synonyms, revised species concepts and current ranges unresolved.",
                        "The IUCN 2006 appendix uses full Sri Lanka national-island scope, not exact observed WPI port buffers.",
                        "Do not infer the 2006 non-E names are all nonendemic; summary and E-flag counts do not align.",
                        "Presence/non-presence in a pre-selected HOSTS-compatible butterfly set has a different denominator than Sri Lanka fauna species richness."]
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,default=Path("_srilanka_historical_endemic_selection/E_MARKED_2006_VS_225_20261117.json"))
    a=ap.parse_args()
    x=analyze()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:x[k] for k in ("evidence","original_selected_consumer_panel","exact_historical_name_overlap")},indent=2))
if __name__=="__main__":
    main()

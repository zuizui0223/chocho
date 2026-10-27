#!/usr/bin/env python3
"""Real Lepidoptera interceptions across original South African import commodities.

Source: Saccaggi et al. Scientific Data 2021, DOI 10.1038/s41597-021-00869-z,
original Figshare COLLECTION 10.6084/m9.figshare.c.5180681.
Original CSV sample 26,295 rows, 92 fields; paper cleaned subset 25,279
NOT recoverable by simply calling all raw source lines "study samples".
Original 878 taxonomic metadata rows identify historical and current status.
First sheet 'Insecta: Lepidoptera': NA=UNTESTED, 0=TESTED NEGATIVE,
other nonempty=POSITIVE taxonomic LABEL. No stage or wild establishment
can be inferred for unlabeled insects. Source is ALL Lepidoptera, not butterflies.
No individual consignment ID, importer/exporter code or raw row exported.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlparse

from audit_southafrica_plant_import_figshare_source_20261027 import source_get
from analyze_hosts_location_wcvp_origin_20261027 import read_rows

SCHEMA="chocho_SA_original_plant_import_lepidoptera_interception_routes_20261027_v01"
URLS={
 "sample":"https://ndownloader.figshare.com/files/25166393",
 "metadata":"https://ndownloader.figshare.com/files/26235020",
}
SOURCE_SHA={
 "sample":"c12d4af0c96a247d520f71ba6c2dfa67eacb33b53dd5aa329265c850d3fe15da",
 "metadata":"0f7ad1b9b44dad51a53d55a1f6bad193d623b0b8beb7853eb0b6f025dfcf640c",
}
NEEDED=("catalogNumber","year","sample.source","crop.commodity","Insecta: Lepidoptera")
REGIMES=("audit","intervention","PEQ","MISSING_REGIME","UNREGISTERED_REGIME")
PERIODS=("1994-2004","2005-2019","MISSING_OR_OTHER_YEAR")
CLASSES=("LIVE_PROPAGATION","SEED_GRAIN","FRESH_OR_DRY_PRODUCE","OTHER_OR_UNRESOLVED")
NATIVE_PANEL="data/frozen/figure_sources/anthropogenic_species_metrics.csv"
EXPECTED_STATUS={"TESTED_NEGATIVE":15673,"NOT_TESTED":10443,"DETECTED_LEPIDOPTERA":178,
                 "MISSING_ASSESSMENT":1}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def load_primary_file(alias):
    content,actual,ctype=source_get(URLS[alias])
    if sha(content)!=SOURCE_SHA[alias]:
        raise ValueError("Public original Figshare source byte SHA changed: "+alias)
    try:
        source_text=content.decode("utf-8-sig")
        decoding="utf-8-sig"
    except UnicodeDecodeError:
        source_text=content.decode("latin-1")
        decoding="latin-1 byte-preserving fallback (original byte SHA retained)"
    rows=list(csv.DictReader(io.StringIO(source_text),delimiter=";"))
    if alias=="sample" and len(rows)!=26295:
        raise ValueError("Original raw laboratory source 26,295-row grain changed")
    if alias=="metadata" and len(rows)!=878:
        raise ValueError("Original taxon metadata 878-row grain changed")
    return rows,{"source_download":URLS[alias],"member_SHA256":sha(content),
                 "source_bytes":len(content),"original_source_rows":len(rows),
                 "field_names":list(rows[0]),"original_text_decode":decoding}


def lepidoptera_result(value):
    val=str(value or "").strip()
    if val.upper()=="NA":
        return "NOT_TESTED"
    if val=="0":
        return "TESTED_NEGATIVE"
    if val=="":
        return "MISSING_ASSESSMENT"
    return "DETECTED_LEPIDOPTERA"


def commodity_class(raw):
    """Transported material form is not proof of hosted insect/egg/larva."""
    v=str(raw or "").strip()
    if not v:
        return "OTHER_OR_UNRESOLVED"
    low=v.lower()
    if (low.startswith(("cuttings","cutting","in vitro","plantlets",
                        "live plant","rooted plant")) or low in {
        "plants","bulbs","rhizomes","tubers","roots","plantlets","cuttings",
        "tissue culture","planting material"}):
        return "LIVE_PROPAGATION"
    if low in {"seeds","seed","grains","grain","beans","bean","nuts","nut"}:
        return "SEED_GRAIN"
    if (low.startswith(("fresh ","dried ","dry ","frozen ")) or
         low in {"fruit","vegetables","foliage","cut flowers","flower buds",
                 "dried fruit"}):
        return "FRESH_OR_DRY_PRODUCE"
    return "OTHER_OR_UNRESOLVED"


def sample_period(value):
    try:
        yr=int(float(str(value)))
    except (ValueError,TypeError):
        return "MISSING_OR_OTHER_YEAR"
    if 1994<=yr<=2004:return "1994-2004"
    if 2005<=yr<=2019:return "2005-2019"
    return "MISSING_OR_OTHER_YEAR"


def detection_count(count):
    pos=count.get("DETECTED_LEPIDOPTERA",0)
    neg=count.get("TESTED_NEGATIVE",0)
    unt=count.get("NOT_TESTED",0)
    missing=count.get("MISSING_ASSESSMENT",0)
    n=pos+neg
    return {
        "original_raw_samples":pos+neg+unt+missing,
        "lepidoptera_tested_samples":n,
        "detected_Lepidoptera_sample_rows":pos,
        "original_tested_negative_samples":neg,
        "lepidoptera_NOT_TESTED_samples":unt,
        "missing_source_lepidoptera_assessment":missing,
        "conditional_detection_per_100_TESTED":100*pos/n if n else None,
    }


def aggregate(samples):
    counts=Counter()
    per_source=defaultdict(Counter)
    per_commodity=defaultdict(Counter)
    per_type=defaultdict(Counter)
    per_joint=defaultdict(Counter)
    per_joint_period=defaultdict(Counter)
    original_commodity=Counter()
    original_positive_names=Counter()
    per_year=defaultdict(Counter)
    original_sample_ids=set()
    for row in samples:
        if not set(NEEDED)<=set(row):
            raise ValueError("Original Figshare 92-column sample metadata changed")
        source_val=(row["sample.source"] or "").strip()
        regime=source_val if source_val in REGIMES else (
            "MISSING_REGIME" if not source_val else "UNREGISTERED_REGIME")
        commodity=(row["crop.commodity"] or "").strip() or "<EMPTY>"
        typ=commodity_class(commodity)
        v=(row["Insecta: Lepidoptera"] or "").strip()
        status=lepidoptera_result(v)
        period=sample_period(row.get("year"))
        orig_year=(row.get("year") or "").strip() or "<MISSING>"
        counts[status]+=1
        per_source[regime][status]+=1
        per_commodity[commodity][status]+=1
        per_type[typ][status]+=1
        per_joint[(regime,typ)][status]+=1
        per_joint_period[(regime,typ,period)][status]+=1
        per_year[orig_year][status]+=1
        original_commodity[commodity]+=1
        if status=="DETECTED_LEPIDOPTERA":
            original_positive_names[v]+=1
        ident=(row.get("catalogNumber") or "").strip()
        if ident and ident in original_sample_ids:
            counts["duplicated_nonblank_catalogNumber"]+=1
        if ident:original_sample_ids.add(ident)
    if {k:counts[k] for k in EXPECTED_STATUS}!=EXPECTED_STATUS:
        raise ValueError("Original public raw Lepidoptera test NA/0/positive source margins changed")
    if sum(original_commodity.values())!=26295 or len(original_commodity)!=31:
        raise ValueError("Original 30 commodity labels plus 1 missing changed")
    def output_counts(d):
        return [dict({"group":group},**detection_count(v)) for group,v in
                sorted(d.items(),key=lambda x:str(x[0]))]
    return {
        "all_source_rows":detection_count(counts),
        "source_duplicate_nonblank_catalog_number_counts":
            counts["duplicated_nonblank_catalogNumber"],
        "by_sample_source":output_counts(per_source),
        "by_original_commodity":output_counts(per_commodity),
        "by_predeclared_material_class":output_counts(per_type),
        "by_sampling_regime_and_material_class":[
            {"sampling_regime":s,"material_class":c,**detection_count(v)}
            for (s,c),v in sorted(per_joint.items())
        ],
        "by_regime_class_and_era":[
            {"sampling_regime":s,"material_class":c,"era":period,
             **detection_count(v)}
            for (s,c,period),v in sorted(per_joint_period.items())
        ],
        "by_original_source_year":output_counts(per_year),
        "source_commodity_value_and_predeclared_class":[
            {"original_commodity_label":raw,
             "predeclared_material_class":commodity_class(raw),
             **detection_count(per_commodity[raw])}
            for raw in sorted(per_commodity)
        ],
        "all_original_distinct_lepidoptera_positive_labels":len(original_positive_names),
        "original_positive_Lepidoptera_literal_labels":[
            {"original_reported_taxon_stage_label":name,
             "positive_sample_rows":n} for name,n
            in original_positive_names.most_common()
        ],
        "year_2005_registration_change_paper":True,
        "reported_paper_FINAL_cleaned_samples":25279,
        "note_original_raw_26295_not_identical_to_paper_final_25279":True,
        "biological_precision":"Order Lepidoptera moths and butterflies; some labels are larval, many species/higher taxa, no direct linkage to wild establishment.",
    }


def taxon_metadata(rows, positive_labels, butterfly_names):
    if not {"order","scientificName","SA.occurenceStatus.historical",
            "SA.occurenceStatus.current","Taxon"}<=set(rows[0]):
        raise ValueError("Original 14-column source taxon metadata changed")
    lepid=[r for r in rows if (r.get("order") or "").strip()=="Lepidoptera"]
    names=Counter(str(r.get("scientificName") or "").strip() for r in lepid)
    historic=Counter((r.get("SA.occurenceStatus.historical") or "").strip() for r in lepid)
    present=Counter((r.get("SA.occurenceStatus.current") or "").strip() for r in lepid)
    later=[r for r in lepid if (r.get("SA.occurenceStatus.historical") or "").strip().lower()=="absent" and
           (r.get("SA.occurenceStatus.current") or "").strip().lower().startswith("present")]
    detected_exact_binomials=set()
    ambiguous=0
    # Original taxon label free text may include family and cf.; keep exact
    # two-token prefix only when Latin binomial canonical and unambiguous.
    for raw in positive_labels:
        text=raw.strip()
        if not text or "cf." in text or "aff." in text:
            ambiguous+=1
            continue
        m=re.fullmatch(r"([A-Z][a-z]+ [a-z][a-z-]+)(?: \([A-Z][A-Za-z\- ]+\))?",text)
        if not m or m.group(1).startswith(("Lepidoptera ","Pyralidae ","Tortricidae ")):
            ambiguous+=1
            continue
        detected_exact_binomials.add(m.group(1))
    original_intersection=sorted(detected_exact_binomials & butterfly_names)
    lepid_present=sorted(set(names)&detected_exact_binomials)
    return {
        "total_original_taxon_metadata_rows":len(rows),
        "Lepidoptera_metadata_taxa_rows":len(lepid),
        "source_original_taxonomic_metadata_distinct_names":len(names),
        "Lepidoptera_historical_SA_status":dict(historic),
        "Lepidoptera_current_SA_status":dict(present),
        "source_later_established_Lepidoptera_rows":{
            "reported_historically_absent_and_later_present_count":len(later),
            "observed_taxa_names_from_original_metadata":[
                r.get("scientificName","") for r in later],
            "note":"Historical-to-current status changed in metadata, not proof interception caused establishment."
        },
        "unambiguous_original_source_positive_binomials":len(detected_exact_binomials),
        "ambiguous_or_higher_taxon_source_positive_labels":ambiguous,
        "unambiguous_positive_names_matching_original_metadata":lepid_present,
        "exact_positive_species_names_matching_frozen_239_butterflies":original_intersection,
        "butterfly_taxon_population_frozen":len(butterfly_names),
    }


def audit(source_table,metadata_table,original239_source):
    orig=read_rows(original239_source)
    butterfly={r["species"] for r in orig}
    if len(orig)!=239 or len(butterfly)!=239:
        raise ValueError("The original 239 butterfly source cohort changed")
    grouped=aggregate(source_table)
    metadata=taxon_metadata(
        metadata_table,
        {r["original_reported_taxon_stage_label"] for r in
         grouped["original_positive_Lepidoptera_literal_labels"]},
        butterfly
    )
    # Compare live / seeds only in EXACT identical sampling regimes; do not
    # present PEQ live as 'effect' against audit seeds.
    matched=[]
    for regime in REGIMES:
        original={r["material_class"]:r for r in grouped[
            "by_sampling_regime_and_material_class"] if
                  r["sampling_regime"]==regime}
        a=original.get("LIVE_PROPAGATION",{})
        b=original.get("SEED_GRAIN",{})
        ta=a.get("lepidoptera_tested_samples",0)
        tb=b.get("lepidoptera_tested_samples",0)
        pos=a.get("detected_Lepidoptera_sample_rows",0)+b.get("detected_Lepidoptera_sample_rows",0)
        status=("SOURCE_DIAGNOSTIC_ONLY_PRESPECIFIED_INFORMATION_GATE" if
                ta>=30 and tb>=30 and pos>=5 else
                "INSUFFICIENT_COMMODITY_OVERLAP_WITHIN_SAMPLING_REGIME")
        matched.append({
            "sampling_regime":regime,
            "live_propagation_tested_samples":ta,
            "seed_grain_tested_samples":tb,
            "live_propagation_positive_samples":a.get("detected_Lepidoptera_sample_rows",0),
            "seed_grain_positive_samples":b.get("detected_Lepidoptera_sample_rows",0),
            "live_minus_seed_detected_percentage_points":(
                100*(a["detected_Lepidoptera_sample_rows"]/ta-
                     b["detected_Lepidoptera_sample_rows"]/tb)
                if ta and tb else None),
            "status":status,
            "not_causal_transport_trade_effect":True,
        })
    return {
        "schema":SCHEMA,
        "status":"SOURCE_COMPLETE_RAW_LAB_SAMPLE_LEPIDOPTERA_COMMODITY_CONDITIONAL_DETECTION",
        "original_source":{
            "article_DOI":"10.1038/s41597-021-00869-z",
            "figshare_COLLECTION_DOI":"10.6084/m9.figshare.c.5180681",
            "original_sample_figshare_file_id":25166393,
            "original_taxon_figshare_file_id":26235020,
            "raw_sample_file_SHA256":SOURCE_SHA["sample"],
            "taxon_metadata_file_SHA256":SOURCE_SHA["metadata"],
            "original_239_butterfly_metric_panel_SHA256":hashlib.sha256(
                Path(original239_source).read_bytes()).hexdigest(),
        },
        "samples":grouped,
        "Lepidoptera_taxon_metadata":metadata,
        "within_identical_sampling_regime_live_vs_seed_diagnostic":matched,
        "published_samples_are_cleaned_subset_not_original_raw":{
            "original_source_raw_sample_rows":26295,
            "published_paper_final_sample_rows":25279,
            "unreplicated_exclusions":1016,
        },
        "limits":[
            "These are laboratory-tested importer specimens, not all imported consignments, not rates per tonne or per international shipment.",
            "PEQ, intervention and audit have different source sampling selection; pooled commodity class associations confounded.",
            "Taxon NA means NOT tested, distinct from actual 0 tested negative.",
            "Main 178 nonzero source rows encompass Lepidoptera moths and butterflies and higher-stage taxon names; butterfly-only stage/country not identified for most.",
            "Crop commodity may carry incidental adult/larval stage and may not itself be the insect's larval host.",
            "Separate taxon metadata later-establishment status does not causally link intercepted individual to established descendants.",
            "The original 239 butterfly source HOSTS/WCVP host resource portfolio cannot be joined directly to these aggregate source positive stage labels for causal arrival-resource estimates.",
            "Original publicly collected unfiltered source includes 26295 rows, the paper's finalized 25279 is a further curated subset. The author's exact exclusion filter has not yet been reproduced.",
            "This derived report excludes all original sample-level identifiers and raw importer/exporter identities."
        ],
        "computed_utc":datetime.now(timezone.utc).isoformat()
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original239",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    samples,s_source=load_primary_file("sample")
    metadata,t_source=load_primary_file("metadata")
    r=audit(samples,metadata,args.original239)
    r["original_source"]["sample_original_retrieval"]=s_source
    r["original_source"]["taxon_original_retrieval"]=t_source
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    x=r["samples"]
    print(json.dumps({
        "status":r["status"],
        "source_all_samples":x["all_source_rows"],
        "by_source":x["by_sample_source"],
        "by_material":x["by_predeclared_material_class"],
        "by_sample_regime_material":x["by_sampling_regime_and_material_class"],
        "prespecified_within_source_compare":r["within_identical_sampling_regime_live_vs_seed_diagnostic"],
        "taxon":r["Lepidoptera_taxon_metadata"],
        "raw_commodity_categories":x["source_commodity_value_and_predeclared_class"],
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()

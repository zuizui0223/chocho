#!/usr/bin/env python3
"""Prior-1987 *recorded* host-opportunity variability at actual UKBMS transects.

Audits spatial IDENTIFIABILITY only, not any butterfly effect. Each host is
documented in original frozen HOSTS source and exact-name BSBI Plant Atlas.
No plant observation is converted to known absence, even at matched hectad.
Original UKBMS sites, botanical hectads, author-only HOSTS edges, and
geocoordinates are transient in the Actions runner, never in the source result.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,statistics,zipfile
from collections import Counter,defaultdict
from pathlib import Path

from audit_bsbi_historical_hectad_source_20261029 import obtain
from audit_bsbi_nested_hectad_layout_20261029 import NESTED
from audit_bsbi_original_hosts_name_match_20261029 import (
    get_host_edges,canonical,AUTHOR_SHA
)
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk,source_visits,evaluate
from audit_uk_ltmn_butterfly_detection_panel_20261028 import CROSSWALK

SCHEMA="chocho_uk_19_LTMN_sites_preatlas1987_doc_host_BSBI_coverage_20261029_v01"
PERIODS_PRIOR=("to 1929","1930 - 1969","1970 - 1986")
PERIODS_LATER=("1987 - 1999","2000 - 2009","2010 - 2019")
EXPECTED_ORIGINAL_TAXA=3497
MAX_INVALID_FLAGS=0

def to_osgb_hectad(easting,northing):
    """Ordnance Survey GB 100km letters plus 10km easting/northing digits."""
    e=int(easting)//100000
    n=int(northing)//100000
    if not (0<=e<=7 and 0<=n<=13):
        raise ValueError("Original UKBMS OSGB easting/northing grid outside GB")
    l1=(19-n)-(19-n)%5+(e+10)//5
    l2=((19-n)*5)%25+e%5
    l1+=int(l1>7);l2+=int(l2>7)
    return (chr(ord('A')+l1)+chr(ord('A')+l2)+
            str((int(easting)//10000)%10)+str((int(northing)//10000)%10))

def recorded_flag(value):
    """Only explicit source positives count; blank/0 not botanical absence."""
    text=str(value or "").strip().casefold()
    if text in {"","0","false","na","n/a","none","null","-"}:return False
    if text in {"1","true","yes","y"}:return True
    try:
        f=float(text)
        if f>=0 and f.is_integer():return f>0
    except ValueError:pass
    raise ValueError("Original BSBI date-class flag outside documented numeric/bool vocabulary")

def plant_lookup(source:bytes):
    with zipfile.ZipFile(io.BytesIO(source)) as z:
        read=z.read("PA2020_10by10kmSquareData/nameLookup.csv")
        inner=z.read(NESTED)
    rows=list(csv.DictReader(io.StringIO(read.decode("utf-8-sig"))))
    if len(rows)!=EXPECTED_ORIGINAL_TAXA:raise ValueError("Original atlas 3497 taxonomy changed")
    maps=defaultdict(set)
    for r in rows:
        key=str(r["ddbidLookup"]).strip()
        for name in (r["taxonName"],r["canonical"]):
            if canonical(name):maps[canonical(name)].add(key)
    return maps,inner

def matched_historical_source(author_artifact,selected_sites,cross):
    if len(selected_sites)!=19:
        raise ValueError("Original LTMN 19 repeat-visit geography changed")
    host_edges,original_hostqc=get_host_edges(author_artifact)
    source,url,errors=obtain()
    maps,inner=plant_lookup(source)
    lookup_for_species={}
    ambiguous=Counter()
    for sp in CROSSWALK:
        matches=set()
        for plant in host_edges.get(sp,set()):
            taxon_ids=maps.get(canonical(plant),set())
            if len(taxon_ids)==1:matches.update(taxon_ids)
            if len(taxon_ids)>1:ambiguous[sp]+=1
        lookup_for_species[sp]=matches
    hectad_by_site={site:to_osgb_hectad(*cross[site][:2]) for site in selected_sites}
    distinct_hectads=set(hectad_by_site.values())
    all_names_in_site={h:set() for h in distinct_hectads}
    all_present_prior=Counter()
    status_values=Counter()
    taxonfile_read=0;nonzero_later=0
    with zipfile.ZipFile(io.BytesIO(inner)) as z:
        for info in z.infolist():
            if not info.filename.endswith(".csv"):continue
            taxon=Path(info.filename).stem
            taxonfile_read+=1
            with z.open(info) as stream:
                reader=csv.DictReader(io.TextIOWrapper(stream,encoding="utf-8-sig",newline=""))
                header=set(reader.fieldnames or [])
                if not {"hectad",*PERIODS_PRIOR,*PERIODS_LATER}<=header:
                    raise ValueError("Original BSBI historic date-class grid headers changed")
                for r in reader:
                    hectad=str(r["hectad"] or "").strip().upper()
                    if hectad not in distinct_hectads:continue
                    prior=any(recorded_flag(r[c]) for c in PERIODS_PRIOR)
                    later=any(recorded_flag(r[c]) for c in PERIODS_LATER)
                    if prior:
                        all_names_in_site[hectad].add(taxon)
                        all_present_prior[taxon]+=1
                    if later:nonzero_later+=1
    if taxonfile_read!=EXPECTED_ORIGINAL_TAXA:
        raise ValueError("Original BSBI 3497 per-taxon historical files incomplete")
    host_support={}
    host_species_with_grid_variation=0
    host_species_with_nonzero_prior=0
    for sp in sorted(CROSSWALK):
        vals=[]
        for site in sorted(selected_sites):
            h=hectad_by_site[site]
            vals.append(len(lookup_for_species[sp] & all_names_in_site[h]))
        if max(vals)>0:host_species_with_nonzero_prior+=1
        if len(set(vals))>1:host_species_with_grid_variation+=1
        host_support[sp]={
            "original_HOSTS_accepted_names":len(host_edges.get(sp,set())),
            "unambiguously_matched_BSBI_named_taxa":len(lookup_for_species[sp]),
            "source_early_documented_host_hectad_count_min":min(vals),
            "source_early_documented_host_hectad_count_max":max(vals),
            "source_early_documented_host_count_median":statistics.median(vals),
            "unique_site_host_coverage_levels":len(set(vals)),
            "sites_with_one_or_more_prior_documented_hosts":sum(x>0 for x in vals),
            "source_n_sites":len(vals),
            "host_nonrecord_is_not_true_absence":True
        }
    effort=[len(all_names_in_site[hectad_by_site[s]]) for s in selected_sites]
    data_support={
        "original_LTMN_repeated_transects":len(selected_sites),
        "unique_BSBI_historical_10km_hectads":len(distinct_hectads),
        "number_original_butterflies_with_prior_documented_host":host_species_with_nonzero_prior,
        "number_butterflies_with_any_site_host_count_variation":host_species_with_grid_variation,
        "atlas_prior_recorded_all_plant_taxa_per_site_min":min(effort),
        "atlas_prior_recorded_all_plant_taxa_per_site_median":statistics.median(effort),
        "atlas_prior_recorded_all_plant_taxa_per_site_max":max(effort),
        "source_time_classes_used_as_historical_predictors":list(PERIODS_PRIOR),
        "source_1987to2019_later_classes_used_in_predictors":False,
        "later_plant_records_inspected_only_for_source_qc":True,
        "BSBI_original_source_sha256":hashlib.sha256(source).hexdigest(),
        "author_original_HOSTS_SHA256":AUTHOR_SHA,
        "n_original_plant_taxa_10km_files":taxonfile_read,
        "source_names_canonical_ambiguous_count":sum(ambiguous.values()),
        "no_historical_host_stem_occupancy_observed":True,
        "no_human_management_exposure_observed":True,
        "no_source_original_site_ids_hectads_plant_occurrence_rows_exported":True,
    }
    if len(distinct_hectads)<10 or host_species_with_grid_variation<5:
        gate="FAIL_PRIOR_HOST_GEOGRAPHY_SUPPORT"
    else:gate="PASS_PRIOR_DOCUMENTED_HOST_GEOGRAPHY_VARIATION_ONLY"
    return {
        "schema":SCHEMA,
        "status":"BSBI_PRIOR_DOCUMENTED_HOST_HECTAD_SITE_COVERAGE_AUDITED",
        "information_gate":gate,
        "original_source":data_support,
        "host_taxa_by_butterfly":host_support,
        "underlying_data_are_presence_only_and_alias_incomplete":True,
        "host_use_true_absence_or_management_not_identified":True,
        "same_19_sites_UKBMS_outcomes_still_not_predictively_modelled":True
    }

def run(author,output):
    try:
        cross,qc=read_crosswalk()
        visits=source_visits()
        original=evaluate(visits,cross)
        if original["strict_year_and_visit_matched_transects"]!=19:
            raise ValueError("Original UKBMS 19 historical repeating route information gate changed")
        sites={s for s,v in visits.items() if
               s in cross and v[0]>=10 and v[1]>=10 and len(v[2])>=3 and len(v[3])>=3}
        report=matched_historical_source(author,sites,cross)
    except (OSError,ValueError,KeyError,zipfile.BadZipFile) as exc:
        report={"schema":SCHEMA,"status":"SOURCE_MATCHING_INCOMPLETE",
                "blocker":str(exc)[:300],
                "no_original_hectad_source_released":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":report["status"],"information_gate":report.get("information_gate"),
        "source":report.get("original_source"),
        "blocker":report.get("blocker")},ensure_ascii=False))
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--author-artifact-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    result=run(a.author_artifact_zip,a.output)
    if result["status"]!="BSBI_PRIOR_DOCUMENTED_HOST_HECTAD_SITE_COVERAGE_AUDITED":raise SystemExit(2)

#!/usr/bin/env python3
"""Join source-frozen UKBMS transects to dated BSBI wild-host evidence.

Original HOSTS/WCVP accepted butterfly-host edges and BSBI 3497 taxon hectad
records are read transiently; only *aggregate* coverage and support gates
leave the source run. An unobserved wild plant is NOT a confirmed absence.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,zipfile
from collections import Counter,defaultdict
from pathlib import Path

from audit_bsbi_historical_hectad_source_20261029 import MD5
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk
from audit_uk_ltmn_historical_event_join_20261028 import get_source,excel_rows,FIELDS_VISITS
from audit_uk_ltmn_NASA_POWER_historical_baseline_20261028 import nasa_site_effort
from audit_uk_ltmn_butterfly_detection_panel_20261028 import CROSSWALK

SCHEMA="chocho_bsbi_exact_historic_wildhost_source_UK_BMS_20261029_v01"
HOST_FILE="inputs/original/insect_host_accepted.csv"
HOST_SHA="0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9"
BSBI_NAME="PA2020_10by10kmSquareData/nameLookup.csv"
BSBI_ZIP="PA2020_10by10kmSquareData/hectadDistributions.zip"
EARLY=("to 1929","1930 - 1969")
PRE1987=EARLY+("1970 - 1986",)
DATE_COLUMNS=PRE1987+("1987 - 1999","2000 - 2009","2010 - 2019")
MIN_SITES=10
MIN_HECTADS=5
MIN_MATCHED_HOST_TAXA=15
MIN_FOCAL_VARIATION=8

def easting_northing_hectad(east,north):
    """OS National Grid two letters + easting/northing 10km digits; no rounding."""
    try:
        east=float(east); north=float(north)
    except (ValueError,TypeError):return None
    if not(math.isfinite(east) and math.isfinite(north) and 0<=east<700000 and 0<=north<1300000):
        return None
    e100=int(east//100000);n100=int(north//100000)
    if e100>6 or n100>12:return None
    # Official OS National Grid scheme skips letter I.
    letters="ABCDEFGHJKLMNOPQRSTUVWXYZ"
    x=19-n100
    first=x-x%5+(e100+10)//5
    second=(x*5)%25+e100%5
    if not 0<=first<25 or not 0<=second<25:return None
    return (letters[first]+letters[second]
            +str(int(east%100000//10000))
            +str(int(north%100000//10000)))

def native_host_edges(source_root):
    path=Path(source_root)/HOST_FILE
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=HOST_SHA:
        raise ValueError("Author-only original accepted HOSTS insect-host source changed")
    links=defaultdict(set)
    with path.open(encoding="utf-8-sig",newline="") as f:
        rows=csv.DictReader(f)
        if not {"insect_species","accepted_name","accepted_plant_name_id"}<=set(rows.fieldnames or []):
            raise ValueError("WCVP accepted host original source columns changed")
        for row in rows:
            species=row["insect_species"].strip()
            accepted=row["accepted_name"].strip()
            if species in CROSSWALK and accepted and len(accepted.split())==2:
                links[species].add(accepted)
    if len(links)!=len(CROSSWALK):
        raise ValueError("Original 27 focal UK butterflies lack accepted source host identities")
    return links

def exact_bsbi_accepted_files(names, filenames):
    """Require unambiguous exact source binomial, species rank, and source ID."""
    identity=defaultdict(set);name_ids=Counter();qc=Counter()
    for row in names:
        name=str(row["taxonName"] or "").strip()
        key=str(row["ddbidLookup"] or "").strip()
        rank=str(row["rank"] or "").strip().lower()
        qc["lookup_rows"]+=1
        if not key or not name:continue
        if (rank not in ("species","sp","sp.") or len(name.split())!=2):
            qc["nonexact_species_or_rank"]+=1;continue
        filename=f"hectadDistributions/{key}.csv"
        if filename in filenames:
            identity[name].add(filename)
            name_ids[name]+=1
        else:qc["lookup_missing_exact_taxon_csv"]+=1
    single={sp:next(iter(files)) for sp,files in identity.items() if len(files)==1}
    ambiguous={sp for sp,files in identity.items() if len(files)>1}
    qc["unique_exact_species_plant_names"]=len(single)
    qc["ambiguous_exact_species_names"]=len(ambiguous)
    return single,dict(qc)

def site_historical_host_coverage(source_root, original_bsbi, *, return_internal=False):
    exact_host=native_host_edges(source_root)
    cross,cross_qc=read_crosswalk()
    all_events,sites,visits=nasa_site_effort(
        excel_rows(get_source("visits"),FIELDS_VISITS),cross)
    if not 15<=len(sites)<=19:
        raise ValueError(f"Original NASA 1981–1995 and later survey cohort changed: {len(sites)}")
    hectad_by_site={site:easting_northing_hectad(*cross[site][:2]) for site in sites}
    if any(h is None for h in hectad_by_site.values()):
        raise ValueError("Original transect BNG 10km grid failed")
    target_hectads=set(hectad_by_site.values())
    raw=Path(original_bsbi).read_bytes()
    if hashlib.md5(raw).hexdigest()!=MD5:
        raise ValueError("Original BSBI ZIP no longer matches published source")
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        with outer.open(BSBI_NAME) as handle:
            names=list(csv.DictReader(io.TextIOWrapper(handle,encoding="utf-8-sig",newline="")))
        nested=outer.read(BSBI_ZIP)
    if len(names)!=3497:raise ValueError("Original 3497 BSBI source lookup names changed")
    with zipfile.ZipFile(io.BytesIO(nested)) as inner:
        files=set(inner.namelist())
        if len([f for f in files if f.endswith(".csv")])!=3497:
            raise ValueError("Original BSBI 3497 taxon distribution files changed")
        exact_by_name,lookup_qc=exact_bsbi_accepted_files(names,files)
        all_focal_hosts=set().union(*exact_host.values())
        matched_names=all_focal_hosts & set(exact_by_name)
        matched_csv={exact_by_name[n]:n for n in matched_names}
        # Botanical survey effort baseline is counted from ALL taxa regardless
        # of whether they are recorded larval hosts for a focal butterfly.
        effort_pre1970=Counter();effort_pre1987=Counter()
        host_presence=defaultdict(lambda:{"early":set(),"pre1987":set()})
        qc=Counter()
        for item in inner.infolist():
            if item.is_dir() or not item.filename.endswith(".csv"):continue
            source_name=matched_csv.get(item.filename)
            with inner.open(item) as f:
                rd=csv.DictReader(io.TextIOWrapper(f,encoding="utf-8-sig",newline=""))
                if not {"hectad",*DATE_COLUMNS}<=set(rd.fieldnames or []):
                    raise ValueError("Original BSBI date columns changed")
                found=set()
                for row in rd:
                    hectad=row["hectad"].strip()
                    if hectad not in target_hectads:continue
                    if hectad in found:raise ValueError("Original plant file duplicate hectad")
                    found.add(hectad)
                    qc["focal_hectad_taxon_source_rows"]+=1
                    if any(str(row[p]).strip() not in ("0","1") for p in DATE_COLUMNS):
                        raise ValueError("Unexpected original BSBI historical period code")
                    early=any(row[p].strip()=="1" for p in EARLY)
                    pre1987=any(row[p].strip()=="1" for p in PRE1987)
                    if early:effort_pre1970[hectad]+=1
                    if pre1987:effort_pre1987[hectad]+=1
                    if source_name is not None:
                        if early:host_presence[hectad]["early"].add(source_name)
                        if pre1987:host_presence[hectad]["pre1987"].add(source_name)
                qc["original_bsbi_taxon_csv_files_scanned"]+=1
                if source_name:qc["original_matched_host_taxon_csv_files_scanned"]+=1
    per_species={}
    variation=0
    for species in sorted(CROSSWALK):
        hosts=exact_host[species]
        source_hosts=hosts & matched_names
        counts_pre1970=[
            len(host_presence[hectad_by_site[site]]["early"] & source_hosts)
            for site in sorted(sites)]
        counts_pre1987=[
            len(host_presence[hectad_by_site[site]]["pre1987"] & source_hosts)
            for site in sorted(sites)]
        historical_var=len(set(counts_pre1970))>1
        if historical_var:variation+=1
        per_species[species]={
            "original_accepted_host_species":len(hosts),
            "exact_wcvp_to_bsbi_taxa":len(source_hosts),
            "transects_with_any_original_pre1970_documented_host":sum(x>0 for x in counts_pre1970),
            "transects_with_any_original_pre1987_documented_host":sum(x>0 for x in counts_pre1987),
            "pre1970_documented_host_count_variation_between_transects":historical_var,
            "min_pre1970_documented_hosts_per_transect":min(counts_pre1970),
            "max_pre1970_documented_hosts_per_transect":max(counts_pre1970),
            "source_past1987_date_classes_not_used_as_pre1970_evidence":True,
            "zero_recorded_host_does_not_mean_host_absent":True,
        }
    eligible=(len(sites)>=MIN_SITES and len(target_hectads)>=MIN_HECTADS
       and len(matched_names)>=MIN_MATCHED_HOST_TAXA and variation>=MIN_FOCAL_VARIATION)
    result={"schema":SCHEMA,
        "status":"BSBI_HISTORIC_HOST_SITE_SOURCE_OVERLAP_VERIFIED",
        "source":{
           "bsbi_original_MD5":MD5,
           "authorized_original_WCVP_host_links_sha256":HOST_SHA,
           "UKBMS_original_geography_crosswalk_sha256":"e38853b775ecadf39332a8412e2a4a2fa3c83512cf7218e4f1eda172dd83c9cd",
           "host_plant_records_source":"BSBI Plant Atlas 2020 publicly documented wild plants, not managed garden/crop plants",
           "period_pre1970":list(EARLY),
           "secondary_period_pre1987":list(PRE1987),
        },
        "source_ukbms_transect_count":len(sites),
        "source_10km_hectad_cells":len(target_hectads),
        "BSBI_3497_taxa_CSVs_verified":qc["original_bsbi_taxon_csv_files_scanned"],
        "BSBI_original_source_hectad_taxon_rows_in_selected_cells":qc["focal_hectad_taxon_source_rows"],
        "BSBI_botanical_recording_effort_positive_pre1970_minmax":[min(effort_pre1970.get(h,0) for h in target_hectads),max(effort_pre1970.get(h,0) for h in target_hectads)],
        "BSBI_botanical_recording_effort_positive_pre1987_minmax":[min(effort_pre1987.get(h,0) for h in target_hectads),max(effort_pre1987.get(h,0) for h in target_hectads)],
        "focal_original_UK_species":len(CROSSWALK),
        "original_unique_WCVP_accepted_plant_names":len(all_focal_hosts),
        "focal_exact_accepted_BSBI_binomial_host_plant_taxa":len(matched_names),
        "matched_host_species_with_any_original_pre1970_recorded_host_variation":variation,
        "taxon_match_source_quality":lookup_qc,
        "focal_butterfly_source_summary":per_species,
        "pre_model_information_gate":"PASS_WILD_HOST_HISTORICAL_SOURCE_CONTRAST_CANDIDATE" if eligible else "FAIL_WILD_HOST_HISTORICAL_SOURCE_CONTRAST_CANDIDATE",
        "DO_NOT_CLAIM":{
          "BSBI_0_is_confirmed_host_absence":False,
          "observed_plant_is_used_as_larval_host_at_this_site":False,
          "plants_human_maintained_or_cultivated":False,
          "host_effect_on_historic_butterfly_detection_tested":False,
          "host_mediated_climate_tracking_supported":False,
          "plant_taxon_x_hectad_and_UKBMS_site_data_exported":False,
        },
    }
    if return_internal:
        # Internal model-only source series, never included in public JSON.
        return result,{
            "site_to_hectad":hectad_by_site,
            "original_pre1970_site_botanical_effort":{site:effort_pre1970.get(hectad_by_site[site],0) for site in sites},
            "pre1970_observed_bsbi_plant_names_by_site":{site:frozenset(host_presence[hectad_by_site[site]]["early"]) for site in sites},
            "exact_source_host_plant_names_by_butterfly":{sp:frozenset(exact_host[sp]&matched_names) for sp in CROSSWALK},
            "all_source_plant_name_pool":frozenset(matched_names),
        }
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--bsbi-zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    res=site_historical_host_coverage(a.source_root,a.bsbi_zip)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(res,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":res["status"],"10km_cells":res["source_10km_hectad_cells"],
        "original_transects":res["source_ukbms_transect_count"],
        "BSBI_matched_host_taxa":res["focal_exact_accepted_BSBI_binomial_host_plant_taxa"],
        "source_27_species_with_pre1970_host_contrast":res["matched_host_species_with_any_original_pre1970_recorded_host_variation"],
        "source_plant_effort_pre1970_minmax":res["BSBI_botanical_recording_effort_positive_pre1970_minmax"],
        "decision":res["pre_model_information_gate"],
        "original_host_name_map_QC":res["taxon_match_source_quality"],
        },sort_keys=True),flush=True)

#!/usr/bin/env python3
"""Source-reconstruct global butterfly × GSHHG-island × year from original GBIF ZIP.

All 14,637,152 human-observation occurrences in DOI 10.15468/dl.umqmsf were
source-filtered: 225 GBIF butterfly taxonKeys, HUMAN_OBSERVATION, 2001–2024,
precise coordinate uncertainty 0–1,000m, valid coordinates, PRESENT. This script
streams their already completed 2026-10-29 GBIF archive *without* creating any
new download, and GIS-joins exact coordinates against all source GSHHG 2.3.7
high-resolution level-1 landmass polygons. The original private author Parquet
reported 14,425,536 indexed occurrences across 3,384 original sibling IDs, 224
source butterfly species, and 83,573 species×landmass×year rows; those figures
are benchmark diagnostics, NOT assumptions to force outcomes to agree.

The new global 3,696-island frame includes only GSHHG area >=20km² and <7m km²,
chosen independent of future butterfly index. The output 831,600 source
candidate butterfly×island cells includes GBIF recording nonindexes, NOT
ecological absences or known non-colonizations. External historical host
supply, shipping, or island source populations are not inferred.

No raw occurrence coordinates or potentially sensitive individual IDs are
committed; only species×island×year aggregates and source DOI are published.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import shapely
from shapely.strtree import STRtree

from build_gshhg_global_island_universe import source_shape

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
GBIF_DOI="10.15468/dl.umqmsf"
GBIF_KEY="0013264-261023074918325"
GBIF_ARCHIVE_SHA256="a1769481e91cb3002c25f69c6d8aa8e6d8d73a012c9495e8394032af02a9df18"
ZIP_MEMBER="0013264-261023074918325.csv"
SOURCE_RECORDS=14_637_152
BENCHMARK={"records":14_425_536,"rows":83_573,
           "species":224,"landmasses":3384}
PINS={
 "frame":("provenance/reviewer_defenses/results/global_island_risk_20261031/global_island_frame_area20.csv",
          "f0fabefd7b5e15d23c1ce7ad2127c8a03357f2a8"),
 "ports":("provenance/reviewer_defenses/results/global_island_risk_20261031/global_island_port_opportunity.csv",
          "42b13a4e2476d1ded979ff0c5425ae7e4e06864a"),
 "butterfly":("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv",
              "39fe84981cab2b0e47063399f0a5c0c9731e800d"),
 "metrics":("data/frozen/figure_sources/anthropogenic_species_metrics.csv",
            "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"),
 "migration":("provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
              "439dd1034ca16f96ff94222695d98dc39747c2e4"),
}
CHUNK=40000
def read(k):
    path,sha=PINS[k]
    raw=(ROOT/path).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if actual!=sha:
        raise ValueError("Source Git blob changed "+path+": "+actual)
    return list(csv.DictReader(io.StringIO(raw.decode("utf8"))))

def compute(zip_path,shapezip,work):
    archive_h=hashlib.sha256()
    with zip_path.open("rb") as f:
        while True:
            b=f.read(8*1024*1024)
            if not b:break
            archive_h.update(b)
    if archive_h.hexdigest()!=GBIF_ARCHIVE_SHA256:
        raise ValueError("Original GBIF archive bytes not source-identical to published DOI ZIP")
    frame=read("frame")
    ports=read("ports")
    butterflies=read("butterfly")
    metrics=read("metrics")
    migrant=read("migration")
    if len(frame)!=3696 or len(ports)!=3696 or len(butterflies)!=225:
        raise ValueError("Published global geometry / butterfly source coverage changed")
    original_ids={x["landmass_id"] for x in frame}
    if {x["landmass_id"] for x in ports}!=original_ids:
        raise ValueError("Global port frame source mismatches island identities")
    lookup={x["species"]:x for x in butterflies}
    species={x["taxon_key"]:x["species"] for x in butterflies}
    if len(species)!=225 or len(lookup)!=225:
        raise ValueError("Not 225 unique frozen butterfly species")
    port_by={x["landmass_id"]:x for x in ports}
    frame_by={x["landmass_id"]:x for x in frame}
    metrics_by={x["species"]:x for x in metrics}
    migrant_by={x["species"]:x["author_2025_GBIF_requested_taxon_key_proxy"]
        for x in migrant}
    land,source_gshhg=source_shape(shapezip,work)
    # 144,749 original WGS84 GSHHG level-1 polygon rows; the same authoritative
    # original sibling ID records as the previous private source calculation.
    geometry=np.asarray(land.geometry.to_numpy(),dtype=object)
    ids=np.asarray(land.landmass_id.tolist())
    if len(geometry)<140000:
        raise ValueError("Original GSHHG high-resolution geometry unexpectedly incomplete")
    # Prepared geometries accelerate 14m coordinate intersection tests.
    shapely.prepare(geometry)
    tree=STRtree(geometry)
    del land

    total=0;bad_rows=Counter();no_land=0;ambiguous_cross_landmass=0
    original_counter=Counter()
    global_counter=Counter()
    # Sparse count keyed by (speciesKey, originalGSHHGlandmassID, year).
    # Region-wide mainland records are excluded from candidate-island risk,
    # but still used to validate the published original author's landmass sum.
    buffer_species=[];buffer_year=[];buffer_lon=[];buffer_lat=[]
    index_year_outside_source=0
    def flush():
        nonlocal no_land,ambiguous_cross_landmass
        if not buffer_lat:return
        points=shapely.points(np.asarray(buffer_lon,dtype=float),
                              np.asarray(buffer_lat,dtype=float))
        hits=tree.query(points,predicate="intersects")
        xidx=hits[0].tolist();pidx=hits[1].tolist()
        if not xidx:
            no_land+=len(points)
            buffer_lat.clear();buffer_lon.clear();buffer_species.clear();buffer_year.clear()
            return
        mask=np.ones(len(points),dtype=bool)
        # Shapely point-to-polygon result may include nearby polygon boundaries.
        # Deduplicate a single coordinate to a *single original GSHHG sibling*.
        assigned={}
        ambiguous=set()
        for ip,ig in zip(xidx,pidx):
            island=str(ids[ig])
            if ip in assigned and assigned[ip]!=island:
                ambiguous.add(ip)
            else:assigned[ip]=island
        for ip in ambiguous:assigned.pop(ip,None)
        ambiguous_cross_landmass+=len(ambiguous)
        no_land+=len(points)-len(assigned)-len(ambiguous)
        for i,id in assigned.items():
            species_key=buffer_species[i]
            year=buffer_year[i]
            original_counter[(species_key,id,year)]+=1
            if id in original_ids:
                global_counter[(species_key,id,year)]+=1
        buffer_lat.clear();buffer_lon.clear();buffer_species.clear();buffer_year.clear()
    with zipfile.ZipFile(zip_path) as z:
        if ZIP_MEMBER not in z.namelist():
            raise ValueError("Original GBIF indexed 2026 ZIP member no longer matches DOI metadata")
        with z.open(ZIP_MEMBER) as bf:
            reader=csv.reader(io.TextIOWrapper(bf,encoding="utf8",errors="replace",newline=""),
                              delimiter="\t")
            headers=next(reader)
            col={n:i for i,n in enumerate(headers)}
            required=("speciesKey","decimalLatitude","decimalLongitude",
                      "coordinateUncertaintyInMeters","year","basisOfRecord",
                      "occurrenceStatus","gbifID")
            if any(k not in col for k in required):
                raise ValueError("Original ZIP schema changed "+repr(headers[:60]))
            ksp=col["speciesKey"];ky=col["year"];klat=col["decimalLatitude"];klon=col["decimalLongitude"]
            kprec=col["coordinateUncertaintyInMeters"];kbasis=col["basisOfRecord"]
            kstatus=col["occurrenceStatus"]
            for row in reader:
                total+=1
                if len(row)!=len(headers):
                    bad_rows["wrong_column_count"]+=1;continue
                if row[kbasis]!="HUMAN_OBSERVATION" or row[kstatus]!="PRESENT":
                    bad_rows["not_download_predicate_status"]+=1;continue
                sp=row[ksp]
                if sp not in species:
                    bad_rows["not_source_225_accepted_speciesKey"]+=1;continue
                try:
                    year=int(row[ky]);lat=float(row[klat]);lon=float(row[klon])
                    precision=float(row[kprec])
                except (TypeError,ValueError):
                    bad_rows["missing_numeric_year_coord_or_precision"]+=1;continue
                if not(2001<=year<=2024):
                    index_year_outside_source+=1;continue
                if not(-90<=lat<=90 and -180<=lon<=180 and 0<=precision<=1000
                       and math.isfinite(lat) and math.isfinite(lon)):
                    bad_rows["outside_geo_precision_limits"]+=1;continue
                buffer_species.append(sp);buffer_year.append(year)
                buffer_lon.append(lon);buffer_lat.append(lat)
                if len(buffer_lat)>=CHUNK:flush()
                if total%1_000_000==0:
                    print("ORIGINAL_GBIF_ROWS_PARSED",total,"original_GSHHG_aggregate_cells",len(original_counter),
                          "global_area20_cells",len(global_counter),flush=True)
            flush()
    if total!=SOURCE_RECORDS:
        raise ValueError("Published GBIF ZIP number of parsed occurrence rows changed: "+str(total))
    if any(v<0 for v in original_counter.values()):
        raise ValueError("Impossible source counts")
    original_total=sum(original_counter.values())
    benchmark={
        "matched_GSHHG_occurrences":original_total,
        "species_island_year_rows":len(original_counter),
        "unique_butterfly_species":len({k[0] for k in original_counter}),
        "unique_landmasses":len({k[1] for k in original_counter}),
        "benchmark":BENCHMARK,
        "matches_original_private_Author_Parquet_summary":{
            "total_records":original_total==BENCHMARK["records"],
            "rows":len(original_counter)==BENCHMARK["rows"],
            "species":len({k[0] for k in original_counter})==BENCHMARK["species"],
            "landmasses":len({k[1] for k in original_counter})==BENCHMARK["landmasses"]
        }
    }
    HERE.mkdir(exist_ok=True,parents=True)
    with gzip.open(HERE/"original_all_landmass_species_year_aggregates.csv.gz","wt",encoding="utf8",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["speciesKey","landmass_id","year","indexed_records"])
        for (key,landmass,year),n in sorted(original_counter.items()):
            writer.writerow([key,landmass,year,n])
    pair_year=defaultdict(Counter)
    for (sp,id,year),n in global_counter.items():
        pair_year[(sp,id)][year]+=n
    island_pre=Counter();island_pre_record=Counter()
    for (sp,id),y in pair_year.items():
        if any(year<=2010 for year in y):
            island_pre[id]+=1
            island_pre_record[id]+=sum(n for year,n in y.items() if year<=2010)
    island_early=Counter();island_late=Counter()
    species_early=Counter();species_late=Counter()
    global_eligible_early=0;global_eligible_late=0
    global_early=0;global_late=0
    summary_rows=[]
    with gzip.open(HERE/"global_225x3696_butterfly_island_risk.csv.gz","wt",encoding="utf8",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["species","speciesKey","landmass_id","pre2001_2010_indexed_records",
          "2011_2017_indexed_records","2018_2022_indexed_records",
          "pre2011_indexed_other_species_at_island","pre2011_island_candidate_species_index_count",
          "pre2011_total_candidate_butterfly_GBIF_records",
          "eligible_2011_2017","first_2011_2017","eligible_2018_2022","first_2018_2022",
          "GSHHG_equal_area_km2","reference_latitude","reference_longitude",
          "source_2019_WPI_port_count","WPI_2019_any_indexed_port",
          "GEB_introduced_share_of_contemporary","GEB_number_resolved_hosts",
          "independent_positive_migrant_roster_else_UNKNOWN",
          "GBIF_2026_reindex_not_true_butterfly_colonization"])
        for id in sorted(original_ids,key=int):
            area=frame_by[id]
            port=port_by[id]
            current_first=0;current_late=0;current_eligible=0;current_eligible_late=0
            current_pre=0
            for sp in sorted(species,key=int):
                name=species[sp]
                yy=pair_year.get((sp,id),{})
                npre=sum(n for y,n in yy.items() if y<=2010)
                nmiddle=sum(n for y,n in yy.items() if 2011<=y<=2017)
                nlate=sum(n for y,n in yy.items() if 2018<=y<=2022)
                eligible1=int(npre==0)
                first1=int(npre==0 and nmiddle>0)
                eligible2=int(npre==0 and nmiddle==0)
                first2=int(eligible2 and nlate>0)
                current_eligible+=eligible1;current_eligible_late+=eligible2
                current_first+=first1;current_late+=first2
                if npre>0:current_pre+=1
                species_early[sp]+=first1;species_late[sp]+=first2
                m=metrics_by.get(name,{})
                writer.writerow([
                  name,sp,id,npre,nmiddle,nlate,
                  island_pre[id]-int(npre>0),island_pre[id],island_pre_record[id],
                  eligible1,first1,eligible2,first2,
                  area["area_km2_equal_area"],area["reference_lat"],area["reference_lon"],
                  port["WPI_2019_assigned_port_count"],
                  port["WPI_2019_nearshore_assigned_any"],
                  m.get("introduced_share_of_contemporary",""),
                  m.get("resolved_host_species",""),
                  migrant_by.get(name,"NOT_LISTED_STATUS_UNKNOWN"),
                  True])
            global_eligible_early+=current_eligible;global_eligible_late+=current_eligible_late
            global_early+=current_first;global_late+=current_late
            summary_rows.append({
              "landmass_id":id,
              "area_km2":area["area_km2_equal_area"],
              "lat":area["reference_lat"],"lon":area["reference_lon"],
              "WPI_2019_ports":port["WPI_2019_assigned_port_count"],
              "WPI_2019_any_port":port["WPI_2019_nearshore_assigned_any"],
              "pre2011_species_indexed":current_pre,
              "eligible_2011_2017":current_eligible,
              "first_index_2011_2017":current_first,
              "eligible_2018_2022":current_eligible_late,
              "first_index_2018_2022":current_late
            })
    with (HERE/"global_island_butterfly_index_summary.csv").open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary_rows[0]))
        writer.writeheader();writer.writerows(summary_rows)
    species_rows=[{
      "species":species[sp],"speciesKey":sp,
      "first_index2011_2017_islands":species_early[sp],
      "first_index2018_2022_islands":species_late[sp],
      "GEB_introduced_resource_fraction":metrics_by.get(species[sp],{}).get("introduced_share_of_contemporary",""),
      "migration_review_positive_else_UNKNOWN":migrant_by.get(species[sp],"NOT_LISTED_STATUS_UNKNOWN")
    } for sp in sorted(species,key=int)]
    with (HERE/"global_species_butterfly_index_summary.csv").open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(species_rows[0]))
        writer.writeheader();writer.writerows(species_rows)
    global_port=Counter()
    for r in summary_rows:
        port=r["WPI_2019_any_port"]=="True"
        global_port[(port,"islands")]+=1
        global_port[(port,"eligible1")]+=r["eligible_2011_2017"]
        global_port[(port,"events1")]+=r["first_index_2011_2017"]
        global_port[(port,"eligible2")]+=r["eligible_2018_2022"]
        global_port[(port,"events2")]+=r["first_index_2018_2022"]
    result={
      "schema":"chocho_global_original_14m_GBIF_225_butterfly_x_3696_GSHHG_island_risk_v1",
      "status":"ORIGINAL_PUBLIC_GBIF_DERIVED_GLOBAL_ISLAND_RISK_NOT_CAUSAL_BUTTERFLY_COLONIZATION",
      "checked_at_utc":datetime.now(timezone.utc).isoformat(),
      "original_doi":"https://doi.org/"+GBIF_DOI,
      "original_GBIF_download_key":GBIF_KEY,
      "original_GBIF_zip_sha256":GBIF_ARCHIVE_SHA256,
      "all_downloaded_GBIF_rows":total,
      "original_GSHHG_source_sha256":source_gshhg,
      "original_GSHHG_source_polygon_rows":len(geometry),
      "GBIF_predicate_and_source_ZIP_previously_audited":True,
      "GBIF_rows_excluded_by_safety_checks":dict(bad_rows),
      "GBIF_years_outside_2001_2024":index_year_outside_source,
      "GBIF_points_unassigned_to_original_GSHHG_landmass":no_land,
      "GBIF_points_ambiguous_multiple_original_landmasses":ambiguous_cross_landmass,
      "original_private_panel_summary_reproduction":benchmark,
      "global_GSHHG_area20_islands":len(frame),
      "fixed_butterfly_species":len(species),
      "global_all_butterfly_species_x_island_risk_cells":225*len(frame),
      "global_2011_2017_eligible_indexing_cells":global_eligible_early,
      "global_2011_2017_first_indexed_cells":global_early,
      "global_2018_2022_eligible_indexing_cells":global_eligible_late,
      "global_2018_2022_first_indexed_cells":global_late,
      "observed_GBIF_2026_historical_index_events_by_modern_WPI_port":{
        "with_indexed_port":{k:global_port[(True,k)] for k in ("islands","eligible1","events1","eligible2","events2")},
        "without_indexed_port":{k:global_port[(False,k)] for k in ("islands","eligible1","events1","eligible2","events2")}
      },
      "source_git_blobs":{k:{"path":v[0],"sha":v[1]} for k,v in PINS.items()},
      "limits":[
        "Archive's current GBIF index is 2026-10-29, not an observational census of butterflies' true range.",
        "GBIF nonindexes are not biological absences and cannot identify true establishment or natural colonization.",
        "The archive contains human observations selected on 225 taxa, 2001–2024, and coordinate-uncertainty <=1km; not all global butterfly species.",
        "GSHHG 20km area cutoff may exclude smaller islands used by the original 12028 event chronology.",
        "Modern 2019 WPI port distribution is not historical cargo, path, or intensity; urbanization and recording bias confound comparisons.",
        "Separate survey effort, continent distance, historical climate and species-level island source-range accessibility are not yet adjusted.",
        "GEB globally introduced host portfolio is not actual island-specific historically introduced foodplants.",
        "Migration positive-only literature roster cannot infer nonmigration for unlisted butterflies.",
        "This is all 225 source candidates not post-selected 60 or 35 source butterflies."
      ]
    }
    (HERE/"global_original_14m_source_and_risk_receipt.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:result[k] for k in (
      "status","all_downloaded_GBIF_rows",
      "original_private_panel_summary_reproduction",
      "global_GSHHG_area20_islands","global_all_butterfly_species_x_island_risk_cells",
      "global_2011_2017_eligible_indexing_cells","global_2011_2017_first_indexed_cells",
      "global_2018_2022_eligible_indexing_cells","global_2018_2022_first_indexed_cells",
      "observed_GBIF_2026_historical_index_events_by_modern_WPI_port")},
      ensure_ascii=False,indent=2),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gbif-zip",type=Path,required=True)
    p.add_argument("--gshhg-zip",type=Path,required=True)
    p.add_argument("--shape-work",type=Path,required=True)
    args=p.parse_args()
    compute(args.gbif_zip,args.gshhg_zip,args.shape_work)

if __name__=="__main__":
    main()

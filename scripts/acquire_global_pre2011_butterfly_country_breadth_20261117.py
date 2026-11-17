#!/usr/bin/env python3
"""Source-pinned, country-faceted GBIF 2001–2010 records for 225 butterflies.

The independent exposure is geographical breadth of the butterfly's own
GBIF-indexed historical occurrence, NOT its true/native range or endemism.
All results must retain 225 complete request receipts: no query error or empty
facet is silently mapped to confirmed butterfly absence. Avoid downloading
individual location records; save non-sensitive aggregate country facets only.
Search API has a time-varying index; this is an audited 2026 snapshot, not an
immutable event-time download. HTTP 429 retries are bounded and concurrency
modest. The GBIF occurrence facets are a sampling index, not native status.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,math,random,time,urllib.parse,urllib.request,urllib.error
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[1]
SOURCE="data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
PIN="b6a3c98eb76864052d45d520a55b7d830816126c"
ENDPOINT="https://api.gbif.org/v1/occurrence/search"
QUERY={"year":"2001,2010","hasCoordinate":"true","hasGeospatialIssue":"false",
       "facet":"country","facetLimit":"260","limit":"0"}
AGENT="chocho-island-butterfly-geography/1.0 (academic research, 225 GBIF occurrence country facets; https://github.com/zuizui0223/chocho)"

def gitblob(path):
 data=Path(path).read_bytes()
 return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def facet_request(species,sk,tries=4):
 params=dict(QUERY,speciesKey=str(sk))
 url=ENDPOINT+"?"+urllib.parse.urlencode(sorted(params.items()))
 last_error=None
 for attempt in range(tries):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"})
   with urllib.request.urlopen(req,timeout=38) as res:
    if res.status!=200:raise RuntimeError("Non-200 "+str(res.status))
    data=res.read(5*1024*1024)
    response=json.loads(data)
   if not isinstance(response.get("count"),int) or not isinstance(response.get("facets"),list):
    raise ValueError("GBIF response missing count or facets")
   fs=[f for f in response["facets"] if str(f.get("field","")).upper()=="COUNTRY"]
   if len(fs)!=1 or not isinstance(fs[0].get("counts"),list):
    raise ValueError("GBIF country facet missing in source response")
   seen=set()
   countlist=[]
   for x in fs[0]["counts"]:
    code=str(x.get("name",""))
    count=int(x["count"])
    if len(code)!=2 or not code.isalpha() or code in seen or count<0:
     raise ValueError("Invalid or duplicated ISO-2 country facet "+repr(code))
    seen.add(code)
    countlist.append({"country":code,"indexed_occurrences":count})
   facets_total=sum(z["indexed_occurrences"] for z in countlist)
   if facets_total>response["count"]:
    raise ValueError("Country facet sum exceeds search response count")
   if response["count"]>0 and facets_total==0:
    raise ValueError("GBIF positive records with no country facet; no fabricated zero")
   receipt={
    "species":species,"speciesKey":int(sk),
    "status":"COMPLETE_GBIF_COUNTRY_FACETS",
    "gbif_uri":url,"GBIF_count_indexed_2001_2010":response["count"],
    "country_records_in_facets":facets_total,
    "country_facet_records_are_fraction_of_all":(facets_total/response["count"] if response["count"] else None),
    "n_countries_any_country_facet":len(countlist),
    "n_countries_5plus_records":sum(z["indexed_occurrences"]>=5 for z in countlist),
    "n_countries_10plus_records":sum(z["indexed_occurrences"]>=10 for z in countlist),
    "country_facets":countlist,
    "source_response_sha256":hashlib.sha256(data).hexdigest(),
    "acquired_at_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "query_attempts":attempt+1}
   return receipt
  except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError,KeyError,RuntimeError) as exc:
   last_error=str(exc)
   if isinstance(exc,urllib.error.HTTPError) and exc.code in (400,401,403,404):
    break
   if attempt<tries-1:time.sleep((2**attempt)*1.8+random.Random(int(sk)+attempt).random())
 return {"species":species,"speciesKey":int(sk),"status":"GBIF_QUERY_FAILED_NOT_BIOLOGICAL_ABSENCE",
         "error":last_error,"gbif_uri":url,"query_attempts":tries}

def main():
 pa=argparse.ArgumentParser()
 pa.add_argument("--out-dir",type=Path,required=True)
 pa.add_argument("--pilot",action="store_true",help="12 fixed-source butterflies to validate API scope")
 pa.add_argument("--workers",type=int,default=4)
 a=pa.parse_args()
 if a.workers>5 or a.workers<1:raise ValueError("Concurrent GBIF API cap 1-5")
 original=ROOT/SOURCE
 if gitblob(original)!=PIN:
  raise ValueError("Original 225×51 GBIF source identity changed")
 with original.open(encoding="utf-8",newline="") as f:
  rows=list(csv.DictReader(f))
 pairs=defaultdict(set)
 for row in rows:
  pairs[row["species"]].add(row["speciesKey"])
 if len(rows)!=225*51 or len(pairs)!=225 or any(len(v)!=1 for v in pairs.values()):
  raise ValueError("Not original 225 butterflies × 51 landmasses")
 allspecies=sorted((name,int(next(iter(keys)))) for name,keys in pairs.items())
 if a.pilot:
  allspecies=[allspecies[i] for i in (0,13,25,50,70,90,110,130,155,180,205,224)]
 expected=len(allspecies)
 outputs=[]
 with ThreadPoolExecutor(max_workers=a.workers) as executor:
  futures={executor.submit(facet_request,sp,sk):sp for sp,sk in allspecies}
  for future in as_completed(futures):
   z=future.result();outputs.append(z)
   print(f"GBIF_COUNTRY_SOURCE {z['species']} {z['status']} "
         f"{z.get('n_countries_any_country_facet','NA')} countries",flush=True)
 outputs.sort(key=lambda r:r["species"])
 nfail=sum(z["status"]!="COMPLETE_GBIF_COUNTRY_FACETS" for z in outputs)
 success=[z for z in outputs if z["status"]=="COMPLETE_GBIF_COUNTRY_FACETS"]
 result={
  "schema":"chocho_global_225_butterflies_GBIF_2001_2010_georeferenced_country_facets_v1",
  "status":"ALL_REQUESTS_COMPLETE" if nfail==0 else "PARTIAL_GBIF_REQUEST_FAILURES_RETAINED_NOT_AS_ABSENCE",
  "pilot":a.pilot,
  "requested":expected,"completed":len(success),"failed":nfail,
  "original_species_by_landmass_index_blob":PIN,
  "original_butterfly_panel":225,"original_landmass_panel":51,
  "GBIF_occurrence_endpoint":ENDPOINT,"fixed_parameters":QUERY,
  "obtained_at_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
  "data_use":{
   "source":"GBIF occurrence/search facets, GBIF Backbone numeric speciesKey",
   "relative_period":"2001–2010 only, BEFORE butterfly 2011–2022 outcome records in source 51 landmasses",
   "interpretation":"Observed country indexing breadth, not native butterfly range, historical endemicity, abundance or established populations",
   "reproducibility":"Timestamped independent API snapshot, neither a GBIF record-ID download nor source-provider DOI",
   "rate_limit_compliance":"At most 4 concurrent country facet calls, bounded HTTP429 backoff; no screen scraping",
   "rights":"No coordinates or individual butterflies recorded in this archive; counts aggregated by GBIF country facet"},
  "countries_per_butterfly":outputs}
 a.out_dir.mkdir(parents=True,exist_ok=True)
 fname="PILOT_" if a.pilot else "SOURCE_"
 (a.out_dir/(fname+"GBIF2001_2010_225COUNTRY_FACETS_20261117.json")).write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
 with (a.out_dir/(fname+"country_breadth_by_butterfly.csv")).open("w",newline="") as f:
  writer=csv.DictWriter(f,fieldnames=["species","speciesKey","status","GBIF_count_indexed_2001_2010",
                                      "country_records_in_facets","n_countries_any_country_facet",
                                      "n_countries_5plus_records","n_countries_10plus_records"])
  writer.writeheader()
  for z in outputs:
   writer.writerow({k:z.get(k,"") for k in writer.fieldnames})
 print("GBIF_COUNTRY_FACET_SOURCE_COMPLETED="+json.dumps({"requested":expected,"success":len(success),"failed":nfail,
       "no_country_species":sum(z["n_countries_any_country_facet"]==0 for z in success),
       "total_country_count":sum(z["n_countries_any_country_facet"] for z in success)},sort_keys=True),flush=True)
 if nfail:
  raise RuntimeError(f"{nfail} GBIF facets failed; no missing values have been recoded as butterfly range zeros")

if __name__=="__main__":main()

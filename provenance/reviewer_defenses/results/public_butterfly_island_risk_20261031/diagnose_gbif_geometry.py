#!/usr/bin/env python3
"""Isolate GBIF 400 error without changing analytic input cohort or WKT shapes.

Compare same calendar year/facets under no geometry, simple valid rectangle,
source GSHHG original area-preserving island polygon, and one vs six families.
Record HTTP error body and URL character count. Diagnostics are not butterfly
outcomes and do not show biological absence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from datetime import datetime,timezone
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import build_public_island_risk as src

def test_case(name,params):
    url=src.GBIF+"/occurrence/search?"+urlencode(params,doseq=True)
    result={"name":name,"url_chars":len(url),
      "year":params.get("year"),"families":params.get("family"),
      "has_geometry":"geometry" in params,
      "geometry_sha256":hashlib.sha256(str(params.get("geometry","")).encode()).hexdigest()
                    if params.get("geometry") else None,
      "geometry_chars":len(params.get("geometry","")),
    }
    try:
        with urlopen(Request(url,headers={"User-Agent":src.UA,"Accept":"application/json"}),timeout=65) as f:
            payload=f.read(1_200_000).decode("utf8")
            data=json.loads(payload)
            result.update({"http":f.status,"count":data.get("count"),
                           "facets":[x.get("field") for x in data.get("facets",[])],
                           "error":data.get("error")})
    except HTTPError as err:
        result.update({"http":err.code,"error_type":"HTTPError",
          "error_body":err.read(3000).decode("utf8",errors="replace")})
    except (URLError,TimeoutError,ValueError) as err:
        result.update({"error_type":type(err).__name__,"error_body":str(err)[:700]})
    print("CASE",json.dumps(result,ensure_ascii=False),flush=True)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gshhg-zip",type=Path,required=True)
    ap.add_argument("--work",type=Path,required=True)
    args=ap.parse_args()
    source_islands={"169"}
    polygons,audit,digest=src.load_polygons(args.gshhg_zip,args.work,source_islands,1800)
    wkt=polygons["169"]
    families=["Nymphalidae","Pieridae","Lycaenidae","Papilionidae","Hesperiidae","Riodinidae"]
    base={"year":"2012","hasCoordinate":"true","hasGeospatialIssue":"false",
          "occurrenceStatus":"PRESENT","facet":"speciesKey","facetLimit":10,
          "limit":0}
    rectangle="POLYGON((9.9 55.0,10.3 55.0,10.3 55.3,9.9 55.3,9.9 55.0))"
    cases=[
      ("country_no_family",{**base,"country":"DK"}),
      ("country_one_family",{**base,"country":"DK","family":"Nymphalidae"}),
      ("country_six_families",{**base,"country":"DK","family":families}),
      ("country_nonsense_family",{**base,"country":"DK","family":"Nonsenseidae"}),
      ("country_known_aglais_species_taxonKey",{**base,"country":"DK","taxonKey":"1894840"}),
      ("rectangle_one_family",{**base,"geometry":rectangle,"family":"Nymphalidae"}),
      ("rectangle_six_families",{**base,"geometry":rectangle,"family":families}),
      ("original_polygon_one_family",{**base,"geometry":wkt,"family":"Nymphalidae"}),
      ("original_polygon_six_families",{**base,"geometry":wkt,"family":families}),
    ]
    out=[test_case(name,params) for name,params in cases]
    receipt={"schema":"chocho_GBIF_2026_butterfly_polygon_400_diagnosis_v1",
      "checked_at_utc":datetime.now(timezone.utc).isoformat(),
      "source_gshhg_sha256":digest,"geometry_source_id":"169",
      "geometry_audit":audit,
      "diagnostic_only_no_butterfly_absence_inference":True,
      "results":out}
    path=HERE/"gbif_polygon_api_diagnostic.json"
    path.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    if not any(x.get("http")==200 for x in out):
        raise ValueError("No successful GBIF occurrence reference query; source unavailable")

if __name__=="__main__":main()

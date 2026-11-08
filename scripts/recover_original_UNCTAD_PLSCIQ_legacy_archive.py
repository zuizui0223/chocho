#!/usr/bin/env python3
"""Recover 3 historic UNCTAD PLSCIQ/DBnomics source pages and source-only WPI194 join.

Original 2026-10-30 source vintage/status/hashes recorded in chocho historical audit.
Do not label newer/reindexed API bytes as an exact old file recovery.
PLSCIQ is a port-level liner network index, NOT real ship voyages.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,time,unicodedata
import urllib.request,urllib.error
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter,defaultdict
try:
    import pycountry
except ImportError:
    pycountry=None

ROOT=Path(__file__).resolve().parents[1]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
API="https://api.db.nomics.world/v22/series/UNCTAD/PLSCIQ"
HISTORIC_DIR_HASH="5cfa0ae00d0c42645dae30be6a5fad068dc35c5f"
HISTORIC_SHA={
    0:"29a343bb6b1113b1699a9ca1821952cce1cb8f150831921da68e248e78fb3e1b",
    500:"c82e2d2119b0f4b03e4afb18c04f9c4442f19a2bfce96efbc8d30a87aed42933",
    1000:"7eb1eaf4349008b784d9e6d2bd834cadf29c51fa493c2cdd2f83e2d767978255",
}
Q_START="2006-Q1"
Q_STOP="2021-Q3"
EARLY_YEARS=set(range(2006,2011))

def git_blob_sha(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def clean(s):
    return re.sub("[^A-Z0-9]","",unicodedata.normalize("NFKD",str(s)).encode("ascii","ignore").decode().upper())

def country_iso2(s):
    if pycountry is None:
        raise RuntimeError("pycountry required to avoid a guessed port country join")
    try:
        return pycountry.countries.lookup(str(s)).alpha_2
    except LookupError:
        return ""

def original_wpi():
    if git_blob_sha(WPI)!=WPI_BLOB:
        raise ValueError("Exact original WPI port coordinate table Git blob changed")
    with WPI.open("rt",encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=194 or len({r["port_id"] for r in rows})!=194 or len({r["island_id"] for r in rows})!=51:
        raise ValueError("Original 194 biological WPI port / 51 landmass source differs")
    if sum(r["included"].lower()=="true" for r in rows)!=182:
        raise ValueError("Original source marine qualified port count differs")
    return rows

def fetch_page(offset,path):
    url=API+f"?limit=500&offset={offset}&observations=1"
    if path.is_file():
        raw=path.read_bytes()
        return raw,{"status":"LOCAL_CACHE","url":url,"bytes":len(raw)}
    problem=None
    for trial in range(4):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"chocho-UNCTAD-PLSCIQ-original-port-source-recovery/1.0",
                "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=90) as response:
                code=response.status
                raw=response.read(35_000_000)
            if code!=200 or len(raw)>=35_000_000:
                raise ValueError("Public DBnomics response incomplete or too big")
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(raw)
            return raw,{"status":code,"url":url,"bytes":len(raw)}
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ConnectionError) as exc:
            problem=f"{type(exc).__name__}: {str(exc)[:220]}"
            if trial<3:time.sleep(4*(trial+1))
    raise RuntimeError(f"Public UNCTAD/DBnomics source unavailable at offset {offset}: {problem}")

def parse_value(v):
    if v is None or str(v).strip() in ("","NA","NaN","nan","null"):
        return None
    try:
        n=float(v)
        if not math.isfinite(n):
            return None
        return n
    except (ValueError,TypeError):
        return None

def save_csv(path,rows,fields):
    with Path(path).open("wt",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def run(out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("data","results","provenance")):
        raise ValueError("Do not overwrite original butterfly source or previous audit")
    out.mkdir(parents=True,exist_ok=True)
    original=original_wpi()
    sources=[]
    oldhash_matches=[]
    cached=out/"original_2021_source_pages"
    all_docs=[]
    meta=None
    for offset in (0,500,1000):
        raw,fetch=fetch_page(offset,cached/f"offset_{offset}.json")
        sha=hashlib.sha256(raw).hexdigest()
        same=(sha==HISTORIC_SHA[offset])
        doc=json.loads(raw)
        if "dataset" not in doc or "series" not in doc or not isinstance(doc["series"].get("docs"),list):
            raise ValueError("Source UNCTAD dataset/series shape differs from archived direct API")
        batch=doc["series"]["docs"]
        dataset=doc["dataset"]
        if meta is None:
            meta=dataset
        elif dataset.get("dir_hash")!=meta.get("dir_hash"):
            raise ValueError("UNCTAD archive pages have inconsistent provider source vintage")
        if int(doc["series"].get("offset",-1))!=offset:
            raise ValueError("Original page offset not equal to requested page")
        if len(batch)!=(300 if offset==1000 else 500):
            raise ValueError("Original 1300-series source page count changed")
        all_docs.extend(batch)
        oldhash_matches.append(same)
        sources.append({"offset":offset,"count":len(batch),
            "sha256":sha,"matches_original_2026_10_30_sha256":same,
            "old_sha256":HISTORIC_SHA[offset],"fetch":fetch})
        print("UNCTAD_PLSCIQ_SOURCE_PAGE="+json.dumps({
            "offset":offset,"count":len(batch),"sha256":sha,
            "matches_old_cache":same,"provider_dir_hash":dataset.get("dir_hash")},sort_keys=True),flush=True)
    if len(all_docs)!=1300 or len({d["series_code"] for d in all_docs})!=1300:
        raise ValueError("UNCTAD original 1300 port series incomplete")
    labels=meta["dimensions_values_labels"]["port"]
    audit_version={
        "source_dir_hash":meta.get("dir_hash"),
        "is_exact_2026_10_30_original_recovered":all(oldhash_matches) and meta.get("dir_hash")==HISTORIC_DIR_HASH,
        "provider_updated_at":meta.get("updated_at"),
        "provider_indexed_at":meta.get("indexed_at"),
        "provider_label_measure":meta.get("dimensions_values_labels",{}).get("measure",{})
    }
    series=[]
    observations=0
    nonmissing=0
    first_period=None;last_period=None
    for item in all_docs:
        code=item["series_code"]
        pcode=item["dimensions"]["port"]
        title=labels[pcode]
        if ", " not in title:
            raise ValueError("Provider country/port label lacks comma delimiter")
        country,port=title.split(", ",1)
        periods=item["period"]
        values=item["value"]
        if len(periods)!=len(values):
            raise ValueError("Provider series period count does not match values")
        q={}
        for p,x in zip(periods,values):
            observations+=1
            v=parse_value(x)
            q[p]=v
            if v is not None:nonmissing+=1
            first_period=p if first_period is None or p<first_period else first_period
            last_period=p if last_period is None or p>last_period else last_period
        early=[v for p,v in q.items() if p[:4].isdigit() and int(p[:4]) in EARLY_YEARS and v is not None]
        series.append({"series_code":code,"provider_port_code":pcode,"source_country_label":country,
            "source_port_label":port,"ISO2_from_pycountry_if_resolved":country_iso2(country),
            "port_key":clean(port),"quarterly":q,
            "nonmissing_2006_2010":len(early),
            "mean_2006_2010":sum(early)/len(early) if early else None,
            "nonmissing_all_periods":sum(v is not None for v in q.values())})
    if observations!=81900 or first_period!=Q_START or last_period!=Q_STOP:
        raise ValueError("Legacy 2021-edition series support differs from historical audit")
    if audit_version["is_exact_2026_10_30_original_recovered"] and nonmissing!=57312:
        raise ValueError("Exact old cache SHA matched but numerical nonmissing count changed")
    # Anonymous source ID/country and exact unambiguous normalized port-name matches only.
    byprovider=defaultdict(list)
    for s in series:
        if s["ISO2_from_pycountry_if_resolved"]:
            byprovider[(s["ISO2_from_pycountry_if_resolved"],s["port_key"])].append(s)
    byoriginal=Counter((p["country_code"],clean(p["port_name"])) for p in original)
    coverage=[]
    for p in original:
        pair=(p["country_code"],clean(p["port_name"]))
        hits=byprovider.get(pair,[])
        unique=(byoriginal[pair]==1 and len(hits)==1)
        s=hits[0] if unique else None
        coverage.append({
            "original_WPI_id":p["port_id"],
            "original_island_id":p["island_id"],
            "original_WPI_country":p["country_code"],
            "original_WPI_name":p["port_name"],
            "status":("UNIQUE_COUNTRY_AND_CANONICAL_NAME_CANDIDATE_NOT_GEO_VERIFIED" if unique
                      else "AMBIGUOUS_NAME_COUNTRY" if hits else "NO_PROVIDER_NAME_COUNTRY_MATCH"),
            "provider_series_code":s["series_code"] if s else "",
            "provider_port_code":s["provider_port_code"] if s else "",
            "provider_port_name":s["source_port_label"] if s else "",
            "available_2006_2010_quarters":s["nonmissing_2006_2010"] if s else "",
            "mean_2006_2010_PLCSIQ_old_index_base":s["mean_2006_2010"] if s and s["mean_2006_2010"] is not None else "",
            "available_all_quarters":s["nonmissing_all_periods"] if s else "",
            "independent_port_geolocation_verified":False,
            "provider_series_is_not_actual_od_voyage":True,
        })
    save_csv(out/"UNCTAD_PLSCIQ_originalWPI194_name_country_candidates.csv",coverage,list(coverage[0]))
    matched=[r for r in coverage if r["provider_series_code"]]
    receipt={
        "schema":"chocho_UNCTAD_PLSCIQ_2006_2021_legacy_1300_source_recovery_v1",
        "status":"ALL_1300_PORT_SERIES_AND_81900_QUARTER_VALUES_EXAMINED",
        "source_provider":"UNCTAD",
        "indexed_edition":"DBnomics PLSCIQ archive",
        "source_url":API,
        "original_2026_10_30_archive_sha256_pages":HISTORIC_SHA,
        "page_receipts":sources,
        "source_version":audit_version,
        "original_archive_expected_dir_hash":HISTORIC_DIR_HASH,
        "series":len(series),"total_port_quarter_rows":observations,
        "nonmissing_port_quarter_rows":nonmissing,
        "first_period":first_period,"last_period":last_period,
        "original_WPI_source_ports":len(original),
        "source_unique_name_country_WPI_candidates":len(matched),
        "source_early_16of20_matched_ports":sum(int(r["available_2006_2010_quarters"])>=16 for r in matched),
        "source_early_16of20_original_islands":len({r["original_island_id"] for r in matched
            if int(r["available_2006_2010_quarters"])>=16}),
        "all_1300_series_raw_values_not_redistributed_publicly":True,
        "no_geographic_port_identity_verified":True,
        "port_connectivity_index_not_2020_voyages_or_butterfly_transport":True,
        "source_retrieval_utc":datetime.now(timezone.utc).isoformat()
    }
    (out/"UNCTAD_PLSCIQ_source_recovery_receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("UNCTAD_PLSCIQ_RECOVERY_RESULT_JSON="+json.dumps({
        "series":len(series),"observations":observations,"nonmissing":nonmissing,
        "old_pages_identical":audit_version["is_exact_2026_10_30_original_recovered"],
        "WPI194_name_country_candidates":len(matched),
        "early_at_least_16_quarters":receipt["source_early_16of20_matched_ports"],
        "early_original_islands":receipt["source_early_16of20_original_islands"]},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,required=True)
    run(p.parse_args().out_dir)

#!/usr/bin/env python3
"""Source-pinned Illinois/Michigan/Colorado monitoring test of host-niche overlap.

Outcome is adult survey counts per monitored hour, NOT larvae or demographic
population size. An observational lagged association does not prove competition.

No raw source site, event, date, count row, or sensitive coordinates are persisted.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import urllib.request
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
URL = "https://ndownloader.figshare.com/files/50896473"
SOURCE_MD5 = "3161b1c8802a0d67771b2044e18da928"
SOURCE_LENGTH = 55066096
FOCAL_SHA = "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"
PROGRAMS = (
    "Illinois Butterfly Monitoring Network",
    "Michigan Butterfly Network",
    "Colorado Butterfly Monitoring Network",
)
MIN_EVENTS = 5
MIN_SEASON_BINS = 3
MIN_SITE_YEARS = 3
MIN_SHARED_PAIRS = 20
MIN_TRANSITIONS = 300
MIN_SITES = 50


def load_source(path):
    md = hashlib.md5(usedforsecurity=False)
    length = 0
    request = urllib.request.Request(
        URL, headers={"User-Agent": "chocho-niche-dynamics-research/0.1"}
    )
    with urllib.request.urlopen(request, timeout=90) as r, path.open("wb") as out:
        if r.status != 200:
            raise ValueError("Original Figshare source unavailable")
        while True:
            b = r.read(1024 * 1024)
            if not b:
                break
            length += len(b)
            if length > 60_000_000:
                raise ValueError("Original Figshare source exceeds pinned limit")
            md.update(b)
            out.write(b)
    if md.hexdigest() != SOURCE_MD5 or length != SOURCE_LENGTH:
        raise ValueError("Original Figshare archive checksum or size drift")


def frozen_focal():
    path = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    raw = path.read_bytes()
    blob = hashlib.sha1(
        b"blob " + str(len(raw)).encode() + b"\x00" + raw
    ).hexdigest()
    if blob != FOCAL_SHA:
        raise ValueError("Original butterfly input changed")
    with path.open(newline="", encoding="utf-8") as f:
        names = {row["species"] for row in csv.DictReader(f)}
    if len(names) != 239:
        raise ValueError("Original 239 species changed")
    return names


def host_sets(path, focal):
    hosts = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            sp = row["insect_species"].strip()
            plant = row["accepted_plant_name_id"].strip()
            if sp in focal and plant:
                hosts[sp].add(plant)
    missing = focal - hosts.keys()
    if len(focal) != 239 or len(missing) != 0:
        raise ValueError(f"Accepted HOSTS butterfly join lost taxa: {len(missing)}")
    return {k: frozenset(v) for k, v in hosts.items()}


def read_events(src, focal):
    # Each survey (not each species row) owns one effort measure.
    events = {}
    by_program = Counter()
    with zipfile.ZipFile(src) as zip:
        info = [i for i in zip.infolist()
                if not i.filename.startswith("__MACOSX/") and
                i.filename.endswith("/4_res/complete-share.csv")]
        if len(info) != 1:
            raise ValueError("Original count CSV member missing")
        with zip.open(info[0]) as binary:
            reader = csv.DictReader(
                io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
            )
            expected = ("program", "site", "event.id", "year", "doy",
                        "Reported species name", "count", "duration")
            if any(x not in (reader.fieldnames or []) for x in expected):
                raise ValueError("Unexpected original count schema")
            for r in reader:
                p = r["program"].strip()
                if p not in PROGRAMS:
                    continue
                by_program[p] += 1
                site = r["site"].strip()
                event = r["event.id"].strip()
                species = r["Reported species name"].strip()
                year = int(float(r["year"]))
                doy = int(float(r["doy"]))
                if not site or not event or year < 2000 or year > 2020 or not 1 <= doy <= 366:
                    raise ValueError("Invalid original transect date/site")
                count = float(r["count"])
                if not math.isfinite(count) or count < 0:
                    raise ValueError("Invalid butterfly count")
                t = (p, event)
                value = r.get("duration", "").strip()
                try:
                    dur = float(value)
                except ValueError:
                    dur = float("nan")
                if t not in events:
                    events[t] = {
                        "site": site, "year": year, "doy": doy,
                        "duration": dur, "bad_duration": not (math.isfinite(dur) and dur > 0),
                        "counts": defaultdict(float),
                    }
                v = events[t]
                if v["site"] != site or v["year"] != year or v["doy"] != doy:
                    raise ValueError("Ambiguous source survey identity across site/year/day")
                if not math.isfinite(dur) or dur <= 0 or (
                    math.isfinite(v["duration"]) and abs(v["duration"] - dur) > 0.01
                ):
                    v["bad_duration"] = True
                if species in focal:
                    v["counts"][species] += count
    return events, by_program


def calc_panel(program, events, focal, hostgraph):
    site_years = defaultdict(list)
    dropped_bad_duration = 0
    for (p, e), ev in events.items():
        if p != program:
            continue
        if ev["bad_duration"]:
            dropped_bad_duration += 1
            continue
        site_years[(ev["site"], ev["year"])].append(ev)
    eligible_years = {}
    for key, values in site_years.items():
        if len(values) < MIN_EVENTS or len({(e["doy"] - 1)//31 for e in values}) < MIN_SEASON_BINS:
            continue
        dur = sum(e["duration"] for e in values)
        if not math.isfinite(dur) or dur <= 0:
            continue
        counts = defaultdict(float)
        for ev in values:
            for s, count in ev["counts"].items():
                counts[s] += count
        eligible_years[key] = {"effort_hr": dur/60, "counts": counts,
                               "visits": len(values)}
    years_by_site = defaultdict(set)
    for site, year in eligible_years:
        years_by_site[site].add(year)
    valid_sites = {
        st for st, years in years_by_site.items()
        if any(y+1 in years and y+2 in years for y in years)
    }
    eligible_years = {k:v for k,v in eligible_years.items() if k[0] in valid_sites}
    years_by_site = defaultdict(set)
    taxa_by_site = defaultdict(set)
    for (site, year), data in eligible_years.items():
        years_by_site[site].add(year)
        for name, value in data["counts"].items():
            if value > 0:
                taxa_by_site[site].add(name)
    programme_species = sorted(set().union(*taxa_by_site.values()) if taxa_by_site else set())
    if not programme_species:
        return {"status":"NO_STANDARD_TRANSECT_RECORDS","programme":program}
    all_pairs = [(a,b) for n,a in enumerate(programme_species)
                 for b in programme_species[n+1:]]
    shared_pairs = [
        (a,b) for a,b in all_pairs if hostgraph[a] & hostgraph[b]
    ]
    n_imputed = 0
    n_valid_records = 0
    positive = 0
    annual = {}
    for (site, year), d in eligible_years.items():
        for sp in taxa_by_site[site]:
            value = float(d["counts"].get(sp,0))
            if sp not in d["counts"]:
                n_imputed += 1
            else:
                n_valid_records += 1
            if value > 0:
                positive += 1
            annual[(site,year,sp)] = math.log1p(value/d["effort_hr"])

    observations = []
    for site, years in years_by_site.items():
        taxa = sorted(taxa_by_site[site])
        for year in sorted(years):
            if year+1 not in years:
                continue
            log_rate = {
                sp:annual[(site,year,sp)] for sp in taxa
            }
            mean_log_rate=sum(log_rate.values())/len(log_rate)
            annual_hrs=eligible_years[(site,year)]["effort_hr"]
            annual_hrs_next=eligible_years[(site,year+1)]["effort_hr"]
            for sp in taxa:
                same = [log_rate[t] for t in taxa
                        if t != sp and hostgraph[t]&hostgraph[sp]]
                other = [log_rate[t] for t in taxa
                         if t != sp and not(hostgraph[t]&hostgraph[sp])]
                if not same or not other:
                    continue
                outcome=annual[(site,year+1,sp)]-log_rate[sp]
                observations.append({
                    "group":(site,sp),
                    "site":site,
                    "year":year,
                    "x":[log_rate[sp],sum(same)/len(same),sum(other)/len(other),
                         math.log(annual_hrs_next/annual_hrs),mean_log_rate],
                    "y":outcome,
                })
    result={
        "status":"PANEL_PREPARED_NO_ECOLOGICAL_INFERENCE",
        "programme":program,
        "source_years":sorted({y for _,y in eligible_years}),
        "sites_eligible":len(valid_sites),
        "site_years_eligible":len(eligible_years),
        "species_matched":len(programme_species),
        "source_host_sharing_pairs":len(shared_pairs),
        "source_nonsharing_pairs":len(all_pairs)-len(shared_pairs),
        "year_to_next_year_pairs":sum(y+1 in ys for ys in years_by_site.values() for y in ys),
        "site_species_transitions_with_both_comparator_groups":len(observations),
        "site_species_groups_eligible":len({x["group"] for x in observations}),
        "site_year_species_missing_imputed_zeros":n_imputed,
        "site_year_species_observed_entries":n_valid_records,
        "site_year_species_positive":positive,
        "event_entries_dropped_inconsistent_duration":dropped_bad_duration,
        "local_record_site_ids_and_counts_not_persisted":True,
        "source_implied_zero_caution":"Missing species in otherwise valid transect survey assumed 0 only for ever-observed-at-site species",
    }
    gate=(len(shared_pairs)>=MIN_SHARED_PAIRS
          and len(observations)>=MIN_TRANSITIONS
          and len({o["site"] for o in observations})>=MIN_SITES)
    result["analysis_gate_pass"]=gate
    return result,observations


def group_within_matrix(obs):
    years=sorted({r["year"] for r in obs})
    other_years=years[1:]
    groups=defaultdict(list)
    for i,r in enumerate(obs):
        groups[r["group"]].append(i)
    columns=5+len(other_years)
    X=np.zeros((len(obs),columns),dtype=float)
    y=np.array([r["y"] for r in obs])
    for i,r in enumerate(obs):
        X[i,:5]=r["x"]
        if r["year"] in other_years:
            X[i,5+other_years.index(r["year"])]=1.0
    # Frisch–Waugh species x site fixed-effects transformation.
    for indices in groups.values():
        X[indices]-=X[indices].mean(axis=0)
        y[indices]-=y[indices].mean()
    return X,y,groups


def regression(obs):
    X,y,groups=group_within_matrix(obs)
    n,k=X.shape
    if n<=k+100 or np.linalg.matrix_rank(X) != k:
        return {"status":"SINGULAR_FIXED_EFFECT_PANEL"}
    coef=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@coef
    source_sites=sorted({r["site"] for r in obs})
    cl=defaultdict(list)
    for i,r in enumerate(obs):
        cl[r["site"]].append(i)
    bread=np.linalg.pinv(X.T@X)
    meat=np.zeros((k,k))
    for index in cl.values():
        vec=X[index].T@resid[index]
        meat+=np.outer(vec,vec)
    correction=len(cl)/(len(cl)-1)*(n-1)/(n-k)
    V=correction*bread@meat@bread
    se=np.sqrt(np.maximum(0,np.diag(V)))
    delta=coef[1]-coef[2]
    diff_se=math.sqrt(max(0,V[1,1]+V[2,2]-2*V[1,2]))
    return {
        "status":"OBSERVATIONAL_ADULT_COUNT_ASSOCIATION_ONLY",
        "sample_transitions":n,
        "site_clusters":len(cl),
        "species_site_effects":len(groups),
        "n_covariates_after_FEs":k,
        "beta_own_lag":float(coef[0]),"SE_own_sitecluster":float(se[0]),
        "beta_shared_host_peer_lag":float(coef[1]),"SE_shared_sitecluster":float(se[1]),
        "beta_nonshared_peer_lag":float(coef[2]),"SE_nonshared_sitecluster":float(se[2]),
        "contrast_shared_minus_nonshared":float(delta),
        "SE_contrast_sitecluster":float(diff_se),
        "contrast_normal_95_interval":[float(delta-1.96*diff_se),float(delta+1.96*diff_se)],
        "within_site_year_beta_not_causal_competition":True,
        "needs_host_degree_null_and_independent_habitat_weather_controls_before_a_new_ecological_claim":True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--accepted-host-csv",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    focal=frozen_focal()
    h=host_sets(args.accepted_host_csv,focal)
    with tempfile.TemporaryDirectory() as d:
        src=Path(d)/"original.zip"
        load_source(src)
        events,counts=read_events(src,focal)
    programmes=[]
    for prog in PROGRAMS:
        resp=calc_panel(prog,events,focal,h)
        if isinstance(resp,dict):
            programmes.append(resp)
            continue
        report,obs=resp
        if report["analysis_gate_pass"]:
            report["model"]=regression(obs)
        else:
            report["status"]="INELIGIBLE_FOR_PREDECLARED_HOST_NICHE_LAG_TEST"
        programmes.append(report)
    output={
        "schema":"chocho_us_sbus_illinois_host_niche_lagged_count_pilot_20261025_v01",
        "status":"EXPLORATORY_OBSERVATIONAL_COUNTS_NOT_DEMOGRAPHIC_CAUSAL_PROOF",
        "original_source_doi":"10.6084/m9.figshare.27934602.v1",
        "source_archive_md5":SOURCE_MD5,
        "host_database_pin":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
        "accepted_taxonomy_pin":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
        "n_focal_in_source":len(focal),
        "source_programme_rows":{p:counts[p] for p in PROGRAMS},
        "programmes":programmes,
        "interpretation_restrictions":[
            "Shared larval host identities in source literature do not prove local shared larval resource consumption",
            "Annual survey counts per effort are not actual adult population sizes, larval recruitment or guaranteed comparable detection",
            "Imputed survey zeros are conditional on all-species repeated transect assumption and species ever detected at site",
            "Negative lag associations may be habitat filtering, weather, management or regression-to-mean",
            "Does not identify historical plant introductions or causally attribute changing butterfly abundance to introduced plants",
            "No permutation/plant-degree null has yet been fitted: inference remains exploratory association ONLY",
        ],
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":output["status"],
        "programmes":[{
            "programme":p["programme"],
            "eligible":p.get("analysis_gate_pass"),
            "n_shared_pairs":p.get("source_host_sharing_pairs"),
            "sites":p.get("sites_eligible"),
            "transitions":p.get("site_species_transitions_with_both_comparator_groups"),
            "beta_lag_shared":p.get("model",{}).get("beta_shared_host_peer_lag"),
            "beta_lag_nonshared":p.get("model",{}).get("beta_nonshared_peer_lag"),
            "contrast":p.get("model",{}).get("contrast_shared_minus_nonshared"),
        }for p in programmes],
    },indent=2))


if __name__=="__main__":
    main()

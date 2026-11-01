#!/usr/bin/env python3
"""Frozen 51-island HOSTS × Lenzner four-empire × UNCTAD 2010 geography tests.

This answers WHERE historical imperial botanical-geography contact is
associated with existing host/butterfly similarities among *port-indexed*
islands, not whether the historical link carried these specific species.

Original host plant indexes span 2001–2024 and butterfly first-index events
span 2011–2022; not a prospective ecological forecast, nor evidence that
ports actually transferred living plants/eggs/larvae.

Source unit is 51 GSHHG landmass nodes × 194 WPI 2019 port geocodes, pre-
defined 2/10/25km TDWG 369 region-history quality gate, 2010 UNCTAD source
country-pair LSBCI. All species and island-pair outcomes are source-frozen.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import loio_pair_baseline as source

BASE=Path(__file__).resolve().parent
COLONIAL_FIELDS=("years_GB","years_ESP","years_PRT","years_NL")
LOCK=BASE/"island_port_colonial_design_lock_20261101.json"
SHIPPING=BASE/"full_lsbci2010_island_pairs.csv"
MODELS={
    "B0_geography_prior_richness_ports":source.EFFORT,
    "B1_plus_shared_hosts":source.EFFORT+("host",),
    "H0_plus_colonial_history":source.EFFORT+("colonial_log","colonial_binary"),
    "H1_history_plus_hosts":source.EFFORT+("colonial_log","colonial_binary","host"),
    "H2_history_hosts_interaction":source.EFFORT+(
        "colonial_log","colonial_binary","host","host_x_history"),
}
MARITIME_MODELS={
    "S0_2010_country_shipping":source.EFFORT+("shipping_country_2010",),
    "S1_plus_colonial_history":source.EFFORT+(
        "shipping_country_2010","colonial_log","colonial_binary"),
    "S2_plus_host_composition":source.EFFORT+(
        "shipping_country_2010","colonial_log","colonial_binary","host"),
    "S3_host_by_colonial_history":source.EFFORT+(
        "shipping_country_2010","colonial_log","colonial_binary","host","host_x_history"),
}


def source_rows(path):
    with Path(path).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def joined_pairs(km=10):
    if km not in (2,10,25):
        raise ValueError("Only response-blind fixed geocoding thresholds allowed")
    if not LOCK.exists():
        raise RuntimeError("Original island port colonial test lock missing")
    islands=source_rows(BASE/f"island_port_colonial_{km}km.islands.csv")
    id_map={r["landmass_id"]:r for r in islands
            if r["eligible_for_history_test"]=="1"}
    expected={2:40,10:42,25:42}
    if len(id_map)!=expected[km]:
        raise RuntimeError(f"Original 51-island colonial-port source support changed at {km}km: {len(id_map)}")
    shipping={}
    for r in source_rows(SHIPPING):
        key=tuple(sorted((r["island_a"],r["island_b"])))
        if key in shipping:
            raise RuntimeError("Repeated country shipping pair")
        shipping[key]=r
    if len(shipping)!=1275:
        raise RuntimeError("Original UNCTAD country-pair source changed")
    _, original=source.load()
    if len(original)!=1275:
        raise RuntimeError("Original 51-island host-butterfly pair source changed")
    new=[]
    for o in original:
        a,b=o["a"],o["b"]
        if a not in id_map or b not in id_map:
            continue
        ah,bh=id_map[a],id_map[b]
        shared=[min(float(ah[col]),float(bh[col])) for col in COLONIAL_FIELDS]
        history=sum(shared)
        binary=int(any(d>=25 for d in shared))
        entry=dict(o)
        entry.update({
            "colonial_log":math.log1p(history),
            "colonial_binary":float(binary),
            "colonial_shared_duration_years":history,
            "host_x_history":o["host"]*math.log1p(history),
            "a_history_coverage":float(ah["WPI_port_history_coverage_fraction"]),
            "b_history_coverage":float(bh["WPI_port_history_coverage_fraction"]),
            "same_current_country":float(
                ah["current_WPI_country_codes"]==bh["current_WPI_country_codes"]),
            "historic_4_empire_share_label":";".join(
                col.removeprefix("years_") for col,v in zip(COLONIAL_FIELDS,shared) if v>=25),
        })
        shipping_row=shipping[tuple(sorted((a,b)))]
        entry["shipping_country_2010"]=(
            float(shipping_row["lsbci_2010"])
            if shipping_row["status"]=="INDEX_2010_OBSERVED" else None)
        entry["shipping_source_status"]=shipping_row["status"]
        new.append(entry)
    expect_pairs=len(id_map)*(len(id_map)-1)//2
    if len(new)!=expect_pairs:
        raise RuntimeError("Unmatched island pair source: "+str(len(new)))
    return sorted(id_map),new,{
        "eligible_history_islands":len(id_map),
        "source_port_total_on_eligible_islands":sum(int(r["WPI_ports_original"]) for r in id_map.values()),
        "source_matched_port_total":sum(int(r["WPI_ports_colonial_mapped"]) for r in id_map.values()),
        "historic_diversity":dict(Counter(
            r["historic_4_empire_share_label"] for r in new)),
        "shared_focal_empire_pairs":sum(r["colonial_binary"]>0 for r in new),
        "source_2010_country_shipping_pairs":sum(r["shipping_country_2010"] is not None for r in new),
        "source_same_current_country_pairs":sum(r["same_current_country"]>0 for r in new),
    }


def ranks(values):
    x=np.asarray(values,float)
    ix=np.argsort(x,kind="stable")
    result=np.empty(len(x),float)
    j=0
    while j<len(x):
        k=j+1
        while k<len(x) and x[ix[k]]==x[ix[j]]:
            k+=1
        result[ix[j:k]]=(j+1+k)/2
        j=k
    return result


def partial_rank(rows,left,right,controls=()):
    if len(rows)<40:
        raise ValueError("Insufficient island-pair geography")
    ranked={key:ranks([r[key] for r in rows]) for key in {left,right,*controls}}
    def residual(v):
        if not controls:
            return v-v.mean()
        z=np.column_stack([np.ones(len(rows))]+[ranked[c] for c in controls])
        return v-z@np.linalg.lstsq(z,v,rcond=None)[0]
    a,b=residual(ranked[left]),residual(ranked[right])
    if a.std()<1e-10 or b.std()<1e-10:
        return None
    return float(np.corrcoef(a,b)[0,1])


def fit_held_island(rows,ids,models,lam=1.):
    if not rows or lam<0:
        raise ValueError("No valid data or ridge penalty")
    summaries={}
    per_island={}
    for model,features in models.items():
        stats=[]
        for island in ids:
            train=[r for r in rows if island not in (r["a"],r["b"])]
            held=[r for r in rows if island in (r["a"],r["b"])]
            if not held or len(train)<200 or len(held)<10:
                continue
            pred=source.fit_predict(train,held,features,lam)
            errors=[p-r["y"] for p,r in zip(pred,held)]
            stats.append({
                "island_id":island,
                "training_pairs":len(train),"test_pairs":len(held),
                "mse":sum(v*v for v in errors)/len(errors),
                "mae":sum(abs(v) for v in errors)/len(errors),
            })
        if len(stats)<20:
            raise RuntimeError("Too few valid island node holdout folds for "+model)
        summaries[model]={
            "island_equal_MSE":sum(t["mse"] for t in stats)/len(stats),
            "island_equal_MAE":sum(t["mae"] for t in stats)/len(stats),
            "island_node_folds":len(stats),
        }
        per_island[model]=stats
    return summaries,per_island


def contrast(summaries,folds,base,added):
    a={r["island_id"]:r for r in folds[base]}
    b={r["island_id"]:r for r in folds[added]}
    shared=sorted(set(a)&set(b),key=int)
    if len(shared)<20:
        raise RuntimeError("Inadequate paired island fold support")
    delta=[a[i]["mse"]-b[i]["mse"] for i in shared]
    amse=np.mean([a[i]["mse"] for i in shared])
    return {
        "island_equal_mse_reduction_percent":float(100*np.mean(delta)/amse),
        "mean_mse_reduction":float(np.mean(delta)),
        "median_mse_reduction":float(np.median(delta)),
        "island_nodes_improved":sum(x>0 for x in delta),
        "island_nodes_worsened":sum(x<0 for x in delta),
        "n_shared_node_folds":len(shared),
        "diff_range":[min(delta),max(delta)],
    }


def analyze(km=10,lam=1.):
    ids,rows,audit=joined_pairs(km)
    original={
        "host_x_new_butterfly_raw":partial_rank(rows,"host","y"),
        "host_x_new_butterfly_after_prior_distance":partial_rank(
            rows,"host","y",("dist","old")),
        "host_x_new_butterfly_after_history_prior_distance":partial_rank(
            rows,"host","y",("dist","old","colonial_log","colonial_binary")),
        "history_x_hosts_raw":partial_rank(rows,"colonial_log","host"),
        "history_x_hosts_after_distance_prior":partial_rank(
            rows,"colonial_log","host",("dist","old")),
        "history_x_hosts_after_distance_prior_catalogue":partial_rank(
            rows,"colonial_log","host",("dist","old","h0","h1","p0","p1")),
        "history_x_new_butterfly_after_distance_prior":partial_rank(
            rows,"colonial_log","y",("dist","old")),
    }
    model_scores,folds=fit_held_island(rows,ids,MODELS,lam)
    comparisons={
        "host_addition_without_colonial_context":contrast(model_scores,folds,
            "B0_geography_prior_richness_ports","B1_plus_shared_hosts"),
        "colonial_context_without_host":contrast(model_scores,folds,
            "B0_geography_prior_richness_ports","H0_plus_colonial_history"),
        "host_addition_given_colonial_context":contrast(model_scores,folds,
            "H0_plus_colonial_history","H1_history_plus_hosts"),
        "colonial_context_given_host":contrast(model_scores,folds,
            "B1_plus_shared_hosts","H1_history_plus_hosts"),
        "host_by_colonial_history_interaction_increment":contrast(model_scores,folds,
            "H1_history_plus_hosts","H2_history_hosts_interaction"),
    }
    tagged=[r for r in rows if r["shipping_country_2010"] is not None]
    shipping={}
    if len(tagged)>=200:
        shipids=sorted({i for r in tagged for i in (r["a"],r["b"])},key=int)
        shipassoc={
            "host_x_new_butterfly_adjust_history_dist_prior_and_shipping":
                partial_rank(tagged,"host","y",(
                    "dist","old","colonial_log","colonial_binary","shipping_country_2010")),
            "history_x_hosts_adjust_distance_prior_and_country_shipping":
                partial_rank(tagged,"colonial_log","host",(
                    "dist","old","shipping_country_2010")),
            "shipping_x_hosts_adjust_distance_prior_and_colonial":
                partial_rank(tagged,"shipping_country_2010","host",(
                    "dist","old","colonial_log","colonial_binary")),
        }
        sscore,sfold=fit_held_island(tagged,shipids,MARITIME_MODELS,lam)
        shipping={
            "status":"COUNTRY_PAIR_LSBCI_ASSOCIATION_NOT_PORT_TO_PORT_OBSERVED_VOYAGES",
            "scored_country_bilateral_pairs":len(tagged),
            "island_nodes":len(shipids),
            "associations":shipassoc,
            "models":sscore,
            "comparisons":{
                "historical_context_after_country_shipping":contrast(sscore,sfold,
                  "S0_2010_country_shipping","S1_plus_colonial_history"),
                "host_after_shipping_and_colonial_context":contrast(sscore,sfold,
                  "S1_plus_colonial_history","S2_plus_host_composition"),
                "history_host_interaction_after_country_shipping":contrast(sscore,sfold,
                  "S2_plus_host_composition","S3_host_by_colonial_history"),
            },
        }
    else:
        shipping={"status":"UNCTAD_COUNTRY_BILATERAL_DATA_TOO_SPARSE","pairs":len(tagged)}
    # Assess source-geographic sample rather than hiding attrition; no iid p.
    return {
        "schema":"chocho_51island_4empire_port_source_observation_association_v1",
        "status":"RETROSPECTIVE_ISLAND_PORT_BOTANICAL_HISTORY_ASSOCIATION_NOT_MEDIATION_OR_SHIP_MOVEMENT",
        "source_sha256":{
            "island_port_colonial_portmean":sha256(BASE/f"island_port_colonial_{km}km.islands.csv"),
            "frozen_island_pairs":sha256(source.PAIRS),
            "frozen_island_nodes":sha256(source.ISLANDS),
            "2010_UNCTAD_LSBCI_country_pair":sha256(SHIPPING),
            "pre_response_lock":sha256(LOCK),
        },
        "max_WPI_port_distance_km":km,
        "ridge_lambda":lam,
        "island_source":audit,
        "island_pair_count":len(rows),
        "response":"2011-2022 newly indexed butterfly assemblage Jaccard across 2001-2024 pooled port-indexed possible alien host plant sources",
        "historical_contact":"sum min of four Lenzner area-weighted empire years across matched 2019-WPI port geographies; binary if >=25 years both in any empire",
        "pair_rank_associations":original,
        "models":model_scores,
        "model_comparisons":comparisons,
        "country_shipping_sensitivity":shipping,
        "interpretation_limits":[
            "Port source coordinates are current World Port Index geocodes. They do NOT identify operating colonial-era ports.",
            "Only historical GB/ESP/PRT/NL domains in Lenzner, not France or other empires; zero four-empire exposure is not absence of colonial history.",
            "Lenzner area-weighted region years are not actual port-to-port shipping edges or documented plant exchange dates.",
            "Original 51-island positive host/butterfly correlation and other analyses were not discarded by selecting the history-geocoded island subset.",
            "2011-2022 newly GBIF-indexed butterfly records are not documented arrival/establishment.",
            "Port-indexed host presence 2001–2024 can postdate first butterfly records and reflects botanical reporting effort.",
            "Catalogued host species richness and port count are not observation-effort-complete controls; general alien-plant composition unavailable on the exact same island/port unit.",
            "Non-independent island pair observations and overlapping LOIO scores are not 51 statistically independent tests.",
            "Contemporary country bilateral liner index is not a 2001 historical voyage or per-port route.",
            "Potential host×shared history interaction is retrospective, exploratory and not proof of mediation.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--max-port-to-region-km",type=int,default=10)
    p.add_argument("--ridge-lambda",type=float,default=1.)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    result=analyze(args.max_port_to_region_km,args.ridge_lambda)
    s=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding="utf8")
    else:
        print(s,end="")


if __name__=="__main__":
    main()

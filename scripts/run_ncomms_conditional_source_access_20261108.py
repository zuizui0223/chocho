#!/usr/bin/env python3
"""Conditional original butterfly GBIF index test; not biological colonisation.

Use existing exact Git-blob-pinned source data, 225 taxa x 17 original islands,
and indexed prior-period conspecific sources. Compare geodesic closeness and
sea resistance residuals using fixed species AND island positive-event counts.
Source outcomes are indexed GBIF records, never verified immigration/breeding.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import importlib.util
import json
from pathlib import Path
import random
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts/run_ncomms_public17_ocean_access_holdout.py"
SEED = 20261108

def load_original():
    spec = importlib.util.spec_from_file_location("frozen_public17_source", SRC)
    if spec is None or spec.loader is None:
        raise RuntimeError("Frozen public source importer unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, *module.load(ROOT)

def switch(edges, event_set, eligible, rng):
    a = rng.randrange(len(edges))
    b = rng.randrange(len(edges))
    if a == b:
        return False
    u, v = edges[a], edges[b]
    if u[0] == v[0] or u[1] == v[1]:
        return False
    x = (u[0], v[1])
    y = (v[0], u[1])
    if x not in eligible or y not in eligible or x in event_set or y in event_set:
        return False
    event_set.difference_update((u, v))
    event_set.update((x, y))
    edges[a], edges[b] = x, y
    return True

def two_margin_null(positives, eligible, lookup, draws=1250, chains=4, burn=18000, thin=100):
    keys = sorted(next(iter(lookup.values())))
    stats = []
    diagnostics = []
    m0s = Counter(s for s, i in positives)
    m0i = Counter(i for s, i in positives)
    for ch in range(chains):
        rng = random.Random(SEED+1709*ch)
        edges = list(positives)
        edge_set = set(edges)
        accept = trials = 0
        while accept < burn:
            accept += int(switch(edges, edge_set, eligible, rng))
            trials += 1
            if trials > burn * 800:
                raise RuntimeError("Insufficient null-switch moves")
        seq = []
        for _ in range(draws):
            n = 0
            while n < thin:
                n += int(switch(edges, edge_set, eligible, rng))
                trials += 1
                if trials > 800*(burn+thin*draws):
                    raise RuntimeError("Null chain stalled")
            seq.append([sum(lookup[e][k] for e in edges)/len(edges) for k in keys])
        if (Counter(s for s,i in edges) != m0s or
                Counter(i for s,i in edges) != m0i):
            raise AssertionError("Null simulation changed species/island margin")
        arr = np.array(seq, float)
        stats.extend(seq)
        diagnostics.append({
            "chain": ch, "accepted_switches": burn + thin*draws,
            "proposed_switches": trials,
            "acceptance_fraction": (burn+thin*draws)/trials,
            "means": dict(zip(keys,arr.mean(axis=0).tolist()))})
    return np.array(stats, float), keys, diagnostics

def island_only_null(rows, positives, lookup, draws=5000):
    rng = np.random.default_rng(SEED+51)
    by_island = defaultdict(list)
    by_event_island = Counter(i for _,i in positives)
    for z in rows:
        s, i = int(z["speciesKey"]), int(z["landmass_id"])
        by_island[i].append((s, i))
    keys = sorted(next(iter(lookup.values())))
    out = []
    for _ in range(draws):
        sampled = []
        for i, pool in by_island.items():
            n = by_event_island[i]
            if n:
                sampled.extend(pool[int(k)] for k in
                               rng.choice(len(pool),size=n,replace=False))
        out.append([sum(lookup[e][k] for e in sampled)/len(sampled) for k in keys])
    return np.array(out,float),keys

def compare(observed, null, name):
    mu = float(np.mean(null))
    sd = float(np.std(null,ddof=1))
    return {
        "null":name, "observed_mean_standardized_access":float(observed),
        "null_mean":mu, "observed_minus_null":float(observed-mu),
        "null_sd":sd, "null_95pct":list(map(float,np.quantile(null,[0.025,0.975]))),
        "one_sided_higher_access_p_mc":float((1+np.count_nonzero(null>=observed))/(len(null)+1)),
        "two_sided_p_mc":float((1+np.count_nonzero(abs(null-mu)>=abs(observed-mu)))/(len(null)+1)),
        "null_draws":int(len(null)),
    }

def analyze(src, data, period):
    rows, meta = src.prepare(data,period)
    counts = defaultdict(set)
    for row in rows:
        counts[int(row["speciesKey"])].add(int(row["prior_source_islands"]))
    assert all(len(v)==1 for v in counts.values()), "Source status varied within species"
    subset = [r for r in rows if int(r["prior_source_islands"])>0]
    events = [(int(r["speciesKey"]),int(r["landmass_id"]))
              for r in subset if int(r["record_index_case"])==1]
    expected = (864,57) if period=="2011_2017" else (1039,17)
    assert (len(subset),len(events))==expected
    eligible = {(int(r["speciesKey"]),int(r["landmass_id"])) for r in subset}
    assert len(eligible)==len(subset) and len(set(events))==len(events)
    geo = np.array([float(r["features"]["min_log_geodesic_km"]) for r in subset])
    sea = np.array([float(r["features"]["min_log_uniform_effective_resistance"]) for r in subset])
    assert np.isfinite(geo).all() and np.isfinite(sea).all()
    Z = np.stack([np.ones(len(geo)),geo,geo**2],axis=1)
    b = np.linalg.lstsq(Z,sea,rcond=None)[0]
    residual = sea - Z@b
    metric = np.column_stack([-geo,-residual])
    sd = np.std(metric,axis=0)
    assert np.all(sd>1e-7)
    metric = (metric-metric.mean(axis=0))/sd
    names = ["geodesic_closeness","sea_residual_access"]
    lookup = {(int(r["speciesKey"]),int(r["landmass_id"])):
               dict(zip(names,vals.tolist())) for r,vals in zip(subset,metric)}
    observed = {k:float(np.mean([lookup[e][k] for e in events])) for k in names}
    island_null, ki = island_only_null(subset,events,lookup)
    fixed_null, kf, chains = two_margin_null(events,eligible,lookup)
    assert ki==kf==names
    tests = {}
    for c,k in enumerate(names):
        tests[k] = {
          "island_only_null":compare(observed[k],island_null[:,c],"fixed island counts only"),
          "species_and_island_null":compare(observed[k],fixed_null[:,c],"fixed species and island counts"),
        }
    return {
      "period":period, "eligible_risk":len(rows),
      "new_first_index_cases":int(sum(r["record_index_case"] for r in rows)),
      "eligible_with_prior_other_island_index":len(subset),
      "first_index_with_prior_other_island_index":len(events),
      "species_prior_source_count_invariant_across_eligible_islands":True,
      "source_positive_species":len(set(r["speciesKey"] for r in subset)),
      "islands":len(set(r["landmass_id"] for r in subset)),
      "sea_residual_fit_source_only":"log sea resistance ~ quadratic(log distance), before outcomes, standardized eligible-risk values",
      "fixed_null_chains":chains,
      "tests":tests,
    }

def synthetic_test():
    eligible={(s,i) for s in range(1,4) for i in range(1,4)}
    events=[(1,1),(2,2),(3,3)]
    a=Counter(x for x,y in events); b=Counter(y for x,y in events)
    rng=random.Random(42)
    seq=list(events); st=set(seq)
    for _ in range(20000):
        switch(seq,st,eligible,rng)
    assert len(st)==len(seq)==3 and st==set(seq)
    assert Counter(x for x,y in seq)==a and Counter(y for x,y in seq)==b

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    synthetic_test()
    if args.self_test:
        print("PASS synthetic two-margin switch preservation")
        return
    if args.out_dir is None:
        ap.error("--out-dir required")
    src, data, pins = load_original()
    result = {
      "schema":"chocho_conditional_prior_source_access_20261108_v1",
      "status":"EXPLORATORY_POST_HOC_GBIF_INDEXED_RECORD_ASSOCIATION",
      "source_pins":pins,
      "specification":"Prior indexed same-species sources only; fixed positive events per species and island. Observed higher standardized proximity/access tested against 5000 draws.",
      "fixed_margin_method":"2x2 MCMC event-location swaps over only eligible original species x island cells, 4 chains of 1250 draws; 18000 accepted burn-in, 100 accepted switches between draws",
      "not_measured":["actual arrival","biological colonisation","reproduction","historical butterfly or ship movement","species-specific earlier host plant introduction"],
      "periods":{p:analyze(src,data,p) for p in ("2011_2017","2018_2022")}}
    args.out_dir.mkdir(parents=True,exist_ok=True)
    out = args.out_dir/"conditional_source_access_results.json"
    out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"file":str(out),"results":{
        p:{name:{k:{"excess":v["observed_minus_null"],
                       "p_one_sided":v["one_sided_higher_access_p_mc"]}
                       for k,v in v.items()}
             for name,v in result["periods"][p]["tests"].items()}
        for p in result["periods"]}},indent=2))
if __name__=="__main__":
    main()

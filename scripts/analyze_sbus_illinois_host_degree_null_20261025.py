#!/usr/bin/env python3
"""Family-stratified bipartite host-link degree null for US Illinois adult counts.

This is a POSTHOC, SOURCE-CONDITIONAL structural negative control, not a
butterfly competition experiment, never a population causal estimate.
No site IDs or butterfly individual observations are persisted.
"""
from __future__ import annotations
import argparse
import bisect
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import random
import tempfile
import time

import numpy as np

from analyze_sbus_illinois_niche_counts_20261025 import (
    PROGRAMS, calc_panel, frozen_focal, group_within_matrix, host_sets,
    load_source, read_events, regression,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCE="10.6084/m9.figshare.27934602.v1"
HOSTS_SHA="808e0b869f9ec1adf8efff87cf6a395adda103e0"
WCVP_SHA="65bed76bae9d644ccb6ad200c05f9f5071d89e05"
SEED=20261025
DRAWS=199
SWAP_FACTOR=5
OBSERVED_CONTRAST=-0.28350714681067035
OBSERVED_N=9178


def butterfly_families(focal):
    path=ROOT/"data/external/leptraits_consensus_v1.0.csv"
    fam={}
    with path.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f)
        names=set(reader.fieldnames or [])
        if not {"Species","Family"}<=names:
            raise ValueError(f"Original LepTraits lacks taxonomic family fields: {names}")
        for row in reader:
            sp=row["Species"].strip()
            fa=row["Family"].strip()
            if sp in focal:
                if not fa or (sp in fam and fam[sp]!=fa):
                    raise ValueError(f"Conflicting or absent butterfly family: {sp}")
                fam[sp]=fa
    if len(fam)!=239:
        raise ValueError(f"Expected all 239 butterflies assigned exact families; got {len(fam)}")
    return fam


def accepted_edges(path,focal,taxon):
    original_hostset=host_sets(path,focal)
    plantfam={}
    edge_groups=defaultdict(set)
    original=set()
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        need={"insect_species","accepted_plant_name_id","family"}
        if not need<=set(reader.fieldnames or []):
            raise ValueError("Original WCVP/HOSTS sidecar lacks source plant-family labels")
        for row in reader:
            sp=row["insect_species"].strip()
            if sp not in focal:
                continue
            host=row["accepted_plant_name_id"].strip()
            family=row["family"].strip()
            if not host or not family:
                raise ValueError("Unresolved accepted host/family in source graph")
            if host in plantfam and plantfam[host]!=family:
                raise ValueError("Accepted host has inconsistent plant-family designation")
            plantfam[host]=family
            original.add((sp,host))
            edge_groups[(taxon[sp],family)].add((sp,host))
    if {k:frozenset(h for sp,h in original if sp==k) for k in focal}!=original_hostset:
        raise ValueError("Original host pairs do not match primary analysed source")
    edges={group:sorted(rows) for group,rows in edge_groups.items()}
    viable={k:v for k,v in edges.items() if len(v)>=2 and
            len({a for a,b in v})>=2 and len({b for a,b in v})>=2}
    return original_hostset, original,edges,viable,plantfam


def draw_graph(original_hostset, original_pairs, all_groups, viable, plantfam, taxon, replicate):
    # Every draw starts afresh from the same source-locked network, independent RNG.
    host_map={s:set(links) for s,links in original_hostset.items()}
    connections=set(original_pairs)
    groups={k:list(v) for k,v in viable.items()}
    keys=sorted(groups)
    lengths=[len(groups[k]) for k in keys]
    weight_cumulative=[]
    total=0
    for n in lengths:
        total+=n
        weight_cumulative.append(total)
    E=len(connections)
    nviable=sum(lengths)
    if nviable < 100:
        raise ValueError("Insufficient swappable host relationships")
    accepted=0
    rng=random.Random(SEED+7919*(replicate+1))
    for _ in range(SWAP_FACTOR*E):
        pos=rng.randrange(nviable)
        index=bisect.bisect_right(weight_cumulative,pos)
        group=groups[keys[index]]
        start=weight_cumulative[index-1] if index else 0
        idx=pos-start
        other=rng.randrange(len(group))
        if idx==other:
            continue
        a,x=group[idx]
        b,y=group[other]
        if a==b or x==y or (a,y) in connections or (b,x) in connections:
            continue
        connections.remove((a,x))
        connections.remove((b,y))
        connections.add((a,y))
        connections.add((b,x))
        host_map[a].remove(x)
        host_map[a].add(y)
        host_map[b].remove(y)
        host_map[b].add(x)
        group[idx]=(a,y)
        group[other]=(b,x)
        accepted+=1
    final=set(connections)
    moved_links=len(original_pairs-final)
    if len(final)!=len(original_pairs):
        raise ValueError("Fixed graph edge count drift")
    # All butterfly degrees and plant per-family usage must be preserved.
    if any(len(host_map[sp])!=len(original_hostset[sp]) for sp in original_hostset):
        raise ValueError("Butterfly host degree not preserved")
    originally=Counter((host for sp,host in original_pairs))
    after=Counter((host for sp,host in final))
    if after!=originally:
        raise ValueError("Plant popularity/degree not preserved")
    # Explicitly audit taxonomic butterfly family × botanical plant family
    # strata after every independently rebuilt source graph.
    origin_strata=Counter((taxon[sp],plantfam[host]) for sp,host in original_pairs)
    final_strata=Counter((taxon[sp],plantfam[host]) for sp,host in final)
    if final_strata!=origin_strata:
        raise ValueError("Butterfly family by botanical host family margins changed")
    origin_per_butterfly_family=Counter((sp,plantfam[host]) for sp,host in original_pairs)
    final_per_butterfly_family=Counter((sp,plantfam[host]) for sp,host in final)
    if origin_per_butterfly_family!=final_per_butterfly_family:
        raise ValueError("Consumer botanical family host portfolio drift")
    return {sp:frozenset(hs) for sp,hs in host_map.items()},{
        "accepted":accepted,
        "moved_original_edges_fraction":moved_links/max(1,E),
        "moved_viable_edges_fraction":moved_links/max(1,nviable),
        "edges_total":E,
        "viable_edges":nviable,
    }


def slope(obs):
    X,y,groups=group_within_matrix(obs)
    n,k=X.shape
    if n<300 or n<=k+100:
        return None
    beta,_,rank,_=np.linalg.lstsq(X,y,rcond=None)
    if rank<k or not np.isfinite(beta).all():
        return None
    return float(beta[1]-beta[2])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--accepted-host-csv",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    focal=frozen_focal()
    families=butterfly_families(focal)
    hosts,original_pairs,groups,viable,plantfam=accepted_edges(a.accepted_host_csv,focal,families)
    with tempfile.TemporaryDirectory() as temp:
        source=Path(temp)/"original"
        load_source(source)
        events,source_count=read_events(source,focal)
    observed_report,observed_obs=calc_panel(PROGRAMS[0],events,focal,hosts)
    beta=regression(observed_obs)
    if not math.isclose(beta["contrast_shared_minus_nonshared"],OBSERVED_CONTRAST,abs_tol=1e-9) or len(observed_obs)!=OBSERVED_N:
        raise ValueError("Original Illinois source cohort or coefficient drift")
    output={
        "schema":"chocho_illinois_bipartite_host_degree_null_20261025_v01",
        "dataset_doi":SOURCE,
        "status":"NOT_EVALUABLE_UNTIL_SOURCE_AND_CHAIN_QC",
        "source_239_species":239,
        "observed_beta_shared_minus_nonshared":OBSERVED_CONTRAST,
        "observed_n":OBSERVED_N,
        "original_shared_host_pairs":observed_report["source_host_sharing_pairs"],
        "plant_species_link_edges":len(original_pairs),
        "strata_all":len(groups),
        "strata_viable":len(viable),
        "source_preserves":[
            "Every accepted butterfly-to-host link count and botanical host family per butterfly",
            "Every host plant's global butterfly degree",
            "Every butterfly family x botanical plant family edge stratum",
            "Original adult survey counts, source visit/site/year and duration records unchanged",
        ],
        "do_not_interpret_as_causality":True,
    }
    samples=[]
    failures=[]
    started=time.monotonic()
    for i in range(DRAWS):
        randomized,diag=draw_graph(hosts,original_pairs,groups,viable,plantfam,families,i)
        report,ob=calc_panel(PROGRAMS[0],events,focal,randomized)
        coef=slope(ob) if report["analysis_gate_pass"] else None
        item={
            "replicate":i,
            "beta":coef,
            "n":len(ob),
            "n_sites":report["sites_eligible"],
            "shared_pairs":report["source_host_sharing_pairs"],
            "accepted_swaps":diag["accepted"],
            "moved_viable":diag["moved_viable_edges_fraction"],
        }
        samples.append(item)
        if (coef is None or diag["accepted"]<100 or diag["moved_viable_edges_fraction"]<.10
            or abs(len(ob)-OBSERVED_N)/OBSERVED_N>.25):
            failures.append(i)
        if (i+1)%25==0:
            print("NULL_REPLAY",i+1,"of",DRAWS,
                  "current_beta",coef,"bad",len(failures),
                  "elapsed_s",round(time.monotonic()-started,1),flush=True)
    eligible=[s["beta"] for s in samples if s["replicate"] not in failures]
    qualification=len(eligible)>=190 and len(failures)<=9
    result={
        "status":("POSTHOC_SOURCE_CONDITIONED_NETWORK_NULL_QUALIFIED"
                  if qualification else "NOT_EVALUABLE_DEGREE_OR_COHORT_NULL_GATE"),
        "total_requested":DRAWS,
        "complete_numerically_valid":len(eligible),
        "quality_failed_replicates_n":len(failures),
        "accepted_swaps_range":[min(s["accepted_swaps"] for s in samples),
                                max(s["accepted_swaps"] for s in samples)],
        "changed_viable_edge_fraction_range":[min(s["moved_viable"] for s in samples),
                                             max(s["moved_viable"] for s in samples)],
        "valid_transition_count_range":[min(s["n"] for s in samples),
                                        max(s["n"] for s in samples)],
        "host_sharing_pair_count_range":[min(s["shared_pairs"] for s in samples),
                                         max(s["shared_pairs"] for s in samples)],
        "independent_random_graphs":DRAWS,
        "random_seed_start":SEED,
        "observed_beta":OBSERVED_CONTRAST,
        "null_beta_95pct_and_median":(
            np.quantile(eligible,[.025,.5,.975],method="linear").tolist()
            if qualification else None),
        "left_tail_empirical_fraction_exploratory_only":(
            (1+sum(b<=OBSERVED_CONTRAST for b in eligible))/(len(eligible)+1)
            if qualification else None),
        "null_beta_history":eligible if qualification else None,
        "failed_null_replicate_indexes":failures,
        "not_an_experimental_competition_effect":True,
        "interpretation":"Posthoc host-graph degree structure falsification; adult observational data and imputed zeros remain major limitations.",
    }
    output["result"]=result
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("HOST_GRAPH_NULL_FINAL",json.dumps({
        "status":result["status"],
        "obs":OBSERVED_CONTRAST,
        "ci":result["null_beta_95pct_and_median"],
        "failed":len(failures),
        "changed_edges":result["changed_viable_edge_fraction_range"],
    }),flush=True)


if __name__=="__main__":
    main()

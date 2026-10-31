#!/usr/bin/env python3
"""Global matched native/current flora -> butterfly larval-resource beta diversity.

Inputs are EXACT pinned HOSTS + WCVP-v13 sidecars, already used in chocho GEB.
Outputs represent potential HOST RESOURCES of the fixed 239 butterflies.
They DO NOT represent butterfly occupancy, historical plant introduction dates
or the effect of colonial trade.

Primary: all accepted vascular flora as comparison against a 239-butterfly
fixed-diet resource-availability assemblage at the SAME WGSRPD level-3 units.
Null: rewire recorded host edges within plant families, preserving BOTH
butterfly host counts and each plant's consumer degree, as well as all flora
geographies, then compare mean signed butterfly-beta shifts. Randomization is
a diagnostic, not evidence of independently randomized historical introductions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import random

import numpy as np


EXPECTED_NATIVE_UNITS = 26530
EXPECTED_CURRENT_UNITS = 41083
EXPECTED_CONSUMERS = 239
SOURCE = {
    "hosts_commit": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
    "wcvp_commit": "65bed76bae9d644ccb6ad200c05f9f5071d89e05",
    "wcvp_region_level": "WGSRPD3",
}


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f)


def checksum(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def jaccard_bits(a, b):
    union = (a | b).bit_count()
    return ((a & b).bit_count() / union) if union else None


def rank(values):
    arr = np.asarray(values, float)
    ix = np.argsort(arr, kind="stable")
    out = np.empty(len(ix), float)
    i = 0
    while i < len(ix):
        j = i + 1
        while j < len(ix) and arr[ix[j]] == arr[ix[i]]:
            j += 1
        out[ix[i:j]] = (i + j + 1) / 2
        i = j
    return out


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3:
        return None
    x, y = rank(a[valid]), rank(b[valid])
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def load_panel(panel):
    x = list(read_rows(panel))
    species = [r["species"].strip() for r in x]
    if len(species) != EXPECTED_CONSUMERS or len(set(species)) != len(species):
        raise RuntimeError(f"expected frozen 239 butterfly names, got {len(species)}")
    return species


def load_interactions(interactions, panel):
    index = {n: i for i, n in enumerate(panel)}
    host_family = {}
    edges = set()
    unresolved = Counter()
    for r in read_rows(interactions):
        b = r["insect_species"].strip()
        if b not in index:
            continue
        p = r["accepted_plant_name_id"].strip()
        family = r["family"].strip()
        if not p or not family:
            unresolved["missing_accepted_host_or_family"] += 1
            continue
        if p in host_family and host_family[p] != family:
            raise RuntimeError("host plant has conflicting family: " + p)
        host_family[p] = family
        edges.add((index[b], p))
    if any(not any(e[0] == i for e in edges) for i in range(len(panel))):
        raise RuntimeError("fixed panel contains butterfly with zero mapped hosts")
    return sorted(edges), host_family, dict(unresolved)


def load_region_data(native_path, current_path, edges):
    # First scan contemporary to determine exact categorical encoding.
    plants = set()
    region_codes = set()
    for row in read_rows(current_path):
        p = row["accepted_plant_name_id"].strip()
        r = row["area_code_l3"].strip()
        if p and r:
            plants.add(p); region_codes.add(r)
    if not (1000 <= len(plants) <= 500000):
        raise RuntimeError(f"unexpected number of accepted vascular plant IDs {len(plants)}")
    regions = sorted(region_codes)
    if not (250 <= len(regions) <= 450):
        raise RuntimeError(f"unexpected WGSRPD3 region coverage {len(regions)}")
    pi = {name: i for i, name in enumerate(sorted(plants))}
    ri = {code: i for i, code in enumerate(regions)}
    width = (len(pi) + 7)//8
    bytearrays = {kind: [bytearray(width) for _ in regions] for kind in ("native","current")}
    plant_reg = {kind: defaultdict(int) for kind in ("native", "current")}
    hostset = {p for _,p in edges}
    totals = {}
    for kind, path in (("native",native_path), ("current",current_path)):
        seen = set() if kind=="native" else None
        rows = 0
        for line in read_rows(path):
            p, region = line["accepted_plant_name_id"].strip(), line["area_code_l3"].strip()
            if p not in pi or region not in ri:
                raise RuntimeError(f"{kind} plant/region missing from contemporary: {p}/{region}")
            j = ri[region]
            t = pi[p]
            bytearrays[kind][j][t >> 3] |= 1 << (t & 7)
            if p in hostset:
                plant_reg[kind][p] |= 1 << j
            rows += 1
            if seen is not None:
                seen.add((p,region))
        totals[kind] = rows
        if seen is not None and len(seen) != rows:
            raise RuntimeError("native source has duplicate accepted plant×region rows")
    flora_bits = {k:[int.from_bytes(x, byteorder="little")
                     for x in v] for k,v in bytearrays.items()}
    for i in range(len(regions)):
        if flora_bits["native"][i] & ~flora_bits["current"][i]:
            raise RuntimeError("native plant flora not subset of current at region " + regions[i])
    # HOSTS may include accepted WCVP host taxa with no non-extinct,
    # non-doubtful current regional distribution. Keep those edges as
    # unresolved geographic opportunity, rather than inventing presence.
    return regions, flora_bits, plant_reg, totals, len(plants)


def load_exact_host_footprints(native_host_path, current_host_path, edges, regions):
    """Use the original HOSTS-matched *host* sidecars for the resource layer.

    Full WCVP all-flora sheets include only accepted-rank species; some
    HOSTS-resolved accepted IDs are not represented in that independently
    filtered flora table. Merging both layers on all-flora IDs would silently
    discard mapped host regions, producing wrong GEB totals.
    """
    ri = {name:i for i,name in enumerate(regions)}
    hosts = {p for _,p in edges}
    host_reg = {"native": defaultdict(int), "current": defaultdict(int)}
    counts = {}
    for kind, path in (("native",native_host_path), ("current",current_host_path)):
        total = 0
        used = 0
        out_of_region = Counter()
        for row in read_rows(path):
            p = row["accepted_plant_name_id"].strip()
            unit = row["area_code_l3"].strip()
            total += 1
            if p not in hosts:
                continue
            if unit not in ri:
                out_of_region[unit] += 1
                continue
            host_reg[kind][p] |= 1 << ri[unit]
            used += 1
        if out_of_region:
            raise RuntimeError("HOSTS sidecar has regions absent from all-flora WGSRPD3: "+
                               repr(dict(out_of_region)))
        counts[kind] = {"total_host_distribution_rows":total,
                        "focal_host_distribution_rows":used,
                        "focal_hosts_with_distribution":len(host_reg[kind])}
    for p in hosts:
        if host_reg["native"].get(p,0) & ~host_reg["current"].get(p,0):
            raise RuntimeError("Host-native footprint outside host-contemporary footprint for "+p)
    return host_reg, counts


def build_butterfly_by_region(edges, plant_reg, n_regions, n_consumers):
    """Turn each butterfly's union of accepted host ranges into region bitsets."""
    footprints = {kind: [0] * n_consumers for kind in ("native","current")}
    for i, p in edges:
        for kind in footprints:
            footprints[kind][i] |= plant_reg[kind].get(p, 0)
    region = {kind:[0]*n_regions for kind in footprints}
    for kind in footprints:
        for butterfly, footprint in enumerate(footprints[kind]):
            active = footprint
            bit = 1 << butterfly
            while active:
                lowest = active & -active
                r = lowest.bit_length() - 1
                region[kind][r] |= bit
                active ^= lowest
    return region


def count_species_units(region_bits):
    return sum(x.bit_count() for x in region_bits)


def pair_rows(regions, flora, resource):
    pairs = []
    for i,a in enumerate(regions):
        for j in range(i+1,len(regions)):
            b = regions[j]
            pn = jaccard_bits(flora["native"][i],flora["native"][j])
            pc = jaccard_bits(flora["current"][i],flora["current"][j])
            rn = jaccard_bits(resource["native"][i],resource["native"][j])
            rc = jaccard_bits(resource["current"][i],resource["current"][j])
            pairs.append({
                "region_a":a,"region_b":b,
                "plant_native_jaccard":pn,"plant_current_jaccard":pc,
                "plant_delta":pc-pn if pn is not None and pc is not None else None,
                "resource_native_jaccard":rn,"resource_current_jaccard":rc,
                "resource_delta":rc-rn if rn is not None and rc is not None else None,
            })
    return pairs


def describe_delta(vals):
    x=np.asarray([v for v in vals if v is not None and math.isfinite(v)],float)
    if len(x)==0:
        return {"n_valid":0}
    tol = 1e-12
    return {
        "n_valid":len(x), "mean":float(np.mean(x)), "median":float(np.median(x)),
        "negative_fraction":float(np.mean(x< -tol)),
        "zero_fraction":float(np.mean(abs(x)<=tol)),
        "positive_fraction":float(np.mean(x>tol)),
        "q025":float(np.quantile(x,.025)),
        "q975":float(np.quantile(x,.975)),
        "note":"Pairs reuse regions; distribution across pairs is not a confidence interval",
    }


def mean_resource_delta(resource):
    n=len(resource["native"])
    total=0.0
    valid=0
    for i in range(n):
        for j in range(i+1,n):
            a=jaccard_bits(resource["native"][i],resource["native"][j])
            b=jaccard_bits(resource["current"][i],resource["current"][j])
            if a is not None and b is not None:
                total+=b-a
                valid+=1
    return (total/valid if valid else None), valid


def mean_resource_delta_common_region_set(resource, selected, min_pairs=100):
    """Compare pairs drawn from the SAME native-resource-supported regions.

    If every selected region has nonempty native and contemporary
    assemblages, every induced pair has a defined Jaccard value. This avoids
    comparing nulls over subtly different lists of empty-vs-empty pairs.
    """
    total=0.
    count=0
    for ai,i in enumerate(selected):
        if resource["native"][i]==0 or resource["current"][i]==0:
            raise RuntimeError("Common-support region missing native/current resources")
        for j in selected[ai+1:]:
            x=jaccard_bits(resource["native"][i],resource["native"][j])
            y=jaccard_bits(resource["current"][i],resource["current"][j])
            if x is None or y is None:
                raise RuntimeError("Unexpected undefined Jaccard on common support")
            total+=y-x
            count+=1
    if count<min_pairs:
        raise RuntimeError("No independent comparable common-support region-pair set")
    return total/count,count


def swap_family_edges(edges, family, rng, attempts):
    # Reversible degree-preserving bipartite rewiring within host families.
    pool=defaultdict(list)
    active=list(edges)
    for ix,(b,p) in enumerate(active):
        pool[family[p]].append(ix)
    byb=set(active)
    families=[v for v in pool.values() if len(v)>=2]
    weights=[len(v) for v in families]
    accepted=0
    for _ in range(attempts):
        ix_group=rng.choices(range(len(families)),weights=weights,k=1)[0]
        group=families[ix_group]
        u,v=rng.sample(group,2)
        b1,p1=active[u];b2,p2=active[v]
        if b1==b2 or p1==p2 or (b1,p2) in byb or (b2,p1) in byb:
            continue
        byb.remove((b1,p1));byb.remove((b2,p2))
        byb.add((b1,p2));byb.add((b2,p1))
        active[u]=(b1,p2);active[v]=(b2,p1)
        accepted+=1
    return active,accepted


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--butterfly-panel",type=Path,required=True)
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--allflora-native-csv",type=Path,required=True)
    ap.add_argument("--allflora-current-csv",type=Path,required=True)
    ap.add_argument("--native-host-distribution-csv",type=Path,required=True)
    ap.add_argument("--current-host-distribution-csv",type=Path,required=True)
    ap.add_argument("--official-wgsrpd3-lookup-json",type=Path,required=True)
    ap.add_argument("--null-draws",type=int,default=199)
    ap.add_argument("--output-prefix",type=Path,required=True)
    args=ap.parse_args()
    if not (0 <= args.null_draws <= 999):
        ap.error("--null-draws must be between 0 and 999")
    panel=load_panel(args.butterfly_panel)
    edges,family,unresolved=load_interactions(args.insect_host_csv,panel)
    regions,flora,flora_mapped_host_regions,flora_rows,n_plants=load_region_data(
        args.allflora_native_csv,args.allflora_current_csv,edges)
    host_regions,host_distribution_rows=load_exact_host_footprints(
        args.native_host_distribution_csv,args.current_host_distribution_csv,
        edges,regions)
    resources=build_butterfly_by_region(edges,host_regions,len(regions),len(panel))
    sums={k:count_species_units(resources[k]) for k in ("native","current")}
    masked_resources=build_butterfly_by_region(
        edges,flora_mapped_host_regions,len(regions),len(panel))
    allflora_masked_sums={
        k:count_species_units(masked_resources[k]) for k in ("native","current")}
    if sums != {"native":EXPECTED_NATIVE_UNITS, "current":EXPECTED_CURRENT_UNITS}:
        raise RuntimeError("GEB source checksum failed using EXACT host distribution: "+
                           repr(sums)+"; all-flora masked "+
                           repr(allflora_masked_sums)+
                           "; required native=26530/current=41083")
    official=json.loads(args.official_wgsrpd3_lookup_json.read_text(encoding="utf8"))
    if official.get("schema") != "chocho_wgsrpd3_level1_lookup_v0.1":
        raise RuntimeError("Official WGSRPD3 lookup source changed")
    official_codes=sorted(official.get("units",{}))
    if len(official_codes)!=369 or len(set(official_codes))!=369:
        raise RuntimeError("Expected fixed 369 official level-3 codes")
    index={code:i for i,code in enumerate(regions)}
    extra=sorted(set(regions)-set(official_codes))
    absent=sorted(set(official_codes)-set(regions))
    if len(extra)>5 or len(absent)>5:
        raise RuntimeError("Large WGSRPD3 source mismatch: extraneous="+repr(extra)+
                           "; unrepresented="+repr(absent))
    official_indices=[index.get(code) for code in official_codes]
    def to_official(seq):
        return [seq[i] if i is not None else 0 for i in official_indices]
    flora_official={k:to_official(flora[k]) for k in ("native","current")}
    resources_official={k:to_official(resources[k]) for k in ("native","current")}
    official_resource_units={k:count_species_units(resources_official[k])
                             for k in ("native","current")}
    # Null uses exactly the same 369 valid codes; absent communities are
    # empty (undefined pairwise Jaccard), not fabricated species presences.
    pairs=pair_rows(official_codes,flora_official,resources_official)
    plant=describe_delta([p["plant_delta"] for p in pairs])
    consumer=describe_delta([p["resource_delta"] for p in pairs])
    subset=[p for p in pairs if p["plant_delta"] is not None
            and p["resource_delta"] is not None]
    rho=spearman([p["plant_delta"] for p in subset],
                 [p["resource_delta"] for p in subset])
    # No claim that raw Jaccard differences between different taxon universes
    # are comparable as effect sizes: the null tests the *consumer mapping*.
    null=[]
    accepts=[]
    null_assemblages=[]
    null_units=[]
    null_valid_pairs=[]
    rng=random.Random(31031)
    if args.null_draws:
        if not edges:
            raise RuntimeError("zero insect-host interactions")
        for i in range(args.null_draws):
            e,accepted=swap_family_edges(edges,family,rng,attempts=5*len(edges))
            m_full=build_butterfly_by_region(e,host_regions,len(regions),len(panel))
            m={kind:[m_full[kind][j] if j is not None else 0
                     for j in official_indices]
               for kind in ("native","current")}
            x,n=mean_resource_delta(m)
            if x is None:
                raise RuntimeError("null all pairs undefined")
            null.append(x);accepts.append(accepted)
            null_assemblages.append(m)
            null_valid_pairs.append(n)
            null_units.append({
                "native":count_species_units(m["native"]),
                "current":count_species_units(m["current"]),
            })
        if np.median(accepts)<len(edges)*.1:
            raise RuntimeError("too few valid within-family swaps; null is too constrained")
    # Common region-level support, fixed once across the observed network
    # and ALL degree-preserving null draws; protects against changes in
    # the number of undefined empty-empty region pairs.
    common_indices=[
        i for i in range(len(official_codes))
        if resources_official["native"][i] and resources_official["current"][i]
        and all(m["native"][i] and m["current"][i] for m in null_assemblages)
    ]
    common_observed=None
    common_null=[]
    common_pair_count=0
    if args.null_draws:
        if len(common_indices)<250:
            raise RuntimeError("Too few shared nonempty native/current regions for common-support null: "+
                               str(len(common_indices)))
        common_observed,common_pair_count=mean_resource_delta_common_region_set(
            resources_official,common_indices)
        for m in null_assemblages:
            v,n=mean_resource_delta_common_region_set(m,common_indices)
            if n!=common_pair_count:
                raise RuntimeError("Common support pair count drift")
            common_null.append(v)
    observed=consumer.get("mean")
    null_summary={
        "n_draws":args.null_draws,
        "scope":"fixed exact plant geographies; within-family edge-swaps preserve each insect and host-plant degree",
        "observed_mean_consumer_delta":observed,
        "median_valid_swaps":float(np.median(accepts)) if accepts else None,
        "median":float(np.median(null)) if null else None,
        "q025":float(np.quantile(null,.025)) if null else None,
        "q975":float(np.quantile(null,.975)) if null else None,
        "p_greater_equal":(1+sum(x>=observed for x in null))/(len(null)+1) if null else None,
        "p_less_equal":(1+sum(x<=observed for x in null))/(len(null)+1) if null else None,
        "null_all_valid_pair_count_min_max":[min(null_valid_pairs),max(null_valid_pairs)] if null_valid_pairs else None,
        "common_support_regions":len(common_indices) if args.null_draws else None,
        "common_support_pair_count":common_pair_count if args.null_draws else None,
        "common_support_observed_mean":common_observed,
        "common_support_null_median":float(np.median(common_null)) if common_null else None,
        "common_support_null_q025":float(np.quantile(common_null,.025)) if common_null else None,
        "common_support_null_q975":float(np.quantile(common_null,.975)) if common_null else None,
        "common_support_p_low":(1+sum(x<=common_observed for x in common_null))/(len(common_null)+1) if common_null else None,
        "common_support_p_high":(1+sum(x>=common_observed for x in common_null))/(len(common_null)+1) if common_null else None,
        "null_native_resource_units_mean":float(np.mean([x["native"] for x in null_units])) if null_units else None,
        "null_current_resource_units_mean":float(np.mean([x["current"] for x in null_units])) if null_units else None,
        "null_added_resource_units_mean":float(np.mean([x["current"]-x["native"] for x in null_units])) if null_units else None,
        "note":"Post-hoc MCMC-style degree-preserving network rewiring. Draws share fixed source; not independent biological replicates.",
    }
    report={
        "schema":"chocho_global_WGSRPD3_flora_consumer_beta_v1",
        "status":"SOURCE_LOCKED_POTENTIAL_LARVAL_RESOURCE_COMPOSITION_NOT_OCCUPANCY",
        "source_commit_pins":SOURCE,
        "source_sha256":{k:checksum(v) for k,v in {
          "butterfly_panel":args.butterfly_panel,
          "insect_host_csv":args.insect_host_csv,
          "allflora_native_csv":args.allflora_native_csv,
          "allflora_current_csv":args.allflora_current_csv,
          "native_host_distribution_csv":args.native_host_distribution_csv,
          "current_host_distribution_csv":args.current_host_distribution_csv,
          "official_wgsrpd3_lookup_json":args.official_wgsrpd3_lookup_json}.items()},
        "analysis_units":{"wgsrpd3_regions":len(official_codes),
          "raw_wcvp_region_codes":len(regions),
          "unexpected_wcvp_codes_not_in_official_369":extra,
          "official_codes_without_contemporary_flora":absent,
          "official_369_resource_units":official_resource_units,
          "all_region_pairs":len(pairs),"selected_butterflies":len(panel),
          "distinct_accepted_flora_species":n_plants,
          "resolved_panel_insect_host_edges":len(edges),
          "resolved_panel_host_plants":len({p for _,p in edges}),
          "allflora_plant_region_rows":flora_rows,
          "host_sidecar_distribution_rows":host_distribution_rows,
          "allflora_masked_resource_units":allflora_masked_sums,
          "exact_host_vs_allflora_masked_gap":{
            k:sums[k]-allflora_masked_sums[k] for k in ("native","current")},
          "native_butterfly_resource_units":sums["native"],
          "current_butterfly_resource_units":sums["current"],
          "added_butterfly_resource_units":sums["current"]-sums["native"],
          "added_percent":100*(sums["current"]-sums["native"])/sums["native"],
          "unresolved_host_links":unresolved,
        },
        "plant_flora_beta_delta":plant,"butterfly_resource_beta_delta":consumer,
        "plant_delta_vs_resource_delta_spearman":rho,
        "paired_valid_region_pairs":len(subset),
        "network_structure_null":null_summary,
        "interpretation_limits":[
          "Contemporary means WCVP current native+introduced, not a dated historical year.",
          "Native-only is a counterfactual subset of present-day mapped ranges, not observed pre-colonial flora.",
          "Neither butterfly range occupancy nor realized larval consumption is inferred.",
          "Raw Jaccard shifts across flora and a 239-consumer resource set are not a scale-comparable amplification ratio.",
          "The null keeps plant distributions and family/degree structure, but not all observational/phylogenetic biases.",
          "All pairs reuse regions; do not calculate naïve independent-pair p-values.",
          "Historical empires, trade routes, climate, and GBIF butterfly observations have NOT been fitted.",
          "WGSRPD3 regions are heterogeneous in size and not all are independent.",
          "WCVP distribution may contain codes outside the official frozen 369; only official 369 appear in pair results.",
          "GEB 26530 and 41083 source checks are on unfiltered original HOSTS sidecars, not the official-region-masked subset.",
        ],
    }
    output=args.output_prefix
    output.parent.mkdir(parents=True,exist_ok=True)
    (output.with_suffix(".json")).write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with gzip.GzipFile(filename="",mode="wb",fileobj=output.with_suffix(".pairs.csv.gz").open("wb"),mtime=0,compresslevel=7) as f:
        wrap=io.TextIOWrapper(f,encoding="utf8",newline="")
        w=csv.DictWriter(wrap,fieldnames=list(pairs[0]))
        w.writeheader();w.writerows(pairs);wrap.flush()
        wrap.detach()
    if args.null_draws:
        with output.with_suffix(".null.csv").open("w",encoding="utf8",newline="") as f:
            w=csv.writer(f)
            w.writerow(["draw","accepted_within_family_degree_swaps",
                        "mean_consumer_delta_over_each_null_valid_pairs",
                        "each_null_valid_pairs","mean_delta_fixed_common_region_pairs",
                        "fixed_common_region_pairs","null_native_resource_units",
                        "null_current_resource_units"])
            for k, val in enumerate(null):
                w.writerow([k+1,accepts[k],val,null_valid_pairs[k],
                            common_null[k],common_pair_count,
                            null_units[k]["native"],null_units[k]["current"]])
    # Reusable source-locked butterfly resource opportunity states, NOT
    # GBIF butterfly occurrence. Species-column ordering is in the receipt;
    # a 239-bit hex vector per official region avoids distributing HOSTS
    # observation rows or original WCVP plant tables.
    report["resource_matrix"]={
        "file":output.name+".matrix.csv.gz",
        "encoding":"bit i corresponds to resource_matrix.species_order[i]",
        "species_order":panel,
        "all_regions":len(official_codes),
        "native_species_x_regions":count_species_units(resources_official["native"]),
        "current_species_x_regions":count_species_units(resources_official["current"]),
    }
    # Rewrite receipt with source-locked species order added.
    output.with_suffix(".json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with gzip.GzipFile(filename="",mode="wb",fileobj=output.with_suffix(".matrix.csv.gz").open("wb"),mtime=0,compresslevel=7) as f:
        wrap=io.TextIOWrapper(f,encoding="utf8",newline="")
        w=csv.writer(wrap)
        w.writerow(["wgsrpd3_code","native_bits_hex","added_bits_hex","current_bits_hex"])
        for i,code in enumerate(official_codes):
            n=resources_official["native"][i]
            c=resources_official["current"][i]
            if n & ~c:
                raise RuntimeError("Native-resource set outside current set for "+code)
            w.writerow([code,format(n,"x"),format(c & ~n,"x"),format(c,"x")])
        wrap.flush();wrap.detach()
    with output.with_suffix(".region.csv").open("w",encoding="utf8",newline="") as f:
        w=csv.writer(f)
        w.writerow(["wgsrpd3_code","plant_native_species","plant_current_species",
                    "butterfly_native_resource_species","butterfly_current_resource_species"])
        for j,code in enumerate(official_codes):
            w.writerow([code,flora_official["native"][j].bit_count(),
                        flora_official["current"][j].bit_count(),
                        resources_official["native"][j].bit_count(),
                        resources_official["current"][j].bit_count()])
    print(json.dumps({
        "units":report["analysis_units"],"plants":plant,
        "consumer_resources":consumer,"delta_spearman":rho,
        "structural_null":null_summary,
    },indent=2))


if __name__=="__main__":
    main()

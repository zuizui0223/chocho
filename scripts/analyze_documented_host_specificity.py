#!/usr/bin/env python3
"""Do *documented larval hosts* explain consumer-resource beta beyond other aliens?

SOURCE-LOWER-BOUND, not biological nonhost identification: casefold exact
binomials in pinned HOSTS source and GloNAF alien-flora source, for the frozen
239 butterfly panel. Unmatched/unlisted taxa are 'not documented as hosts for
this panel', NEVER confirmed nonhosts.

Same 284 WGSRPD3 regions from matched colonial history and butterfly resource
panel, original source coordinates/polygons, full colonial/climate/native
floristic controls, and GloNAF uncatalogued flora background.
The 199 identity-based null samples alternative non-focal-host plants with the
same source geographic prevalence classes as the observed documented hosts.
Nothing here dates plant movements or confirms butterfly occupancy.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import random
import zipfile

import numpy as np

from analyze_colonial_legacy_resource_beta import (
    load as load_colonial_butterfly,orthogonal_controls,
    residual,correlation,
)
from analyze_matched_glonaf_consumer_history import (
    produce_matched_regions,
)
from build_lenzner_tdwg_colonial_crosswalk import load_csv_member

ROOT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
LOCK=ROOT/"documented_host_specificity_lock_v1.json"
GLO_NAF=Path("_temp_historical_sources/lenzner_6640264.zip")
TDWG=Path("_external/wgsrpd/geojson/level3.geojson")
HOSTS=Path("_external/HOSTS/resource.csv")
PANEL=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
OUT=ROOT/"documented_host_specificity_v1.json"
BINS=(1,2,5,10,25,50,100,200,285)


def binomial(genus,epithet):
    g=(genus or "").strip()
    e=(epithet or "").strip()
    if not g or not e or "\ufffd" in g+e:
        return None
    v=(e if e.casefold().startswith(g.casefold()+" ") else g+" "+e)
    if len(v)<5 or " sp." in v.casefold():
        return None
    return v.casefold()


def focal_hosts(raw_hosts=HOSTS,panel_path=PANEL):
    with panel_path.open("r",encoding="utf8",newline="") as f:
        panel={r["species"].strip() for r in csv.DictReader(f)}
    if len(panel)!=239:
        raise RuntimeError("Frozen 239 consumer panel changed")
    names=set()
    per_species=defaultdict(set)
    summary=Counter()
    with raw_hosts.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        necessary={"Insect Genus","Insect Species","Hostplant Genus","Hostplant Species"}
        if not necessary.issubset(reader.fieldnames or []):
            raise RuntimeError("HOSTS source schema changed")
        for r in reader:
            species=(r["Insect Genus"].strip()+" "+r["Insect Species"].strip())
            if species not in panel:
                continue
            summary["focal_insect_host_rows"]+=1
            name=binomial(r["Hostplant Genus"],r["Hostplant Species"])
            if name:
                names.add(name)
                per_species[species].add(name)
                summary["focal_rows_with_identifiable_host_binomial"]+=1
    summary.update({"focal_host_names":len(names),
                    "focal_butterflies_with_named_hosts":len(per_species)})
    if len(names)<500:
        raise RuntimeError("Too few focal HOSTS names for source join")
    return names,dict(summary)


def glonaf_index(sourcezip=GLO_NAF):
    with zipfile.ZipFile(sourcezip) as z:
        rows,damaged=load_csv_member(z,"GloNAF.csv")
    names=set()
    for row in rows:
        v=binomial(row.get("genus"),row.get("Species"))
        if v: names.add(v)
    if len(rows)!=470962:
        raise RuntimeError("GloNAF file no longer source-locked")
    return {name:i for i,name in enumerate(sorted(names))},{
        "total_source_rows":len(rows),"distinct_glonaf_names":len(names),
        "source_decode_replacements":damaged,
    }


def prevalence_bins(full_masks,n_taxa,bin_edges=BINS):
    freq=np.zeros(n_taxa,np.int32)
    for mask in full_masks:
        while mask:
            bit=mask & -mask
            freq[bit.bit_length()-1]+=1
            mask^=bit
    bins=defaultdict(list)
    for i,f in enumerate(freq):
        if f:
            bins[bisect_left(bin_edges,int(f))].append(i)
    return freq,bins


def pair_jaccard(masks,pair_i,pair_j,undefined_zero=False):
    values=np.empty(len(pair_i),float)
    unsupported=0
    for k,(i,j) in enumerate(zip(pair_i,pair_j)):
        a,b=masks[int(i)],masks[int(j)]
        u=(a|b).bit_count()
        if u:
            values[k]=(a&b).bit_count()/u
        else:
            unsupported+=1
            values[k]=0.0 if undefined_zero else float("nan")
    return values,unsupported


def host_freq_matched_null(mask,full_masks,bins,freq,seed=31031,draws=199):
    """Independent pseudohost plant sets matched by observed TDWG prevalence bin.

    Source host 'degree'/geographic prevalence is matched in *bins*, not exact
    count or family. All original region-level plant richness and specific
    native butterfly-resource geography stay fixed. The pseudohost pool omits
    documented focal HOSTS matches.
    """
    if draws<19:
        raise ValueError("Minimum 19 null draws")
    freq_profile=Counter()
    k=mask
    while k:
        bit=k & -k
        idx=bit.bit_length()-1
        if freq[idx]>0:
            freq_profile[bisect_left(BINS,int(freq[idx]))]+=1
        k^=bit
    all_ids={index for group in bins.values() for index in group}
    true_ids=set()
    a=mask
    while a:
        bit=a & -a;true_ids.add(bit.bit_length()-1);a^=bit
    pools={group:[k for k in candidates if k not in true_ids]
           for group,candidates in bins.items()}
    deficits={group:count-len(pools.get(group,[]))
              for group,count in freq_profile.items()
              if len(pools.get(group,[]))<count}
    if deficits:
        raise RuntimeError("Cannot prevalence-match documented hosts without reuse: "+repr(deficits))
    rng=random.Random(seed)
    generated=[]
    for _ in range(draws):
        selected=[]
        for b,count in sorted(freq_profile.items()):
            selected.extend(rng.sample(pools[b],count))
        drawmask=0
        for i in selected:
            drawmask |= 1<<i
        generated.append(drawmask)
    return generated,dict(freq_profile),len(all_ids),len(true_ids)


def evaluate(draws=199,seed=31031):
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock.get("schema")!="chocho_same_284_GloNAF_documented_HOSTS_lower_bound_lock_v1":
        raise RuntimeError("Specificity lock changed")
    data=load_colonial_butterfly()
    assert len(data["ids"])==284
    source_masks,_,source_audit=produce_matched_regions(GLO_NAF,TDWG,(0.05,))
    full_masks=[source_masks[.05][code] for code in data["ids"]]
    names,name_counts=focal_hosts()
    idx,glonaf_counts=glonaf_index()
    if len(idx)!=source_audit["glonaf"]["distinct_species"]:
        raise RuntimeError("GloNAF mask taxon index does not match source parser")
    known=sorted(names.intersection(idx))
    original_unknown=sorted(names-set(idx))
    true_mask=sum(1<<idx[name] for name in known)
    f,bins=prevalence_bins(full_masks,len(idx))
    masks_host=[mask & true_mask for mask in full_masks]
    masks_other=[mask & ~true_mask for mask in full_masks]
    source_active=sum(bool(f[idx[name]]) for name in known)
    if source_active<20:
        raise RuntimeError("Insufficient focal HOSTS × GloNAF active exact taxa: "+str(source_active))
    all_ia,all_ib=data["ia"],data["ib"]
    valid=np.fromiter((bool(masks_host[int(i)]|masks_host[int(j)])
                      for i,j in zip(all_ia,all_ib)),bool,count=len(all_ia))
    if valid.sum()<20000:
        raise RuntimeError("Insufficient shared region-pair support after known-host classification")
    ia,ib=all_ia[valid],all_ib[valid]
    host_jac,host_undefined=pair_jaccard(masks_host,ia,ib)
    other_jac,other_undefined=pair_jaccard(masks_other,ia,ib,undefined_zero=True)
    if host_undefined or not np.all(np.isfinite(host_jac)):
        raise RuntimeError("Fixed observed host support was not nonempty")
    alien_counts=[x.bit_count() for x in full_masks]
    count_min=np.log1p(np.minimum(
        np.asarray(alien_counts)[ia],np.asarray(alien_counts)[ib]))
    count_max=np.log1p(np.maximum(
        np.asarray(alien_counts)[ia],np.asarray(alien_counts)[ib]))
    base=[x[valid] for x in
          data["control_sets"]["colonial_vs_resource_beta_after_flora"]]
    colonial=data["exposure"][valid]
    response=data["resource_delta"][valid]
    common_basis=orthogonal_controls(base+[colonial,count_min,count_max,other_jac])
    yy=residual(response,common_basis)
    hx=residual(host_jac,common_basis)
    observed=correlation(hx,yy)
    if observed is None:
        raise RuntimeError("Known host Jaccard not identifiable")
    other_basis=orthogonal_controls(base+[colonial,count_min,count_max,host_jac])
    n_correlation=correlation(residual(other_jac,other_basis),
                              residual(response,other_basis))

    fake_masks,prevalence_profile,n_candidates,n_known=host_freq_matched_null(
        true_mask,full_masks,bins,f,seed=seed,draws=draws)
    null=[]
    absent_fakes=[]
    for drawmask in fake_masks:
        ms=[mask & drawmask for mask in full_masks]
        x,missing=pair_jaccard(ms,ia,ib,undefined_zero=True)
        v=correlation(residual(x,common_basis),yy)
        if v is None:
            raise RuntimeError("Prevalence-matched fake-host predictor constant")
        null.append(v)
        absent_fakes.append(missing)
    return {
        "schema":"chocho_fixed239_glonaf_documented_larval_host_specificity_v1",
        "status":"SOURCE_MATCHED_EXACT_NAME_HOST_CANDIDATE_TEST_NOT_CAUSAL_OR_TRUE_NONHOST",
        "source_pins":{
            "original_host_file":"globalbioticinteractions/HOSTS 808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "source_glonaf_zip_md5":"16f5585d2500676951a4b664bbf0ffde",
            "lookup_region_count":len(data["ids"]),
            "method_lock_sha256":hashlib.sha256(LOCK.read_bytes()).hexdigest(),
        },
        "source_alignment":{
            "host_input":name_counts,"glonaf":glonaf_counts,
            "original_GloNAF_mapped_audit":source_audit,
            "focal_HOSTS_host_names_exact_GloNAF":len(known),
            "source_ACTIVE_focal_host_names":source_active,
            "focal_host_names_not_matched_GloNAF":len(original_unknown),
            "source_same_region_pairs":len(all_ia),
            "eligible_known_host_region_pairs":int(valid.sum()),
            "excluded_both_no_documented_host":int((~valid).sum()),
            "mean_documented_host_Jaccard":float(np.mean(host_jac)),
            "mean_not_documented_plant_Jaccard":float(np.mean(other_jac)),
            "other_plant_both_empty_zero_handled":int(other_undefined),
        },
        "adjusted_partial_rank_associations":{
            "documented_focal_HOSTS_host_Jaccard_given_uncatalogued_aliens":observed,
            "uncatalogued_GloNAF_alien_Jaccard_given_documented_hosts":n_correlation,
        },
        "geographical_prevalence_matched_null":{
            "draws":draws,"seed":seed,
            "host_source_region_frequency_bins":prevalence_profile,
            "eligible_GloNAF_taxa":n_candidates,
            "excluded_true_hosts_from_pseudohost_draws":n_known,
            "null_median":float(np.median(null)),
            "null_q025":float(np.quantile(null,.025)),
            "null_q975":float(np.quantile(null,.975)),
            "observed":observed,
            "upper_tail":(1+sum(q>=observed for q in null))/(draws+1),
            "two_sided_abs_tail":(1+sum(abs(q)>=abs(observed) for q in null))/(draws+1),
            "mean_fake_undefined_pairs_zero_handled":float(np.mean(absent_fakes)),
            "max_fake_undefined_pairs_zero_handled":int(max(absent_fakes)),
        },
        "interpretation_limits":[
            "Only exact source BINOMIAL matches; nonmatch means uncatalogued here, never proven nonhost.",
            "HOSTS and GloNAF taxonomy resolution differences may conceal genuinely compatible hosts.",
            "The source host-list was collected with heterogeneous study intensity and later sampling.",
            "Prevalence bins match the number of GloNAF regions where taxa are recorded, not precise phylogenetic plant identity or family.",
            "The reference treats both-empty fake-host pair similarity as zero; observed documented-host comparisons require nonempty union.",
            "Same 284 regions reuse pairs and do not identify true butterfly occupancy or larval utilization.",
            "No historical dates proving host plants preceded butterflies, nor proof of port/colonial transport.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--draws",type=int,default=199)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    result=evaluate(args.draws)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in result.items()
                      if k!="interpretation_limits"},indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

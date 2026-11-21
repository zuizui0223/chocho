#!/usr/bin/env python3
"""Source-locked Figure 3B documentation-coverage sensitivity; NOT an effort-causal model.

No external dependencies. This audit uses the original 239-species frozen
HOSTS/WCVP panel and the 126-island by-species resource table. The number of
distinct resolved HOSTS plants is ONLY a documentation-richness proxy, not the
number of independent publications, actual recording effort, or host quality.
Both outcome and predictor are computed from linked HOSTS/WCVP inputs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

SOURCE_BLOBS = {
    "metrics": "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
    "island": "a230f6d214c96f874a3ab184ee9823cc9023c1bb",
}
N = 239
SEED = 20261121

def blob_sha(path):
    data=Path(path).read_bytes()
    return hashlib.sha1(("blob "+str(len(data))).encode()+bytes([0])+data).hexdigest()

def readcsv(path):
    with open(path,newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def ranks(v):
    order=sorted(range(len(v)),key=lambda i:v[i])
    out=[0.]*len(v)
    i=0
    while i<len(v):
        j=i+1
        while j<len(v) and v[order[j]]==v[order[i]]:
            j+=1
        val=(i+1+j)/2
        for k in range(i,j):out[order[k]]=val
        i=j
    return out

def corr(x,y):
    n=len(x);mx=sum(x)/n;my=sum(y)/n
    a=[v-mx for v in x];b=[v-my for v in y]
    den=(sum(v*v for v in a)*sum(v*v for v in b))**0.5
    return sum(v*w for v,w in zip(a,b))/den if den>0 else float("nan")

def spearman(x,y):
    return corr(ranks(x),ranks(y))

def residual(y,X):
    k=len(X[0]);A=[[0.]*(k+1) for _ in range(k)]
    for yy,row in zip(y,X):
        for i in range(k):
            A[i][k]+=row[i]*yy
            for j in range(k):
                A[i][j]+=row[i]*row[j]
    for i in range(k):
        p=max(range(i,k),key=lambda j:abs(A[j][i]))
        A[i],A[p]=A[p],A[i]
        if abs(A[i][i])<1e-9:return [float("nan")]*len(y)
        d=A[i][i]
        for j in range(i,k+1):A[i][j]/=d
        for r in range(k):
            if r==i:continue
            f=A[r][i]
            for j in range(i,k+1):A[r][j]-=f*A[i][j]
    coef=[row[k] for row in A]
    return [v-sum(b*a for b,a in zip(coef,row)) for v,row in zip(y,X)]

def partial_rank(rows,controls=("host_count","family_count"),outcome="added"):
    x=ranks([r["native_global"] for r in rows])
    y=ranks([r[outcome] for r in rows])
    c=[ranks([r[name] for r in rows]) for name in controls]
    X=[[1.]+[z[i] for z in c] for i in range(len(rows))]
    return corr(residual(x,X),residual(y,X))

def tier(n):
    if n<=1:return "1"
    if n==2:return "2"
    if n<=4:return "3-4"
    if n<=9:return "5-9"
    return "10+"

def percentiles(values):
    x=sorted(z for z in values if math.isfinite(z))
    if not x:return [None,None]
    def quantile(p):
        a=(len(x)-1)*p;lo=int(a);hi=min(lo+1,len(x)-1)
        return x[lo]*(1-(a-lo))+x[hi]*(a-lo)
    return [round(quantile(.025),6),round(quantile(.975),6)]

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--metrics",default="data/frozen/figure_sources/anthropogenic_species_metrics.csv")
    p.add_argument("--islands",default="results/island_resource_gap_filling/by_species.csv")
    p.add_argument("--out-dir",default="results/geb_v05_fig3b_documentation_proxy_20261121")
    p.add_argument("--resamples",type=int,default=2000)
    a=p.parse_args()
    assert a.resamples>=199
    inputs={"metrics":a.metrics,"island":a.islands}
    shas={k:blob_sha(path) for k,path in inputs.items()}
    assert shas==SOURCE_BLOBS,("Unexpected changed source inputs",shas)
    m=readcsv(a.metrics); b=readcsv(a.islands)
    assert len(m)==N and len({r["species"] for r in m})==N
    bi={r["species"]:r for r in b if r["group"]=="island"}
    assert len(bi)==N and all(int(v["cells"])==126 for v in bi.values())
    assert {r["species"] for r in m}==set(bi)
    rows=[]
    for r in m:
        name=r["species"];i=bi[name];n=int(r["resolved_host_species"])
        rows.append(dict(
            species=name,
            native_global=int(r["native_resource_units"]),
            added=int(i["added"]),
            native_island=int(i["native"]),
            host_count=n,
            family_count=float(r["host_family_count"]),
            adequate=r["host_taxonomy_lower_bound_adequate"]=="True",
            host_tier=tier(n),
        ))
    assert sum(r["added"] for r in rows)==5741
    original=spearman([r["native_global"] for r in rows],[r["added"] for r in rows])
    assert abs(original-0.5788964068803357)<1e-10
    adjusted=partial_rank(rows)
    adjusted_hostonly=partial_rank(rows,("host_count",))
    adequate=[r for r in rows if r["adequate"]]
    onefamily=[r for r in rows if r["family_count"]==1]
    positive=[r for r in rows if r["added"]>0]
    rng=random.Random(SEED)
    bs=[]
    for _ in range(a.resamples):
        sample=[rows[rng.randrange(N)] for _ in range(N)]
        bs.append(partial_rank(sample))
    leaveone=[partial_rank(rows[:k]+rows[k+1:]) for k in range(N)]
    tiers=defaultdict(list)
    for j,r in enumerate(rows):tiers[r["host_tier"]].append(j)
    perm=[]
    for _ in range(a.resamples):
        y=[r["added"] for r in rows]
        for ids in tiers.values():
            vals=[y[j] for j in ids]
            rng.shuffle(vals)
            for j,v in zip(ids,vals):y[j]=v
        new=[dict(row,added=y[j]) for j,row in enumerate(rows)]
        perm.append(partial_rank(new))
    p_two=(1+sum(abs(t)>=abs(adjusted) for t in perm if math.isfinite(t)))/(1+sum(math.isfinite(t) for t in perm))
    gain_gap=[dict(r,frac_added=r["added"]/(126-r["native_island"])) for r in rows]
    result={
        "schema":"chocho_geb_v05_fig3B_source_239_species_documentation_proxy_sensitivity_v1",
        "status":"SOURCE_239_EXACT;EXPLORATORY_DOCUMENTATION_PROXY_NOT_INDEPENDENT_STUDY_EFFORT",
        "provenance":{"source_git_blob_sha":shas,"seed":SEED,"resamples":a.resamples},
        "n_species":N,"island_regions":126,
        "original_spearman_native_global_vs_added_island":round(original,9),
        "spearman_recorded_host_species_vs_native_global":round(spearman([r["host_count"] for r in rows],[r["native_global"] for r in rows]),6),
        "spearman_recorded_host_species_vs_added_island":round(spearman([r["host_count"] for r in rows],[r["added"] for r in rows]),6),
        "partial_spearman_controlling_documented_host_species_only":round(adjusted_hostonly,6),
        "partial_spearman_controlling_documented_host_species_and_families":round(adjusted,6),
        "species_bootstrap_percentile95_exploratory":percentiles(bs),
        "within_host_documentation_tier_permutation_two_sided_p_exploratory":round(p_two,6),
        "leave_one_species_out_partial_rho_range":[round(min(leaveone),6),round(max(leaveone),6)],
        "host_count_tier_sizes":dict(sorted(Counter(r["host_tier"] for r in rows).items())),
        "adequate_taxonomy_subset":{"n":len(adequate),"partial_rho":round(partial_rank(adequate),6)},
        "one_documented_family_subset":{"n":len(onefamily),"partial_rho_host_count_only":round(partial_rank(onefamily,("host_count",)),6)},
        "positive_island_gain_subset":{"n":len(positive),"partial_rho":round(partial_rank(positive),6)},
        "gap_fill_fraction_sensitivity":{"partial_rho_controls_host_count_and_family_count":round(partial_rank(gain_gap,outcome="frac_added"),6)},
        "limitations":[
          "Distinct resolved known-host plant species is a documentation richness proxy, not the independent number of HOSTS reports or publications.",
          "Both global native resource breadth and island added units derive from the same HOSTS/WCVP network; residual association does not prove independent ecological advantage.",
          "Species bootstrap ignores phylogeny and source selection, and tier-wise permutations assume within-tier exchangeability. Both are descriptive stress tests, not confirmatory population p values.",
          "All source cohorts, taxa and original positive/adverse conclusions remain unchanged; no inferred butterfly establishment or host fitness."
        ]
    }
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/"FIG3B_DOCUMENTATION_PROXY_SENSITIVITY_20261121.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    report=f"""# Figure 3B: documentation proxy sensitivity (239 source butterflies)

**Scope: source-exact, preliminary sensitivity; not an observed field research-effort correction.**

Original Spearman native global resource range vs absolute added island resources: **{original:+.3f}**.
Resolved documented host-plant species vs island gain: **{result['spearman_recorded_host_species_vs_added_island']:+.3f}**.
After partial-ranking out documented-host plant species and host-family count: **{adjusted:+.3f}**.
Exploratory resampled 95% interval **{result['species_bootstrap_percentile95_exploratory']}**; tiered permutation p **{p_two:.4f}**.
All 239 single-species removals **{result['leave_one_species_out_partial_rho_range']}**.
More strictly taxonomically adequate subset n={len(adequate)} partial ρ={result['adequate_taxonomy_subset']['partial_rho']:+.3f}.
Exactly one documented host-plant family subset n={len(onefamily)} partial ρ={result['one_documented_family_subset']['partial_rho_host_count_only']:+.3f}.

**Interpretation:** The original positive Figure 3B association is not erased by adjusting for the number
of distinct accepted HOSTS host plants and host families. However, this proxy does not capture
independent sampling effort, bibliography size, temporal source coverage or phylogenetic
dependence, and HOSTS/WCVP is the common data generator of both axes. Do NOT call this
proof of a causal or independent ecological mechanism. Preserve the original ρ=+0.579
and give this diagnostic in SI, not as a new main result.

Source blobs: \`{shas}\`. Deterministic seed: {SEED}; resamples: {a.resamples}.
"""
    (out/"README_JA.md").write_text(report,encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":main()

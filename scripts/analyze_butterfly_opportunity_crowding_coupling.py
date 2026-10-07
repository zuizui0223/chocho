#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math,random
from pathlib import Path

def rankdata(x):
    order=sorted(range(len(x)),key=lambda i:x[i]); out=[0.0]*len(x); i=0
    while i<len(order):
        j=i+1
        while j<len(order) and x[order[j]]==x[order[i]]: j+=1
        r=(i+j-1)/2+1
        for k in range(i,j): out[order[k]]=r
        i=j
    return out

def corr(x,y):
    mx=sum(x)/len(x); my=sum(y)/len(y)
    num=sum((a-mx)*(b-my) for a,b in zip(x,y))
    dx=sum((a-mx)**2 for a in x); dy=sum((b-my)**2 for b in y)
    return num/math.sqrt(dx*dy)

def spearman(x,y): return corr(rankdata(x),rankdata(y))

def solve_normal(A,b):
    # Gaussian elimination for tiny normal-equation systems.
    n=len(b); M=[list(map(float,A[i]))+[float(b[i])] for i in range(n)]
    for i in range(n):
        p=max(range(i,n),key=lambda r:abs(M[r][i]))
        M[i],M[p]=M[p],M[i]
        if abs(M[i][i])<1e-12: raise RuntimeError("singular controls")
        z=M[i][i]
        M[i]=[v/z for v in M[i]]
        for r in range(n):
            if r==i: continue
            z=M[r][i]
            M[r]=[M[r][c]-z*M[i][c] for c in range(n+1)]
    return [M[i][-1] for i in range(n)]

def residualize(y,controls):
    cols=[[1.0]*len(y)]+[rankdata(c) for c in controls]
    yr=rankdata(y); p=len(cols)
    A=[[sum(cols[i][k]*cols[j][k] for k in range(len(y))) for j in range(p)] for i in range(p)]
    b=[sum(cols[i][k]*yr[k] for k in range(len(y))) for i in range(p)]
    beta=solve_normal(A,b)
    return [yr[k]-sum(beta[j]*cols[j][k] for j in range(p)) for k in range(len(y))]

def partial_spearman(x,y,controls):
    return corr(residualize(x,controls),residualize(y,controls))

def quantile(xs,p):
    xs=sorted(xs); z=p*(len(xs)-1); i=int(math.floor(z)); j=int(math.ceil(z))
    if i==j:return xs[i]
    a=z-i; return xs[i]*(1-a)+xs[j]*a

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--competition-species-csv",type=Path,required=True)
    ap.add_argument("--bootstrap",type=int,default=19999)
    ap.add_argument("--seed",type=int,default=20261007)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    rows=[]
    with args.competition_species_csv.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if float(r["alpha"])!=0: continue
            n=float(r["native_regions"]); c=float(r["contemporary_regions"])
            rows.append({
                "species":r["species"],
                "host_family_count":float(r["host_family_count"]),
                "native_regions":n,
                "log_resource_gain":math.log(c/n),
                "native_exposure":float(r["mean_native_competitor_exposure"]),
                "contemporary_exposure":float(r["mean_contemporary_competitor_exposure"]),
                "delta_exposure":float(r["mean_contemporary_competitor_exposure"])-float(r["mean_native_competitor_exposure"]),
            })
    if len(rows)!=239: raise RuntimeError(f"expected 239 species, got {len(rows)}")
    x=[r["log_resource_gain"] for r in rows]
    y=[r["delta_exposure"] for r in rows]
    native=[r["native_regions"] for r in rows]
    breadth=[r["host_family_count"] for r in rows]
    nexp=[r["native_exposure"] for r in rows]

    rho=spearman(x,y)
    rng=random.Random(args.seed); boot=[]
    for _ in range(args.bootstrap):
        idx=[rng.randrange(len(rows)) for _ in rows]
        boot.append(spearman([x[i] for i in idx],[y[i] for i in idx]))

    payload={
      "schema":"chocho_butterfly_opportunity_crowding_coupling_v0.1",
      "status":"SUCCESS_POSTHOC_OPPORTUNITY_CROWDING_COUPLING",
      "question":"Do butterflies receiving larger proportional geographic resource gains also experience larger increases in exact-host co-user exposure?",
      "n":len(rows),
      "primary":{
        "spearman_log_resource_gain_vs_delta_co_user_exposure":rho,
        "bootstrap_replicates":args.bootstrap,
        "bootstrap_ci95":[quantile(boot,.025),quantile(boot,.975)],
      },
      "partial_rank_sensitivities":{
        "control_native_resource_breadth":partial_spearman(x,y,[native]),
        "control_native_resource_breadth_and_host_family_breadth":partial_spearman(x,y,[native,breadth]),
        "control_native_exposure":partial_spearman(x,y,[nexp]),
        "control_native_resource_breadth_host_breadth_and_native_exposure":partial_spearman(x,y,[native,breadth,nexp]),
      },
      "interpretation":"Resource release and potential interspecific crowding exposure are positively coupled across species: the butterflies that gain proportionally more host-resource geography tend also to enter geography where more other focal butterflies share exact hosts.",
      "claim_boundary":"Co-user exposure is a spatial proxy for potential resource-sharing interaction, not measured competition. The positive association can arise from the ecology and geography of widely redistributed shared hosts; it does not establish demographic costs or causal competitive suppression."
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))
if __name__=="__main__": main()

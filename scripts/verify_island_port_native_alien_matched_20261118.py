#!/usr/bin/env python3
"""Reproduce *paired* WCVP native/introduced whole-flora port regressions.

Original, immutable 126-island WCVP/WPI regional table, NOT a sample of 126
independent geographical natural experiments. No butterfly, AIS, or vessel
transport data enter these regressions. Python standard library only.

One identically specified design matrix:
  log1p(native or introduced accepted species richness)
     ~ log2(1+WPI2023 port count) + log(region area) + TDWG Level-1 dummies
A paired contrast regresses log1p(introduced)-log1p(native) on the same design;
its HC3 interval accounts for the correlated two responses. Source-selected
123-island resource-active cohort is a separate, declared sensitivity.
HC3 describes island-specific heteroskedasticity, NOT phylogenetic or spatial
independence, historical shipping, or a zero native effect.
"""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGIONS = ROOT / "results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv"
ACTIVE = ROOT / "results/three-gate-port-introduced-host-taxa/distinct_introduced_host_taxa_and_port_by_island.csv"

def require(ok, message):
    if not ok:
        raise AssertionError(message)

def transpose(a):
    return list(zip(*a))

def invert(matrix):
    n = len(matrix)
    augmented = [list(row)+[float(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for i in range(n):
        pivot = max(range(i,n), key=lambda k:abs(augmented[k][i]))
        require(abs(augmented[pivot][i])>1e-11,"Singular design")
        augmented[i],augmented[pivot]=augmented[pivot],augmented[i]
        fac=augmented[i][i]
        augmented[i]=[z/fac for z in augmented[i]]
        for j in range(n):
            if j==i: continue
            fac=augmented[j][i]
            augmented[j]=[z-fac*w for z,w in zip(augmented[j],augmented[i])]
    return [r[n:] for r in augmented]

def fit(rows, response):
    levels=sorted({int(r["level1"]) for r in rows})
    require(levels==list(range(1,10)),"TDWG Level-1 source labels changed")
    x=[[1.,float(r["log2ports"]),float(r["log_area"])]+
       [float(int(r["level1"])==lv) for lv in levels[1:]] for r in rows]
    if response=="difference":
        y=[math.log1p(float(r["introduced"]))-math.log1p(float(r["native"])) for r in rows]
    else:
        y=[math.log1p(float(r[response])) for r in rows]
    n,p=len(x),len(x[0])
    xtx=[[sum(x[i][j]*x[i][k] for i in range(n)) for k in range(p)] for j in range(p)]
    xty=[sum(x[i][j]*y[i] for i in range(n)) for j in range(p)]
    inv=invert(xtx)
    coef=[sum(z*w for z,w in zip(inv[j],xty)) for j in range(p)]
    vari=0.0
    for xi,yi in zip(x,y):
        residual=yi-sum(a*b for a,b in zip(xi,coef))
        lever=sum(xi[j]*inv[j][k]*xi[k] for j in range(p) for k in range(p))
        require(0<=lever<1,"Hat leverage out of range")
        # sandwich: [ (X'X)^-1 X' diag(e_i^2/(1-h_i)^2) X (X'X)^-1 ]_11
        influence=sum(inv[1][j]*xi[j] for j in range(p))
        vari+=(residual/(1.-lever)*influence)**2
    se=math.sqrt(vari)
    return {"n":n,"covariates":p,"beta":coef[1],
            "HC3_se":se,"HC3_95":[coef[1]-1.96*se,coef[1]+1.96*se]}

def main():
    with REGIONS.open(encoding="utf-8-sig",newline="") as file:
        rows=list(csv.DictReader(file))
    with ACTIVE.open(encoding="utf-8-sig",newline="") as file:
        active={row["code"] for row in csv.DictReader(file)}
    require(len(rows)==126 and len({r["code"] for r in rows})==126,
            "126 source-defined botanical islands")
    require(len(active)==123,"123 source-defined resource-active islands")
    require(sum(int(r["native"]) for r in rows)==184997 and
            sum(int(r["introduced"]) for r in rows)==39024,
            "Original WCVP native/alien richness sums")
    require(all(abs(float(r["log2ports"])-math.log2(1+int(r["wpi_port_count"])))<1e-9
                for r in rows),"Original modern port exposure")
    report={"schema":"island_port_native_alien_matched_20261118",
            "design":"exact same exposure, area, 8 geographic fixed-effect dummies",
            "estimand":"observational log1p(species richness) slope per log2(1+ports)",
            "model_scopes":{}}
    for key,pop in (("all_126",rows),("resource_active_123",[r for r in rows if r["code"] in active])):
        outcome={k:fit(pop,k) for k in ("native","introduced","difference")}
        require(abs(outcome["introduced"]["beta"]-outcome["native"]["beta"]-
                    outcome["difference"]["beta"])<1e-10,"Paired response identity")
        report["model_scopes"][key]=outcome
    a=report["model_scopes"]["all_126"]
    for key,want in (("native",.0954313525414),("introduced",.28544306145485),
                     ("difference",.1900117089134)):
        require(abs(a[key]["beta"]-want)<1e-8,"Fixed original 126 regression: "+key)
    require(a["native"]["HC3_95"][0]<0<a["native"]["HC3_95"][1] and
            a["introduced"]["HC3_95"][0]>0 and
            a["difference"]["HC3_95"][0]<0<a["difference"]["HC3_95"][1],
            "Qualifying uncertainty: native and contrast HC3 intervals include zero")
    report["limitations"]=[
        "No direct historical shipping movements, provenance or causal plant import",
        "Cannot infer a statistical zero effect for native flora",
        "Paired alien-minus-native coefficient difference is uncertain",
        "Island HC3 intervals do not account for full spatial correlation",
        "WCVP regional present-day floras and WPI current ports are cross-sectional",
        "This is not the 126 native-host-resource-adjusted alien beta+0.221 nor 98 botanical-effort alien beta+0.2445"
    ]
    print(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True))

if __name__=="__main__":
    main()

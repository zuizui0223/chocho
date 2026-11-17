#!/usr/bin/env python3
"""Source-locked 51-island nested uniform-ocean circuit vs introduced candidate hosts.

Uses EXACT original GitHub Actions archive 42193019385 (artifact 13490383689)
and the outcome-independent longitude-block assignment from original run
42325176725 (artifact 13554313981). NO species re-index, shore remapping,
new butterfly data or post-hoc island exclusion. B, B+U, B+H, B+U+H are fitted
to the same later and earlier butterfly community-Jaccard samples.

The results describe predictive overlap, NOT biological plant mediation,
confirmed insect settlement, directional ocean currents or vessel transport.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

PIN_PAIRS="133dc6d6a91203d09385a4ffd4db0b712dc3309130d1bf2a09bcee6bafa36bf0"
PIN_PREP="5641711fc7ab2bdcd4626e1f483faf46c945a5d11466f0deada092db456a0363"
PIN_BLOCK="2e4b8b1d4bda18f7193f843378da124b597ead0f6c49751a547cab10f65e2ef3"
MODELS=("B","BU","BH","BUH")
CONTRASTS=(("B","BU"),("B","BH"),("B","BUH"),("BH","BUH"),("BU","BUH"))

def verify_sha(path,expected):
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest!=expected:raise ValueError(f"Original source SHA256 mismatch {path}: {digest}")
    return digest

def predict(x,y,train,test):
    z=x[train]; v=x[test];mu=z.mean(axis=0);sd=z.std(axis=0)
    sd=np.where(sd<1.e-12,1.,sd)
    X=np.column_stack((np.ones(len(z)),(z-mu)/sd))
    V=np.column_stack((np.ones(len(v)),(v-mu)/sd))
    beta=np.linalg.lstsq(X,y[train],rcond=None)[0]
    out=np.clip(V@beta,0.,1.)
    if not np.isfinite(out).all():raise ValueError("Nonfinite prediction")
    return out

def scores(tab,target,design,folds,method):
    y=tab[target].to_numpy(dtype=float)
    rows=[]; preds=[]
    used=np.zeros(len(tab),int)
    for fold_id,train,test in folds:
        if np.any(train&test):raise ValueError("Leakage between train and test")
        if len(set(tab.loc[train,"island_a"])&set(tab.loc[test,"island_a"]) | set(tab.loc[train,"island_b"])&set(tab.loc[test,"island_b"])) and method=="double_unseen":
            raise ValueError("A held-out geographic island is present in training")
        used+=test.astype(int)
        row={"target":target,"design":method,"fold":fold_id,
             "train_pairs":int(train.sum()),"test_pairs":int(test.sum())}
        for m in MODELS:
            yhat=predict(design[m],y,train,test)
            e=y[test]-yhat
            row["MSE_"+m]=float(np.mean(e**2))
            row["MAE_"+m]=float(np.mean(abs(e)))
            for ind,pred,obs in zip(np.flatnonzero(test),yhat,y[test]):
                preds.append({"target":target,"design":method,"fold":fold_id,
                              "source_pair_index":int(ind),"model":m,
                              "observed":float(obs),"predicted":float(pred),
                              "squared_error":float((obs-pred)**2)})
        rows.append(row)
    if method=="LOIO":
        if len(rows)!=51 or not np.array_equal(used,np.full(len(tab),2,int)):
            raise ValueError("Original 51 island left-out support not reproduced")
        avg={m:float(np.mean([z["MSE_"+m] for z in rows])) for m in MODELS}
    else:
        if len(rows)!=15 or not np.array_equal(used,np.ones(len(tab),int)):
            raise ValueError("Original 15 geographic fold pair support not reproduced")
        pf=pd.DataFrame(preds)
        avg={m:float(pf.loc[pf.model==m,"squared_error"].mean()) for m in MODELS}
    increment={a+"_to_"+b:float(100*(avg[a]-avg[b])/avg[a]) for a,b in CONTRASTS}
    n_improved={a+"_to_"+b:int(sum(z["MSE_"+b]<z["MSE_"+a] for z in rows))
                for a,b in CONTRASTS}
    output={"target":target,"design":method,"n_source_islands":51,
            "n_unique_source_pairs":1275,"n_fold_records":len(rows),
            "MSE":avg,"relative_MSE_reduction_pct":increment,
            "folds_improved":n_improved}
    if method=="LOIO":
        rng=np.random.default_rng(20261117+(1 if target=="new_butterfly_jaccard" else 2))
        mat=np.asarray([[v["MSE_"+m] for m in MODELS] for v in rows])
        boot={k:[] for k in increment}
        for _ in range(4999):
            q=mat[rng.integers(0,len(mat),len(mat))].mean(0)
            for a,b in CONTRASTS:
                k=a+"_to_"+b
                boot[k].append(float(100*(q[MODELS.index(a)]-q[MODELS.index(b)])/q[MODELS.index(a)]))
        output["descriptive_island_bootstrap95"]={k:[float(x) for x in np.quantile(v,[.025,.975])] for k,v in boot.items()}
        # Not a causal mediated fraction; it is retained only as loss overlap.
        delta=avg["B"]-avg["BU"]
        output["uniform_increment_retained_conditional_on_hosts"]=(
             float((avg["BH"]-avg["BUH"])/delta) if abs(delta)>1.e-12 else None)
    return output,rows,preds

def main():
    ap=argparse.ArgumentParser()
    for field in ("pairs","prep","longitude-blocks","output"):
        ap.add_argument("--"+field,required=True,type=Path)
    a=ap.parse_args()
    verify_sha(a.pairs,PIN_PAIRS)
    verify_sha(a.prep,PIN_PREP)
    verify_sha(a.longitude_blocks,PIN_BLOCK)
    df=pd.read_csv(a.pairs,dtype={"island_a":str,"island_b":str})
    prep=json.loads(a.prep.read_text())
    region=pd.read_csv(a.longitude_blocks,dtype={"island_id":str})
    ids=sorted(prep["by_island"],key=int)
    if (len(ids),len(df),len(region),prep["original_ports_audited"],
        prep["marine_ports_included"])!=(51,1275,51,194,182):
        raise ValueError("Original source cohort changed")
    if set(region.island_id)!=set(ids) or sorted(region.longitude_block.value_counts())!=[10,10,10,10,11]:
        raise ValueError("Original source-only geographic block population changed")
    pairs={tuple(sorted((l,r),key=int)) for l,r in zip(df.island_a,df.island_b)}
    required={tuple(sorted((l,r),key=int)) for i,l in enumerate(ids) for r in ids[i+1:]}
    if pairs!=required or not np.allclose(df.port_distance_km,df.distance_km_verified,rtol=1e-9):
        raise ValueError("Original complete unique pair coverage or distances changed")
    counts={k:int(v["original_ports"]) for k,v in prep["by_island"].items()}
    if sum(counts.values())!=194 or any(z<=0 for z in counts.values()):
        raise ValueError("Original 194 ports not recovered")
    d=np.log1p(df.port_distance_km.to_numpy(float))
    pa=np.log1p(df.island_a.map(counts).to_numpy(float))
    pb=np.log1p(df.island_b.map(counts).to_numpy(float))
    u=df.log1p_normalized_uniform_circuit.to_numpy(float)
    h=df.host_jaccard.to_numpy(float)
    label=dict(zip(region.island_id,region.longitude_block))
    A=df.island_a.map(label).to_numpy(int);B=df.island_b.map(label).to_numpy(int)
    folds=[]
    for s in range(5):
        for t in range(s,5):
            group=((A==s)|(A==t)),((B==s)|(B==t))
            train=~group[0]&~group[1]
            test=((A==s)&(B==t))|((A==t)&(B==s))
            folds.append((f"B{s}_B{t}",train,test))
    res={};allfold=[];allpred=[]
    for target in ("new_butterfly_jaccard","earlier_butterfly_jaccard"):
        x=[d,d*d]
        if target=="new_butterfly_jaccard":
            x.append(df.earlier_butterfly_jaccard.to_numpy(float))
        x.extend((pa+pb,abs(pa-pb)))
        base=np.column_stack(x)
        models={"B":base,"BU":np.column_stack((base,u)),
                "BH":np.column_stack((base,h)),
                "BUH":np.column_stack((base,u,h))}
        source_folds=[]
        for id in ids:
            test=(df.island_a.eq(id)|df.island_b.eq(id)).to_numpy()
            train=~test
            if (test.sum(),train.sum())!=(50,1225):
                raise ValueError("Original island LOIO fold support changed")
            source_folds.append((id,train,test))
        res[target]={}
        for mode,loc in (("LOIO",source_folds),("double_unseen",folds)):
            o,fr,pr=scores(df,target,models,loc,mode)
            res[target][mode]=o
            allfold.extend(fr);allpred.extend(pr)
    later=res["new_butterfly_jaccard"]["LOIO"]["relative_MSE_reduction_pct"]
    early=res["earlier_butterfly_jaccard"]["LOIO"]["relative_MSE_reduction_pct"]
    for observed,expected in ((later["B_to_BU"],5.0818868452225594),
                              (later["BU_to_BUH"],2.2759926064420992),
                              (early["B_to_BU"],-0.0318640755121882),
                              (early["BU_to_BUH"],5.3700262016993241)):
        if abs(observed-expected)>1.e-8:
            raise ValueError(f"Prior source analysis did not reproduce: {observed} vs {expected}")
    output=a.output
    output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(allfold).to_csv(output/"nested_uniform_host_fold_losses.csv",index=False)
    pd.DataFrame(allpred).to_csv(output/"nested_uniform_host_predictions.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    receipt={"schema":"chocho_51_source_uniform_ocean_and_candidate_host_joint_LOIO_double_unseen_v1",
             "status":"SOURCE_PINNED_DIAGNOSTIC_NOT_MEDIATION_OR_COLONISATION",
             "runs":{"uniform_source":42193019385,"longitude_source":42325176725},
             "n_islands":51,"n_source_pairs":1275,
             "source_SHA256":{"pairs":PIN_PAIRS,"prep":PIN_PREP,"blocks":PIN_BLOCK},
             "models":{"B":"distance, distance squared, original WPI port sums and differences, pre-target butterfly similarity for later response only","BU":"B + undirected original GSHHG uniform marine circuit","BH":"B + candidate introduced host Jaccard pooled 2001-2024","BUH":"B + uniform circuit + host Jaccard"},
             "results":res,
             "limits":["Pooled candidate introduced hosts are not WCVP whole introduced flora nor proven feeding","2001-2024 pooled host records can postdate 2011-2022 butterflies",
             "Original 51 landmasses are different from 104 TDWG botanical regions and 22 EuDiS trait-covered regions",
             "The LOIO design reuses the other island, double-unseen tests 15 source-only longitude-block pair folds",
             "The two-unseen B/BH/BU/BUH baseline does not include early-plant indexing or marginal candidate-host richness controlled in another 3-model source comparison",
             "Uniform electrical circuit is undirected, not NOAA current or shipborne host transfer",
             "Loss overlap does not estimate causal plant-mediated transport; first GBIF indexing is not established butterfly arrival"] }
    (output/"verified_nested_uniform_host_51island.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("SOURCE_VERIFIED_51_NESTED_UNIFORM_HOST",json.dumps({
      k:{m:res[k][m]["relative_MSE_reduction_pct"] for m in ("LOIO","double_unseen")}
      for k in res},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()

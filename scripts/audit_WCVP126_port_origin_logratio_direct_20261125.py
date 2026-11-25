#!/usr/bin/env python3
"""Direct, same-island introduced/native flora log-ratio port diagnostic."""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t
OUTER="3a9bce3a344068da68e293d66676aa94f038c38c0efad0935c5bc4568368a92c"
INNER="d9af86284a8cec0b82257c096094ad21c91724c8b2ae0abeeffef0d004a8a577"
def sha(b):return hashlib.sha256(b).hexdigest()
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def fit(d,geog,eff):
    y=np.log1p(d.introduced_plant)-np.log1p(d.native_plant)
    X=pd.DataFrame({"const":1.,"ports":np.log2(1+d.wpi_port_count),"area":np.log(d.area_km2)},index=d.index)
    if geog:
        X["distance"]=np.log1p(d.continent_distance_km)
        X["abs_lat"]=d.abs_centroid_latitude
    if eff:X["effort"]=np.log1p(d.gbif_all_plant_indexed_records)
    X=pd.concat([X,pd.get_dummies(d.level1,prefix="L1",drop_first=True,dtype=float)],axis=1)
    need(np.isfinite(X).all().all() and np.linalg.matrix_rank(X)==X.shape[1],"full rank")
    hc=sm.OLS(y,X).fit(cov_type="HC3")
    cl=sm.OLS(y,X).fit(cov_type="cluster",cov_kwds={"groups":d.level1})
    b=float(cl.params["ports"]);se=float(cl.bse["ports"]);df=int(d.level1.nunique()-1)
    return {"n":len(d),"beta":b,"HC3_p":float(hc.pvalues["ports"]),
            "cluster_df":df,"cluster_t_p":float(2*t.sf(abs(b/se),df)),
            "cluster_t_95":[float(b-t.ppf(.975,df)*se),float(b+t.ppf(.975,df)*se)],
            "multiplicative_ratio_per_log2_port_unit":float(np.exp(b))}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--archive",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();raw=a.archive.read_bytes();need(sha(raw)==OUTER,"outer artifact hash")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        need(z.testzip() is None,"outer CRC")
        inner=z.read("GEB_V05_240_INTERNAL_STAGING_NOT_FOR_RELEASE.zip")
    need(sha(inner)==INNER,"inner hash")
    with zipfile.ZipFile(io.BytesIO(inner)) as z:
        need(z.testzip() is None,"inner CRC")
        def find(end):
            names=[n for n in z.namelist() if n.endswith(end)];need(len(names)==1,"member "+end)
            b=z.read(names[0])
            return pd.read_csv(io.BytesIO(b),dtype={"code":str,"level1":str}),sha(b)
        f,fh=find("s36/island_taxonomic_resolution_counts.csv")
        e,eh=find("botanical_effort_flora_island_inputs.csv")
    need(len(f)==126 and f.code.is_unique and f.level1.nunique()==9,"126 source")
    ef=f.merge(e[["code","gbif_all_plant_indexed_records","n_native_flora","n_introduced_flora","WPI_ports"]],on="code",validate="one_to_one")
    need(len(ef)==98 and ef.level1.nunique()==9,"98 effort source")
    need((ef.native_plant==ef.n_native_flora).all() and (ef.introduced_plant==ef.n_introduced_flora).all() and (ef.wpi_port_count==ef.WPI_ports).all(),"source agreement")
    arms=[("original98_effort",ef,False,True),("all126_basic",f,False,False),
          ("all126_plus_geography",f,True,False),("common98_geography_without_effort",ef,True,False),
          ("common98_geography_plus_effort",ef,True,True)]
    result={k:fit(d,g,eff) for k,d,g,eff in arms}
    need(np.isclose(result["original98_effort"]["beta"],0.24291347717458123,atol=1e-10,rtol=1e-8),"original log ratio coefficient")
    need(np.isclose(result["all126_basic"]["beta"],0.19001170891340988,atol=1e-10,rtol=1e-8),"paired unstandardized port slope difference")
    out={"schema":"chocho_WCVP126_port_origin_logratio_direct_v1",
         "status":"SOURCE_REPRODUCED_DIRECT_ORIGIN_RATIO_AND_CLUSTER_SENSITIVITY",
         "archive_id":15003375426,"outer_sha256":OUTER,"inner_sha256":INNER,
         "derived_input_sha256":{"126_region_panel":fh,"98_effort_panel":eh},
         "response":"log(1+introduced species count)-log(1+native species count)",
         "port_predictor":"log2(1+WPI ports), unstandardized",
         "all_models_controls":"log area and nine Level1 effects; plus log1p distance, absolute latitude and log1p all-plant indexed record volume where indicated",
         "original98_permutation_p_source_reported":0.0036,
         "original98_permutation_p_recomputed":False,
         "models":result,
         "limits":["New sensitivity reuses frozen derived panels, not new botanical records.",
          "Original permutation p is retained as a source result; clustered t8 and HC3 address different uncertainty assumptions and do not reproduce that permutation null.",
          "No positive direct port-origin contrast claimed under conservative cluster inference; no causal ship or mediation conclusion."]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2)+"\n")
    print("DIRECT_PORT_ORIGIN_RATIO_PASS",json.dumps(out))
if __name__=="__main__":main()

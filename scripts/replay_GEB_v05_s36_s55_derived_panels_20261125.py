#!/usr/bin/env python3
"""Recompute archived derived-panel fits using the exact original model functions."""
from __future__ import annotations
import argparse,hashlib,json,math,sys,textwrap
from pathlib import Path
from types import SimpleNamespace
import pandas as pd

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--tdwg-geojson",required=True,type=Path)
    a=p.parse_args()
    root=Path(__file__).resolve().parent
    sys.path.insert(0,str(root/"scripts"))
    import audit_fixed239_host_taxonomic_resolution_20261125 as h
    import audit_gbif224_island_observed_host_isolation_20261125 as o
    pins=json.loads((root/"SOURCE_PINS.json").read_text())
    assert hashlib.sha256(a.tdwg_geojson.read_bytes()).hexdigest()==pins["tdwg_geometry"]["sha256"]
    checks={}
    maxdiff={}
    def compare(expected,actual,label,path=""):
        if isinstance(expected,dict):
            assert isinstance(actual,dict) and set(expected)==set(actual),(label,path,"dict keys")
            for k in expected:compare(expected[k],actual[k],label,path+"/"+k)
        elif isinstance(expected,list):
            assert isinstance(actual,list) and len(expected)==len(actual),(label,path,"list length")
            for i,(x,y) in enumerate(zip(expected,actual)):compare(x,y,label,path+"/"+str(i))
        elif isinstance(expected,(float,int)) and not isinstance(expected,bool):
            assert isinstance(actual,(float,int)) and math.isclose(expected,actual,rel_tol=1e-8,abs_tol=1e-10),(label,path,expected,actual)
            maxdiff[label]=max(maxdiff.get(label,0.),abs(expected-actual))
            checks[label]=checks.get(label,0)+1
        else:
            assert expected==actual,(label,path,expected,actual)
            checks[label]=checks.get(label,0)+1
    hf=pd.read_csv(root/"s36/island_taxonomic_resolution_counts.csv",dtype={"code":str,"level1":str})
    he=json.loads((root/"s36/taxonomic_resolution_report.json").read_text())
    assert len(hf)==126 and hf.code.is_unique and hf.level1.nunique()==9
    assert (hf.native_species==hf.native_host_239).all() and (hf.current_species==hf.current_host_239).all()
    # Exact source code block after the source gate, through the fitted hierarchy,
    # composition, range and recording-effort sensitivities. No reimplemented fits.
    src=Path(h.__file__).read_text()
    start="    w,graph=make_weights(f,a.tdwg_geojson)"
    end="    nearfar=f.sort_values("
    assert src.count(start)==1 and src.count(end)==1
    block=textwrap.dedent(src[src.index(start):src.index(end)])
    scope=vars(h).copy()
    scope.update(f=hf,a=SimpleNamespace(tdwg_geojson=a.tdwg_geojson))
    exec(compile(block,"original_hierarchy_statistical_block","exec"),scope)
    for k,var in (("full","full"),("no_Fabaceae_or_Poaceae_butterflies","reduced"),
                  ("scratch_comparable_accepted_source","comparable"),("effort_common98_models","effort_models"),
                  ("single_native_region_proxy","end"),("graph","graph")):
        compare(he[k],scope[var],"S3.6",k)
    of=pd.read_csv(root/"s55/island_aggregates.csv",dtype={"code":str,"level1":str})
    oe=json.loads((root/"s55/recorded_island_host_isolation.json").read_text())["panels"]["all224"]
    assert len(of)==126 and of.code.is_unique and set(of.code)==set(hf.code)
    assert tuple(int(of["all1_"+k].sum()) for k in ("N","I","Z"))==(916,193,62)
    assert int((of.all1_I>0).sum())==69 and int(of.all2_I.sum())==143
    co={"all126":of,"original169_islands":of[of.primary169],
        "all_period_min20_days":of[of.all_days>=20],"late_period_min20_days":of[of.late_days>=20]}
    sums={name:{tag:o.describe(d,tag) for tag in ("all1","all2","all3","late1","late2","late3","prior1")} for name,d in co.items()}
    groups={name:{q:{tag:o.describe(d[d.distance_quartile==q],tag) for tag in ("all1","all2","late1","late2")} for q in ("near","far")} for name,d in co.items()}
    fits={name:{tag:o.models(d,tag,period) for tag,period in (("all1","all"),("all2","all"),("late1","late"),("late2","late"))} for name,d in co.items()}
    for k,v in (("summary",sums),("quartiles",groups),("models",fits)):compare(oe[k],v,"S5.5_all224",k)
    print("S36_S55_EXTRACTED_REPLAY_PASS",json.dumps({"checks":checks,"max_absolute_numeric_difference":maxdiff,
      "rtol":1e-8,"atol":1e-10,"scope":"derived-panel statistical replay; without-monarch and raw inputs not independently rebuilt here"}))
if __name__=="__main__":main()

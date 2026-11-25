#!/usr/bin/env python3
"""Exploratory fixed-horizon directional current x native source trait check.
No new ocean trajectories, no seed survival/arrival inference. The response is
the family-matched trait-positive minus trait-negative WCVP listing rate.
"""
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.stats import t as student_t
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups
import statsmodels.api as sm

G=Geod(ellps="WGS84")
STATUS="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE"
SEASONS=[907,919,931]
RESP=["native_excess_pp","alien_excess_pp","native_minus_alien_excess_pp"]
def distance(lon,lat,x,y):
    _,_,m=G.inv(np.broadcast_to(lon,np.shape(x)),np.broadcast_to(lat,np.shape(y)),x,y)
    return np.asarray(m)/1000
def scalar(v):
    return float(v) if np.isfinite(v) else None
def parse_args():
    p=argparse.ArgumentParser()
    for n in ["tracks","matrix","preflight","geography","coverage","paired","traits","out"]:
        p.add_argument("--"+n,type=Path,required=True)
    return p.parse_args()
def run():
    a=parse_args();a.out.mkdir(parents=True,exist_ok=True)
    tracks=json.loads(a.tracks.read_text())
    matrix=pd.read_csv(a.matrix)
    prep=json.loads(a.preflight.read_text())
    assert len(tracks)==280 and len({(x["source"],x["original_NOAA_source_time_index"]) for x in tracks})==280
    assert len(matrix)==15750 and not matrix.duplicated(["source","destination"]).any()
    ev=matrix[matrix.model_physical_source_evaluation_status==STATUS].copy()
    assert len(ev)==7313
    assert (ev.NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km>0).sum()==12
    pre={x["region"]:x for x in prep["126_source_region_release_eligibility"]}
    assert len(pre)==126
    geo=pd.read_csv(a.geography)
    assert len(geo)==7875 and not geo.duplicated(["source_region","destination_region"]).any()
    cols=["source_region","destination_region","centroid_geodesic_km","same_TDWG_Level1","source_representative_lat","source_representative_lon","dest_representative_lat","dest_representative_lon"]
    assert set(cols)<=set(geo)
    points={}
    for r in geo[cols].itertuples(index=False):
        for code,lat,lon in [(r.source_region,r.source_representative_lat,r.source_representative_lon),(r.destination_region,r.dest_representative_lat,r.dest_representative_lon)]:
            if not np.isfinite([lat,lon]).all():continue
            if code in points:assert np.allclose(points[code],[lat,lon],atol=1e-8,rtol=0),code
            points[code]=(lat,lon)
    gg=geo[cols].rename(columns={"source_region":"source","destination_region":"destination"})
    rev=gg.rename(columns={"source":"destination","destination":"source","source_representative_lat":"dest_representative_lat","source_representative_lon":"dest_representative_lon","dest_representative_lat":"source_representative_lat","dest_representative_lon":"source_representative_lon"})
    gg=pd.concat([gg,rev],ignore_index=True)
    assert len(gg)==15750 and not gg.duplicated(["source","destination"]).any()
    paired=pd.read_csv(a.paired)
    fixed22=set(paired.loc[paired["class"]=="island","code"])
    assert len(fixed22)==22
    cov=pd.read_csv(a.coverage)
    assert len(cov)==126
    traits=pd.read_csv(a.traits)
    assert not traits.duplicated(["source","destination"]).any()
    data=ev.merge(gg,on=["source","destination"],validate="one_to_one").merge(traits,on=["source","destination"],validate="one_to_one")
    data=data.replace([np.inf,-np.inf],np.nan).dropna(subset=["centroid_geodesic_km","source_representative_lat","dest_representative_lat","dest_representative_lon"]+RESP)
    data=data[data.centroid_geodesic_km>0].copy()
    data["same_level"]=data.same_TDWG_Level1.astype(str).str.lower().isin(["true","1","1.0","yes"]).astype(int)
    data["abs_lat_diff"]=abs(data.source_representative_lat-data.dest_representative_lat)
    # Source launch distance from original fixed GSHHG wet anchor, not point-to-port.
    eligible=[]
    for tr in tracks:
        src=tr["source"];idx=tr["original_NOAA_source_time_index"]
        assert idx in SEASONS
        off=pre[src]["nearest_valid_offshore_source_km_per_season"][SEASONS.index(idx)]
        assert off is not None and 0<=off<=200.0001
        hist={q["day"]:q for q in tr["trajectory_every_5days"]}
        assert 0 in hist
        if tr["complete_elapsed_days"]<30:continue
        assert 30 in hist
        q,p=hist[0],hist[30]
        assert np.isfinite([q["lat"],q["lon"],p["lat"],p["lon"]]).all()
        eligible.append((src,idx,float(off),q,p))
    # Fixed 30-day endpoint distance to the same frozen WCVP representative point.
    # A signed reduction is continuous proximity, not landfall or seed arrival.
    def progress(frame,offshore,minreleases):
        values=[]
        for src,g in frame.groupby("source",sort=True):
            ts=[x for x in eligible if x[0]==src and x[2]<=offshore]
            if len(ts)<minreleases:continue
            y=g.dest_representative_lat.to_numpy();x=g.dest_representative_lon.to_numpy()
            D=[];P=[];N=[]
            for _,idx,off,q,p in ts:
                d0=distance(q["lon"],q["lat"],x,y)
                d30=distance(p["lon"],p["lat"],x,y)
                D.append(d0);P.append(d0-d30)
                # A deterministic reverse endpoint with identical net displacement,
                # not a simulated reverse-current trajectory or valid water route.
                az,_,m=G.inv(q["lon"],q["lat"],p["lon"],p["lat"])
                rx,ry,_=G.fwd(q["lon"],q["lat"],az+180,m)
                N.append(d0-distance(rx,ry,x,y))
            o=g.copy()
            o["progress_km"]=np.mean(P,axis=0)
            o["reverse_endpoint_progress_km"]=np.mean(N,axis=0)
            o["initial_distance_km"]=np.mean(D,axis=0)
            o["valid_releases_30d"]=len(ts)
            values.append(o)
        return pd.concat(values,ignore_index=True) if values else frame.iloc[:0].copy()
    def fit(frame):
        # Geographic control includes source & destination fixed intercepts,
        # fixed cubic log-distance, same floristic realm and latitude contrast.
        # Equal source weighting, two-way node clustering, conservative small-cluster t.
        n=len(frame);ns=frame.source.nunique();nd=frame.destination.nunique()
        counts={"n_pairs":n,"source_regions":ns,"dest_regions":nd}
        if n<100 or ns<10 or nd<20:
            return {**counts,"status":"INSUFFICIENT_COVERAGE_NO_FIT"}
        p=frame.progress_km.to_numpy()
        if np.std(p)<1e-9:return {**counts,"status":"NO_PROGRESS_VARIATION_NO_FIT"}
        ds=np.log1p(frame.initial_distance_km.to_numpy());ds=(ds-ds.mean())/ds.std()
        src=pd.get_dummies(frame.source,drop_first=True,dtype=float)
        dst=pd.get_dummies(frame.destination,drop_first=True,dtype=float)
        X0=np.column_stack([np.ones(n),ds,ds**2,ds**3,frame.same_level,frame.abs_lat_diff,src.to_numpy(),dst.to_numpy()])
        # No full-rank requirement on redundant controls; residualization is unique.
        w=1/frame.groupby("source").source.transform("size").to_numpy()
        w=w/w.mean();sw=np.sqrt(w)
        def residual(v):return v-X0@np.linalg.lstsq(X0*sw[:,None],v*sw,rcond=None)[0]
        rp=residual(p/100)
        if np.sum(w*rp**2)<1e-8:return {**counts,"status":"NO_CONDITIONAL_PROGRESS_VARIATION_NO_FIT"}
        cl1=pd.factorize(frame.source)[0];cl2=pd.factorize(frame.destination)[0]
        df=min(ns,nd)-1
        out={}
        for name in RESP:
            y=frame[name].to_numpy();ry=residual(y)
            X=np.column_stack([rp])
            fit0=sm.WLS(ry,X,weights=w).fit()
            b=float(fit0.params[0])
            # Two-way cluster sandwich on FWL residuals; finite sample correction
            # is approximate and these are exploratory geographic associations.
            try:
                V=cov_cluster_2groups(fit0,cl1,cl2,use_correction=True)[0]
                se=float(np.sqrt(V[0,0])) if V[0,0]>0 else math.nan
            except Exception:se=math.nan
            ci=[scalar(b-student_t.ppf(.975,df)*se),scalar(b+student_t.ppf(.975,df)*se)]
            # Exact leave-one-source re-residualization on remaining nodes.
            loo=[]
            for i in range(ns):
                keep=cl1!=i
                xx=X0[keep];s=sw[keep]
                z=p[keep]/100;yy=y[keep]
                zp=z-xx@np.linalg.lstsq(xx*s[:,None],z*s,rcond=None)[0]
                yp=yy-xx@np.linalg.lstsq(xx*s[:,None],yy*s,rcond=None)[0]
                denom=np.sum(w[keep]*zp*zp)
                if denom>1e-9:loo.append(float(np.sum(w[keep]*zp*yp)/denom))
            # Directional orientation placebo: deterministic opposite endpoint.
            rr=residual(frame.reverse_endpoint_progress_km.to_numpy()/100)
            denom=np.sum(w*rr**2)
            b_rev=float(np.sum(w*rr*ry)/denom) if denom>1e-9 else None
            out[name]={"beta_pp_per_100km_progress":b,"approx_two_way_cluster_SE":scalar(se),"descriptive95":ci,"approx_two_way_t_p":scalar(2*student_t.sf(abs(b/se),df)) if se>0 else None,"leave_one_source_range":[min(loo),max(loo)] if loo else None,"leave_one_source_positive":sum(x>0 for x in loo),"leave_one_source_count":len(loo),"reverse_endpoint_placebo_beta":b_rev}
        assert abs(out[RESP[0]]["beta_pp_per_100km_progress"]-out[RESP[1]]["beta_pp_per_100km_progress"]-out[RESP[2]]["beta_pp_per_100km_progress"])<1e-8
        return {**counts,"status":"EXPLORATORY_FIT","positive_progress_pairs":int((p>0).sum()),"negative_progress_pairs":int((p<0).sum()),"progress_km_quantiles":np.quantile(p,[0,.25,.5,.75,1]).tolist(),"source_release_counts":frame.groupby("source").valid_releases_30d.first().astype(int).to_dict(),"results":out}
    variants=[("primary_fixed22_50km",50,1,True),("fixed22_50km_two_releases",50,2,True),("fixed22_200km",200,1,True),("all_trait_supported_50km",50,1,False),("all_trait_supported_200km",200,1,False)]
    results={}
    for label,off,mins,fixed in variants:
        frame=data[data.source.isin(fixed22)].copy() if fixed else data.copy()
        f=progress(frame,off,mins)
        results[label]=fit(f)
        print("VARIANT",label,json.dumps(results[label],allow_nan=False),flush=True)
        f.to_csv(a.out/(label+"_derived_pairs.csv"),index=False)
    receipt={"status":"CONTINUOUS_30DAY_CURRENT_TRAIT_ASSOCIATION_DIAGNOSTIC_COMPLETE","original_tracks":280,"tracks_with_exact_30d_endpoint":len(eligible),"original_evaluable_directions":7313,"original_nearshore_positive":12,"fixed_source_frame":sorted(fixed22),"trait_supported_sources_all126":int(cov.eligible.astype(str).str.lower().isin(["true","1"]).sum()),"analysis_definition":"Source-native accepted WCVP plants with known EuDiS trait; family-matched positive minus negative trait sharing in destination native and introduced floras; distance reduction to fixed representative point over first 30 valid days of archived 2005 OSCAR paths. Not arrival.","primary_gate":"10 sources,20 destinations,100 dyads; no relaxing after outcomes","results":results,"limitations":["2005 one-year coarse current and contemporary botanical occurrence; observational not causal.","Conditioning on 30-day survival in fluid model selects trajectories; censored routes excluded, not zero.","Source trait pool mostly European; WCVP botanical regions not individual islands.","Two-way cluster intervals approximate; no independent geographic transfer confirmation.","Source/destination FE and same matched species pools; not independent taxa or arrivals. No real reverse-current simulation.","The primary fixed22 source frame and cutoffs fixed before new output; earlier pilot already known, exploratory."]}
    (a.out/"RESULT.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("FINAL_RESULT",json.dumps(receipt,allow_nan=False),flush=True)
if __name__=="__main__":run()

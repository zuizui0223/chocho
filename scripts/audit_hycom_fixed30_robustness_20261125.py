#!/usr/bin/env python3
"""HYCOM fixed roster: same original family response and geographic model.
Only completed 30-day trajectories, equal date weighting within each source.
A frozen physical/technical gate prevents inference on inadequate support.
"""
import argparse,json,math,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.stats import t as student_t
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups
import statsmodels.api as sm
G=Geod(ellps="WGS84")
RESP=["native_excess_pp","alien_excess_pp","native_minus_alien_excess_pp"]
def distance(lon,lat,x,y):
    _,_,m=G.inv(np.broadcast_to(lon,np.shape(x)),np.broadcast_to(lat,np.shape(y)),x,y)
    return np.asarray(m)/1000
def scalar(v):return float(v) if np.isfinite(v) else None
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


def robustness(frame):
    """Post-result sensitivity audit, not a new confirmatory hypothesis test."""
    n=len(frame);ns=frame.source.nunique();nd=frame.destination.nunique()
    if n<100 or ns<10 or nd<20:return {"status":"INSUFFICIENT_COVERAGE_NO_FIT","n_pairs":n,"source_regions":ns,"dest_regions":nd}
    sources=sorted(frame.source.unique());dests=sorted(frame.destination.unique())
    ds=np.log1p(frame.initial_distance_km.to_numpy());ds=(ds-ds.mean())/ds.std()
    src=pd.get_dummies(frame.source,drop_first=True,dtype=float)
    dst=pd.get_dummies(frame.destination,drop_first=True,dtype=float)
    X0=np.column_stack([np.ones(n),ds,ds**2,ds**3,frame.same_level,frame.abs_lat_diff,src.to_numpy(),dst.to_numpy()])
    w=1/frame.groupby("source").source.transform("size").to_numpy()
    w=w/w.mean();sw=np.sqrt(w)
    # Orthonormal weighted column space, robust to redundant fixed-effect columns.
    U,S,_=np.linalg.svd(X0*sw[:,None],full_matrices=False)
    rank=int((S>S[0]*max(X0.shape)*np.finfo(float).eps).sum())
    Q=U[:,:rank]
    def residual(v):
        z=sw*v
        return z-Q@(Q.T@z)
    p=frame.progress_km.to_numpy()/100
    reverse=frame.reverse_endpoint_progress_km.to_numpy()/100
    x=residual(p);xx=float(x@x)
    if xx<1e-8:return {"status":"NO_CONDITIONAL_PROGRESS_VARIATION_NO_FIT"}
    a=frame.source.to_numpy();b=frame.destination.to_numpy()
    srcidx=pd.factorize(a)[0];dstidx=pd.factorize(b)[0]
    nodes=sorted(set(a)|set(b));unord=["|".join(sorted([s,d])) for s,d in zip(a,b)]
    # Node-overlap covariance: each pair of residualized scores with a shared
    # endpoint contributes once, regardless of source/recipient role.
    def node_var(score):
        u=pd.DataFrame({"node":np.r_[a,b],"score":np.r_[score,score]}).groupby("node").score.sum()
        v=pd.DataFrame({"dyad":unord,"score":score}).groupby("dyad").score.sum()
        return float((np.sum(u**2)-np.sum(v**2))/xx**2)
    # Enumerate the exact orientation-placebo support: each source independently
    # takes actual or opposite-bearing endpoint progress. This is a geometric
    # falsification diagnostic, not an exchangeable physical/causal null.
    base=residual(reverse)
    deltas=np.column_stack([residual(np.where(a==s,p-reverse,0.)) for s in sources])
    bits=((np.arange(1<<ns,dtype=np.uint32)[:,None]>>np.arange(ns,dtype=np.uint32))&1).astype(float)
    orient=base[:,None]+deltas@bits.T
    norms=np.sqrt(np.sum(orient**2,axis=0))
    assert np.allclose(orient[:,-1],x,atol=1e-8,rtol=1e-8)
    out={}
    for name in RESP:
        y=frame[name].to_numpy();ry=residual(y);beta=float(x@ry/xx)
        err=ry-beta*x;score=x*err
        shared=node_var(score)
        ndof=n-rank-1
        # Full design cluster covariance and rank-aware finite sample adjustment
        # are sensitivity diagnostics alongside the original FWL interval.
        model=sm.WLS(y,np.column_stack([X0,p]),weights=w).fit()
        assert abs(model.params[-1]-beta)<1e-7
        c0=cov_cluster_2groups(model,srcidx,dstidx,use_correction=False)[0][-1,-1]
        c1=cov_cluster_2groups(model,srcidx,dstidx,use_correction=True)[0][-1,-1]
        se=float(np.sqrt(c1)) if c1>0 else math.nan
        shared_se=float(np.sqrt(shared*n/ndof)) if shared>0 and ndof>0 else math.nan
        df=min(ns,nd)-1;tcrit=student_t.ppf(.975,df)
        loo={}
        leverage={}
        for s in sources:
            k=a!=s
            x0=X0[k];sw0=sw[k]
            zz=p[k]*sw0;yy=y[k]*sw0
            Uz,Sz,_=np.linalg.svd(x0*sw0[:,None],full_matrices=False)
            rz=int((Sz>Sz[0]*max(x0.shape)*np.finfo(float).eps).sum())
            q=Uz[:,:rz];xr=zz-q@(q.T@zz);yr=yy-q@(q.T@yy)
            loo[s]=float(xr@yr/(xr@xr)) if xr@xr>1e-9 else None
            leverage[s]={"residual_predictor_ss_fraction":float(np.sum(x[a==s]**2)/xx),"partial_numerator_fraction":float(np.sum((x*ry)[a==s])/(x@ry)) if abs(x@ry)>1e-12 else None}
        # Partial-correlation score keeps scale comparable under endpoint flips.
        normy=float(np.linalg.norm(ry))
        scores=(ry@orient)/(normy*norms)
        valid=np.isfinite(scores);assert valid[-1]
        actual=float(scores[-1]);tail=float(np.mean(scores[valid]>=actual-1e-12));twoside=float(np.mean(np.abs(scores[valid])>=abs(actual)-1e-12))
        out[name]={"beta":beta,"full_WLS_two_way_cluster_variance_uncorrected":scalar(c0),"full_WLS_two_way_cluster_SE":scalar(se),"full_WLS_two_way_cluster_95":[scalar(beta-tcrit*se),scalar(beta+tcrit*se)],"full_WLS_two_way_cluster_p":scalar(2*student_t.sf(abs(beta/se),df)) if se>0 else None,"shared_endpoint_variance_uncorrected":scalar(shared),"shared_endpoint_SE_rank_adjusted":scalar(shared_se),"shared_endpoint_95":[scalar(beta-tcrit*shared_se),scalar(beta+tcrit*shared_se)],"shared_endpoint_p_approx":scalar(2*student_t.sf(abs(beta/shared_se),df)) if shared_se>0 else None,"leave_one_source":loo,"source_leverage":leverage,"orientation_placebo":{"configurations":int(valid.sum()),"actual_partial_correlation":actual,"positive_tail_fraction":tail,"two_sided_tail_fraction":twoside,"all_reverse_partial_correlation":scalar(scores[0])}}
    assert abs(out[RESP[0]]["beta"]-out[RESP[1]]["beta"]-out[RESP[2]]["beta"])<1e-8
    return {"status":"POST_RESULT_ROBUSTNESS_COMPLETE","n_pairs":n,"source_regions":ns,"dest_regions":nd,"weighted_controls_rank":rank,"residual_df":n-rank-1,"results":out}

def main():
    p=argparse.ArgumentParser()
    for name in ("physical","tracks","anchors","geography","coverage","paired","traits","out","expected"):p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    phy=json.loads(a.physical.read_text());assert phy["physical_gate_pass"] and phy["technical_gate_pass"]
    tracks=json.loads(a.tracks.read_text());assert len(tracks)==306
    valid=[t for t in tracks if t["day30_valid"]]
    assert len(valid)==phy["day30_valid"] and len({t["source"] for t in valid})==phy["source_regions_valid30"]
    ac=json.loads(a.anchors.read_text())["coastal_anchors_by_TDWG_region"]
    dests={k for k,v in ac.items() if v};assert len(dests)==104
    geo=pd.read_csv(a.geography)
    assert len(geo)==7875 and not geo.duplicated(["source_region","destination_region"]).any()
    cols=["source_region","destination_region","centroid_geodesic_km","same_TDWG_Level1","source_representative_lat","source_representative_lon","dest_representative_lat","dest_representative_lon"]
    assert set(cols)<=set(geo)
    gg=geo[cols].rename(columns={"source_region":"source","destination_region":"destination"})
    rev=gg.rename(columns={"source":"destination","destination":"source","source_representative_lat":"dest_representative_lat","source_representative_lon":"dest_representative_lon","dest_representative_lat":"source_representative_lat","dest_representative_lon":"source_representative_lon"})
    gg=pd.concat([gg,rev],ignore_index=True)
    assert len(gg)==15750 and not gg.duplicated(["source","destination"]).any()
    paired=pd.read_csv(a.paired);fixed22=set(paired.loc[paired["class"]=="island","code"]);assert len(fixed22)==22
    coverage=pd.read_csv(a.coverage)
    supported=set(coverage.loc[coverage.eligible.astype(str).str.lower().isin(["true","1"]),"source"])
    assert len(supported)==17 and supported<=fixed22
    traits=pd.read_csv(a.traits)
    assert len(traits)==2125 and not traits.duplicated(["source","destination"]).any()
    data=gg[gg.source.isin(supported)&gg.destination.isin(dests)].merge(traits,on=["source","destination"],validate="one_to_one")
    data=data.replace([np.inf,-np.inf],np.nan).dropna(subset=["centroid_geodesic_km","source_representative_lat","dest_representative_lat","dest_representative_lon"]+RESP)
    data=data[data.centroid_geodesic_km>0].copy()
    data["same_level"]=data.same_TDWG_Level1.astype(str).str.lower().isin(["true","1","1.0","yes"]).astype(int)
    data["abs_lat_diff"]=abs(data.source_representative_lat-data.dest_representative_lat)
    results={}
    diagnostics={}
    expected=json.loads(a.expected.read_text())
    for label,minimum_dates in [("primary_fixed22_50km_hycom",1),("two_distinct_valid_release_dates",2)]:
        rows=[]
        for src,g in data.groupby("source",sort=True):
            ts=[t for t in valid if t["source"]==src]
            dates=sorted({t["date"] for t in ts})
            if len(dates)<minimum_dates:continue
            y=g.dest_representative_lat.to_numpy();x=g.dest_representative_lon.to_numpy()
            P=[];D=[];N=[]
            for d in dates:
                pp=[];dd=[];nn=[]
                for t in [v for v in ts if v["date"]==d]:
                    q=t["start"];p={"lon":t["end_lon"],"lat":t["end_lat"]}
                    assert q["offset_km"]<=50.00001 and t["elapsed_hours"]==720
                    d0=distance(q["lon"],q["lat"],x,y)
                    d30=distance(p["lon"],p["lat"],x,y)
                    az,_,m=G.inv(q["lon"],q["lat"],p["lon"],p["lat"])
                    rx,ry,_=G.fwd(q["lon"],q["lat"],az+180,m)
                    pp.append(d0-d30);dd.append(d0);nn.append(d0-distance(rx,ry,x,y))
                P.append(np.mean(pp,axis=0));D.append(np.mean(dd,axis=0));N.append(np.mean(nn,axis=0))
            z=g.copy();z["progress_km"]=np.mean(P,axis=0);z["initial_distance_km"]=np.mean(D,axis=0);z["reverse_endpoint_progress_km"]=np.mean(N,axis=0)
            z["valid_releases_30d"]=len(dates);z["valid_tracks_30d"]=len(ts)
            rows.append(z)
        frame=pd.concat(rows,ignore_index=True) if rows else data.iloc[:0].copy()
        result=fit(frame)
        result["valid_release_date_counts"]=frame.groupby("source").valid_releases_30d.first().astype(int).to_dict() if len(frame) else {}
        result["valid_track_counts"]=frame.groupby("source").valid_tracks_30d.first().astype(int).to_dict() if len(frame) else {}
        results[label]=result
        old=expected["results"][label]
        for k in ("n_pairs","source_regions","dest_regions","positive_progress_pairs","negative_progress_pairs"):assert result[k]==old[k],(label,k)
        for k in RESP:assert abs(result["results"][k]["beta_pp_per_100km_progress"]-old["results"][k]["beta_pp_per_100km_progress"])<1e-9,(label,k)
        checks={"all_pairs":robustness(frame)}
        for threshold in (5000,3000,1500):
            checks["initial_distance_le_"+str(threshold)+"km"]=robustness(frame[frame.initial_distance_km<=threshold].copy())
        checks["same_TDWG_Level1"]=robustness(frame[frame.same_level==1].copy())
        diagnostics[label]=checks
        for key,d in checks.items():
            brief={k:v for k,v in d.items() if k!="results"}
            if "results" in d:
                brief["results"]={k:{f:v[f] for f in ("beta","full_WLS_two_way_cluster_95","shared_endpoint_95","orientation_placebo","leave_one_source")} for k,v in d["results"].items()}
            print("HYCOM_ROBUSTNESS",label,key,json.dumps(brief,allow_nan=False),flush=True)
    receipt={"status":"HYCOM_EXPLORATORY_CURRENT_TRAIT_DIAGNOSTIC_COMPLETE","physical_tracks":306,"valid_day30":len(valid),"source_regions_valid30":phy["source_regions_valid30"],"response":"Same source-native accepted plants; within-family positive-minus-negative trait sharing into native and introduced destination flora; paired difference in percentage points.","predictor":"Mean signed distance reduction to fixed recipient representative point at day 30, equal date weights and equal starts within each valid date, per 100 km.","geographic_scope":"Original 17 family-supported source regions within frozen22; 104 original GSHHG anchored destinations; exclude self, missing coordinates and invalid source tracks.","results":results,"interpretation_limits":["New post-pilot exploratory HYCOM physical design, one year 2005, 15m Eulerian model velocities, no seed observations or viability.","Only valid water tracks at 30 days; selection may favor offshore flow and omit true nearshore settlement.","Equal-date weighting, source and destination FE; source-cluster count is small and species/dyads overlap.","Reverse endpoint is a geometric placebo, not reversed-flow simulation.","Do not claim causal native transport or actual shipping. No automatic promotion to GEB manuscript."]}
    audit={"status":"HYCOM_POST_RESULT_ROBUSTNESS_AUDIT_COMPLETE","original_model_coefficients_reproduced":True,"original_result":receipt,"diagnostics":diagnostics,"scope_limits":["All new diagnostics are post-result exploratory, no preregistration claim.","Shared-node covariance handles overlapping geography in either role, not shared taxa across all nonoverlapping dyads.","The orientation flip enumerates geometric opposite endpoints, not physical reversed trajectories or a justified causal randomization null.","Distance/same-level sensitivities retain the original 10-source minimum; failed gates are not effect zero.","No change to manuscript."]}
    (a.out/"HYCOM_ROBUSTNESS_RESULT.json").write_text(json.dumps(audit,indent=2,allow_nan=False)+"\n")
    print("HYCOM_ROBUSTNESS_PASS",json.dumps({"original_coefficients_reproduced":True,"variants":len(diagnostics)}),flush=True)
if __name__=="__main__":main()

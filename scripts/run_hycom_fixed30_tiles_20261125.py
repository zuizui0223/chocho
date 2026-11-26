#!/usr/bin/env python3
"""Fixed HYCOM launch roster, daily 15m fields, 6h RK2, 30d endpoints.
Tiles fetched only on demand; no imputation, no claim of viable seed arrival.
"""
import argparse,collections,datetime as dt,hashlib,json,math,time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
import netCDF4
import numpy as np
from pyproj import Geod
G=Geod(ellps="WGS84")
ENDPOINT="https://ncss.hycom.org/thredds/ncss/grid/GLBu0.08/expt_19.1/2005"
DATES=["2005-05-28","2005-06-27","2005-07-27","2005-08-27","2005-09-26","2005-10-27"]
GRAPH_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
HOURS=6
NSTEPS=120
SIZE=4.
HALO=.24
MAX_TILES=120
def sha(b):return hashlib.sha256(b).hexdigest()
def wrap(x):return (x+180)%360-180
def dist(a,b):return abs(G.inv(a[1],a[0],b[1],b[0])[2])/1000.
def wet(mask,y,x):
    r=int(round((85-y)/.0625-.5));c=int(round((x+180)%360/.0625-.5))%mask.shape[1]
    return 0<=r<mask.shape[0] and bool(mask[r,c]>=0)
def advance(y,x,u,v,sec):
    # Match the original spherical displacement integrator; 6h substeps.
    R=6371008.8
    if abs(math.cos(math.radians(y)))<1e-4:raise ValueError("polar displacement")
    return y+math.degrees(v*sec/R),wrap(x+math.degrees(u*sec/(R*math.cos(math.radians(y)))))
def interpolate(T,Y,X,U,V,t,y,x):
    if not(T[0]<=t<=T[-1] and Y[-1]<=y<=Y[0] and X[0]<=x<=X[-1]):return None
    i=max(0,min(len(T)-2,int(np.searchsorted(T,t,side="right")-1)))
    j=max(0,min(len(Y)-2,int(np.searchsorted(-Y,-y,side="right")-1)))
    k=max(0,min(len(X)-2,int(np.searchsorted(X,x,side="right")-1)))
    td=(t-T[i])/(T[i+1]-T[i]);yd=(Y[j]-y)/(Y[j]-Y[j+1]);xd=(x-X[k])/(X[k+1]-X[k])
    w=np.array([(1-yd)*(1-xd),(1-yd)*xd,yd*(1-xd),yd*xd])
    out=[]
    for field in (U,V):
        a=field[i:i+2,j:j+2,k:k+2].reshape(2,4)
        if not np.isfinite(a).all():return None
        out.append(float((1-td)*np.dot(a[0],w)+td*np.dot(a[1],w)))
    if max(abs(out[0]),abs(out[1]))>5:raise ValueError("implausible velocity")
    return out
def fixture():
    T=np.array([0.,86400.]);Y=np.array([1.,0.]);X=np.array([10.,11.])
    U=np.full((2,2,2),.2);V=np.zeros_like(U)
    assert np.allclose(interpolate(T,Y,X,U,V,43200,.5,10.5),[.2,0])
    U[0,0,0]=np.nan;assert interpolate(T,Y,X,U,V,43200,.5,10.5) is None
    p=advance(.5,10.5,.2,0,21600)
    assert p[1]>10.5 and abs(p[0]-.5)<1e-12
    assert abs(dist((.5,10.5),p)-4.32)<.02
class Fields:
    def __init__(self,date):
        self.date=date
        self.t0=dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc)
        self.days=[(self.t0+dt.timedelta(days=i)).strftime("%Y-%m-%d") for i in range(31)]
        self.cache={};self.receipts=[];self.failures={}
    def fetch_box(self,w,e,s,n,center):
        args=[("var","water_u"),("var","water_v"),("north",str(n)),("south",str(s)),("west",str(w)),("east",str(e)),("time_start",self.days[0]+"T00:00:00Z"),("time_end",self.days[-1]+"T00:00:00Z"),("vertCoord","15"),("accept","netcdf4"),("addLatLon","true")]
        url=ENDPOINT+"?"+urlencode(args)
        last=None
        for attempt in range(3):
            try:
                with urlopen(Request(url,headers={"User-Agent":"chocho-HYCOM-fixed30-reproducible-research/1.0"}),timeout=100) as r:
                    b=r.read(16_000_001);assert r.status==200 and len(b)<=16_000_000
                assert b.startswith(b"CDF") or b.startswith(b"\x89HDF\r\n\x1a\n")
                with netCDF4.Dataset("tile.nc",memory=b) as ds:
                    def name(opts):return next(k for k in opts if k in ds.variables)
                    def arr(k):
                        q=ds[k][:];return np.asarray(q.filled(np.nan) if np.ma.isMaskedArray(q) else q,dtype=float)
                    tn=name(("time","MT","Time"));dn=name(("depth","Depth","Depth1","depth1"));yn=name(("lat","latitude","Latitude"));xn=name(("lon","longitude","Longitude"))
                    d=arr(dn).reshape(-1);assert len(d)==1 and abs(d[0]-15)<.01
                    times=netCDF4.num2date(ds[tn][:],ds[tn].units,calendar=getattr(ds[tn],"calendar","standard"))
                    assert [str(x)[:10] for x in times]==self.days
                    assert all(str(x)[11:19]=="00:00:00" for x in times)
                    dims=[ds[k].dimensions[0] for k in (tn,dn,yn,xn)]
                    def vel(k):
                        q=ds[k];assert len(q.dimensions)==4 and set(q.dimensions)==set(dims)
                        return np.transpose(arr(k),[q.dimensions.index(x) for x in dims])[:,0]
                    Y=arr(yn);X=arr(xn);U=vel("water_u");V=vel("water_v")
                iy=np.argsort(-Y);Y=Y[iy];X=center+(X-center+180)%360-180;ix=np.argsort(X);X=X[ix];U=U[:,iy][:,:,ix];V=V[:,iy][:,:,ix]
                assert U.shape==V.shape==(31,len(Y),len(X))
                assert min(len(Y),len(X))>=2 and np.allclose(np.diff(Y),-.08,atol=.002) and np.allclose(np.diff(X),.08,atol=.002)
                self.receipts.append({"url":url,"sha256":sha(b),"bytes":len(b),"shape":list(U.shape),"valid_uv":int((np.isfinite(U)&np.isfinite(V)).sum())})
                return Y,X,U.astype(np.float32),V.astype(np.float32)
            except Exception as ex:
                last=ex
                if attempt<2:time.sleep(2*(attempt+1))
        raise RuntimeError(str(last)[:300])
    def tile(self,y,x):
        x=wrap(x);j=int(math.floor(y/SIZE));k=int(math.floor((x+180)/SIZE))
        key=(j,k)
        if key in self.failures:return None
        if key not in self.cache:
            if len(self.cache)+len(self.failures)>=MAX_TILES:
                self.failures[key]="TILE_BUDGET_EXHAUSTED";return None
            s=j*SIZE-HALO;n=(j+1)*SIZE+HALO;w=k*SIZE-180-HALO;e=(k+1)*SIZE-180+HALO;center=(w+e)/2
            try:
                # HYCOM longitude wrap requests split at the date line.
                boxes=[(w,e)]
                if w< -180:boxes=[(w+360,180),(-180,e)]
                elif e>180:boxes=[(w,180),(-180,e-360)]
                parts=[self.fetch_box(a,b,s,n,center) for a,b in boxes]
                Y=parts[0][0]
                assert all(np.array_equal(Y,p[0]) for p in parts)
                X=np.concatenate([p[1] for p in parts]);U=np.concatenate([p[2] for p in parts],axis=2);V=np.concatenate([p[3] for p in parts],axis=2)
                ix=np.argsort(X);X=X[ix];U=U[:,:,ix];V=V[:,:,ix]
                keep=np.r_[True,np.diff(X)>1e-5];X=X[keep];U=U[:,:,keep];V=V[:,:,keep]
                assert np.allclose(np.diff(X),.08,atol=.002)
                self.cache[key]=(np.arange(31,dtype=float)*86400,Y,X,U,V,center)
                print("HYCOM_TILE",json.dumps({"date":self.date,"tile":list(key),"shape":list(U.shape),"tiles":len(self.cache)}),flush=True)
            except Exception as ex:
                self.failures[key]=str(ex)[:300]
                print("HYCOM_TILE_FAILURE",json.dumps({"tile":list(key),"error":self.failures[key]}),flush=True)
                return None
        return self.cache[key]
    def vel(self,t,y,x):
        if not math.isfinite(y+x+t) or abs(y)>=80:return None,"OUTSIDE_MODEL_LATITUDE"
        a=self.tile(y,x)
        if a is None:return None,"SOURCE_TILE_UNAVAILABLE"
        T,Y,X,U,V,c=a
        xx=c+(x-c+180)%360-180
        if not(Y[-1]<=y<=Y[0] and X[0]<=xx<=X[-1]):return None,"OUTSIDE_RETRIEVED_TILE_GRID"
        q=interpolate(T,Y,X,U,V,t,y,xx)
        return q,None if q is not None else "MASKED_OR_MISSING_VELOCITY"
def run_track(fields,mask,q,src,index):
    start=(q["lat"],q["lon"]);current=start
    assert q["offset_km"]<=50.00001 and wet(mask,*start)
    stepsec=HOURS*3600;done=0;length=0.;status="COMPLETED_30_DAYS"
    for j in range(NSTEPS):
        t=j*stepsec
        v,err=fields.vel(t,*current)
        if err:status="INITIAL_"+err;break
        half=advance(*current,*v,stepsec/2)
        mid,err=fields.vel(t+stepsec/2,*half)
        if err:status="MIDPOINT_"+err;break
        nxt=advance(*current,*mid,stepsec)
        n=max(1,int(math.ceil(dist(current,nxt)/2.)))
        hit=False
        for k in range(1,n+1):
            a=k/n
            y=current[0]+a*(nxt[0]-current[0]);x=wrap(current[1]+a*wrap(nxt[1]-current[1]))
            if not wet(mask,y,x):hit=True;break
        if hit:status="GSHHG_LAND_CONTACT_BEFORE_DAY30";break
        length+=dist(current,nxt);current=nxt;done=j+1
    return {"source":src,"date":fields.date,"start_index":index,"start":q,"elapsed_hours":done*HOURS,"day30_valid":done==NSTEPS,"status":status,"end_lat":current[0],"end_lon":current[1],"path_length_km":length,"net_displacement_km":dist(start,current)}
def main():
    a=argparse.ArgumentParser();a.add_argument("--preflight",type=Path,required=True);a.add_argument("--graph",type=Path,required=True);a.add_argument("--date",choices=DATES,required=True);a.add_argument("--out",type=Path,required=True);a=a.parse_args()
    fixture()
    r=json.loads(a.preflight.read_text())
    assert r["launch_gate_pass"] and r["source_regions"]==17 and r["total_starts"]==306 and r["all_306_anchor_date_queries_verified"]
    assert r["fixed_dates"]==DATES and sha(a.graph.read_bytes())==GRAPH_SHA
    with np.load(a.graph) as g:mask=g["fine_cell_to_node"]
    assert mask.shape==(2720,5760)
    f=Fields(a.date);rows=[]
    for s in sorted(r["region_launches"],key=lambda x:x["source"]):
        q=[d for d in s["dates"] if d["date"]==a.date];assert len(q)==1 and len(q[0]["starts"])==3
        for i,start in enumerate(q[0]["starts"]):
            z=run_track(f,mask,start,s["source"],i+1);rows.append(z)
            print("HYCOM_TRACK",json.dumps({k:v for k,v in z.items() if k!="start"}),flush=True)
    assert len(rows)==51 and len({(x["source"],x["start_index"]) for x in rows})==51
    a.out.mkdir(parents=True,exist_ok=True)
    payload={"status":"HYCOM_FIXED30_PHYSICAL_TRACKS_NOT_BIOLOGICAL_ARRIVAL","date":a.date,"tracks":rows,"launches":51,"valid_day30":sum(x["day30_valid"] for x in rows),"source_regions_valid_day30":len({x["source"] for x in rows if x["day30_valid"]}),"termination_counts":dict(collections.Counter(x["status"] for x in rows)),"velocity_tile_receipts":f.receipts,"failed_tiles":[{"key":list(k),"error":v} for k,v in f.failures.items()],"source_tile_budget":MAX_TILES,"step_hours":HOURS,"horizon_days":30,"depth_m":15,"raw_velocity_exported":False}
    (a.out/"TRACKS.json").write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")
    print("HYCOM_DATE_RESULT",json.dumps({k:v for k,v in payload.items() if k not in ("tracks","velocity_tile_receipts")}),flush=True)
if __name__=="__main__":main()

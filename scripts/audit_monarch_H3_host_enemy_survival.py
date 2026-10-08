#!/usr/bin/env python3
"""Two-year Greenstein-H3 host x cage effects on monarch adult emergence.

Original ecological study: Diethelm et al. 2026 Dryad 10.5061/dryad.51c59zwgx.
Already-published growth and morphology interaction is NOT an original discovery.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,random,time,urllib.error,urllib.request,zipfile
from pathlib import Path

FILES={2018:("2018_VR_Predation_CG_subsetAfAs.csv",4697215),
       2019:("2019_VR_Predation_CG_subsetAfAs.csv",4697214)}
DOI="10.5061/dryad.51c59zwgx"
BASE="https://datadryad.org"
NEEDED=("milkweed_spp","treatment","survived_binary")

def validated_csv(raw):
    if b"\x00" in raw[:1000]:return False
    try:head=next(csv.reader(io.StringIO(raw.decode("utf-8-sig",errors="replace"))))
    except Exception:return False
    return all(c in head for c in NEEDED)

def get(url):
    request=urllib.request.Request(url,headers={
       "User-Agent":"chocho-butterfly-source-repro/0.1 (github.com/zuizui0223/chocho)",
       "Accept":"text/csv,application/octet-stream,*/*"})
    with urllib.request.urlopen(request,timeout=90) as f:return f.read()

def download(source_dir):
    source_dir.mkdir(parents=True,exist_ok=True)
    source_log={}
    for year,(filename,file_id) in FILES.items():
        urls=[
          f"{BASE}/api/v2/files/{file_id}/download",
          f"{BASE}/downloads/file_stream/{file_id}"
        ]
        errors=[];data=None;selected=None
        for url in urls:
            for rep in range(2):
                try:
                    payload=get(url)
                    if not validated_csv(payload):raise ValueError("download is not original scientific data CSV schema")
                    data=payload;selected=url;break
                except Exception as e:
                    errors.append(f"{url} #{rep+1}: {type(e).__name__}: {str(e)[:180]}")
                    time.sleep(1+rep)
            if data is not None:break
        if data is None:
            # Try DOI-identified full public dataset export without guessing content.
            try:
                payload=get(f"{BASE}/api/v2/datasets/doi%3A10.5061%2Fdryad.51c59zwgx/download")
                with zipfile.ZipFile(io.BytesIO(payload)) as z:
                    candidates=[x for x in z.namelist() if Path(x).name==filename]
                    if len(candidates)!=1:raise RuntimeError(f"full dataset missing exactly one {filename}")
                    raw=z.read(candidates[0])
                if not validated_csv(raw):raise RuntimeError("dataset ZIP member schema mismatch")
                data=raw;selected="Dryad public DOI dataset ZIP"
            except Exception as e:errors.append(f"DOI bulk ZIP: {type(e).__name__}: {str(e)[:180]}")
        if data is None:
            raise RuntimeError("None of the source-verified Dryad public file endpoints worked: "+repr(errors))
        (source_dir/filename).write_bytes(data)
        source_log[year]={"filename":filename,"source":selected,"bytes":len(data),
                          "sha256":hashlib.sha256(data).hexdigest(),"errors_on_alternate_urls":errors}
        print(json.dumps({"year":year,"source":selected,"bytes":len(data)}),flush=True)
    (source_dir/"source_retrieval.json").write_text(json.dumps({"doi":DOI,"sources":source_log},indent=2)+"\n")

def read_data(source_dir):
    all_rows=[];source_audit={}
    for year,(filename,_) in FILES.items():
        with (source_dir/filename).open(newline="",encoding="utf-8-sig") as f:
            reader=csv.DictReader(f)
            if not set(NEEDED)<=set(reader.fieldnames or []):raise RuntimeError("source columns unavailable")
            originals=list(reader)
        excluded={}
        keep=[]
        seen=set()
        for row in originals:
            host=(row.get("milkweed_spp") or "").strip()
            treat=(row.get("treatment") or "").strip().lower()
            outcome=(row.get("survived_binary") or "").strip()
            uid=(row.get("monarch_id") or "").strip()
            if host not in ("Af","As") or treat not in ("cage","roof") or outcome not in ("0","1"):
                key="missing/invalid host_treatment_or_survival";excluded[key]=excluded.get(key,0)+1;continue
            if not uid or uid in seen:raise RuntimeError(f"missing/duplicate original monarch_id in {year}")
            seen.add(uid)
            keep.append({"year":year,"host":host,"cage":int(treat=="cage"),
                         "survival":int(outcome),"monarch_id":uid})
        source_audit[year]={"rows_original":len(originals),"rows_retained":len(keep),
                            "rows_excluded":excluded,"distinct_monarch_ids":len(seen)}
        all_rows.extend(keep)
    return all_rows,source_audit

def empirical(rows):
    years=sorted({r["year"] for r in rows})
    cells={}
    for year in years:
        for host in ("As","Af"):
            for cage in (1,0):
                x=[r["survival"] for r in rows if (r["year"],r["host"],r["cage"])==(year,host,cage)]
                if len(x)<10:raise RuntimeError(f"undersized year x host x treatment cell {year} {host} {cage} {len(x)}")
                cells[(year,host,cage)]={"n":len(x),"adults":sum(x),"p":sum(x)/len(x),"individual":x}
    year_weight={year:sum(v["n"] for k,v in cells.items() if k[0]==year) for year in years}
    N=sum(year_weight.values())
    def diff(d):
        risks={}
        for year in years:
            cage_benefit={host:d[(year,host,1)]-d[(year,host,0)] for host in ("As","Af")}
            risks[year]=cage_benefit["As"]-cage_benefit["Af"]
        return sum(year_weight[y]/N*risks[y] for y in years),risks
    p={k:v["p"] for k,v in cells.items()}
    observed,by_year=diff(p)
    rng=random.Random(20261008)
    draws=[]
    for _ in range(3999):
        b={}
        for key,cell in cells.items():
            a=cell["individual"];n=len(a)
            b[key]=sum(a[rng.randrange(n)] for _ in range(n))/n
        draws.append(diff(b)[0])
    draws.sort()
    def quantile(q):
        z=(len(draws)-1)*q;i=int(z);j=math.ceil(z)
        return draws[i] if i==j else draws[i]*(j-z)+draws[j]*(z-i)
    totals={}
    for host in ("As","Af"):
        for cage in (1,0):
            x=[r["survival"] for r in rows if (r["host"],r["cage"])==(host,cage)]
            totals[f"{host}_{'cage' if cage else 'roof'}"]={"n":len(x),"adult_survivors":sum(x),"risk":sum(x)/len(x)}
    return {"weighted_risk_difference_in_differences":observed,
           "year_specific_difference_in_differences":by_year,
           "bootstrap_ci95":[quantile(.025),quantile(.975)],
           "all_years_direction_consistent":len({(v>0)-(v<0) for v in by_year.values() if v!=0})==1,
           "by_year_host_treatment":[{"year":y,"host":h,"treatment":"cage" if c else "roof",
                      "n":v["n"],"adult_survivors":v["adults"],"survival_rate":v["p"]}
                       for (y,h,c),v in sorted(cells.items())],
           "pooled_survival_by_host_treatment":totals}

def glm_test(rows):
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import chi2
    y=np.asarray([r["survival"] for r in rows],dtype=float)
    year=np.asarray([int(r["year"]==2019) for r in rows],dtype=float)
    host=np.asarray([int(r["host"]=="As") for r in rows],dtype=float)
    cage=np.asarray([r["cage"] for r in rows],dtype=float)
    terms=[np.ones(len(rows)),year,host,cage]
    X0=np.column_stack(terms)
    X1=np.column_stack(terms+[host*cage])
    reduced=sm.GLM(y,X0,family=sm.families.Binomial()).fit()
    full=sm.GLM(y,X1,family=sm.families.Binomial()).fit()
    lr=max(0,2*(full.llf-reduced.llf))
    return {"LR_statistic":float(lr),"chi2_df":1,"LR_p_two_sided":float(chi2.sf(lr,1)),
            "interaction_log_odds_coef":float(full.params[-1]),
            "interaction_se":float(full.bse[-1]),"full_converged":bool(full.converged),
            "model":"binomial GLM y ~ year2019 + As + cage + As:cage"}

def test(source_dir,output):
    data,audit=read_data(source_dir)
    if len(data)<200:raise RuntimeError("common garden source unexpectedly small")
    surv=empirical(data)
    results={
      "schema":"chocho_monarch_H3_host_predator_exposure_adult_survival_v0.1",
      "status":"EXPLORATORY_INDEPENDENT_SOURCE_REANALYSIS_NOT_A_NEW_EXPERIMENT",
      "protocol":"docs/exploratory/MONARCH_H3_HOST_ENEMY_SURVIVAL_PROTOCOL_V01.json",
      "source":"Diethelm et al 2026 Dryad DOI 10.5061/dryad.51c59zwgx",
      "quality_context":"Both Asclepias speciosa and Asclepias fascicularis are H3 high-performance Greenstein 2022 classifications",
      "source_audit":audit,"survival":surv,
      "GLM_interaction":glm_test(data),
      "interpretation_boundary":"The study already reports host-by-treatment growth/development outcomes. Any survival result is a reanalysis of published common-garden data. Predator-exclusion structures may change humidity/temperature; cage is not a pure predation manipulation. This does not establish survival outcomes in any introduced plant region or butterfly geographic colonization."
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(results,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(results,indent=2),flush=True)

def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest="cmd",required=True)
    down=sub.add_parser("download");down.add_argument("--directory",type=Path,required=True)
    testp=sub.add_parser("test");testp.add_argument("--directory",type=Path,required=True)
    testp.add_argument("--output-json",type=Path,required=True)
    args=p.parse_args()
    if args.cmd=="download":download(args.directory)
    else:test(args.directory,args.output_json)
if __name__=="__main__":main()

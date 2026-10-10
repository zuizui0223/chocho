#!/usr/bin/env python3
"""MD5-identical source structure audit of larval instar-count time series.

No newly fitted biological models. Reports whether the 2012 30-cage series can
describe the onset of high-demand instars, without pretending to know food
accessibility at those dates.
"""
import argparse,csv,hashlib,io,json
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

URL="https://ndownloader.figshare.com/files/41146991"
MD5="0e9ff924e4ae38cb78529ec540887c42"
COLUMNS=["plot.no","treat.no","s.density","a.density","date",
         "total.s","total.a","day","s1","s2","s3","s4","s5",
         "spp","sp","cumul.sp","a1","a2","a3","a4","a5",
         "app","ap","cumul.ap"]

def num(x):
    t=str(x or "").strip()
    if not t:return None
    try:return float(t)
    except ValueError:return None

def analyze(raw):
    if hashlib.md5(raw).hexdigest()!=MD5:raise ValueError("original source checksum mismatch")
    lines=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(lines)!=690 or list(lines[0])!=COLUMNS:raise ValueError("source schema changed")
    cages=defaultdict(list)
    nonnumeric=Counter()
    for r in lines:
        cages[str(r["plot.no"])].append(r)
        for k in COLUMNS:
            if k not in ("date",) and num(r[k]) is None:
                nonnumeric[k]+=1
    if len(cages)!=30:raise ValueError("unexpected number of cage IDs")
    samples=[]
    early_late_support=Counter()
    for cage,rr in sorted(cages.items(),key=lambda a:int(a[0])):
        rr=sorted(rr,key=lambda r:num(r["day"]))
        if len(rr)!=23 or len({r["day"] for r in rr})!=23:
            raise ValueError("unexpected repeated days in cage "+cage)
        first=rr[0];last=rr[-1]
        late_count=lambda r:sum(num(r[f"a{j}"]) or 0 for j in (3,4,5))
        early_count=lambda r:sum(num(r[f"a{j}"]) or 0 for j in (1,2))
        early_late_support["cages_with_A_initial"]+=int(num(first["a.density"])>0)
        early_late_support["cages_with_late_A_observed"]+=int(any(late_count(r)>0 for r in rr))
        early_late_support["cages_with_early_A_observed"]+=int(any(early_count(r)>0 for r in rr))
        entry=[num(r["day"]) for r in rr if late_count(r)>0]
        first_pup=[num(r["day"]) for r in rr if (num(r["cumul.ap"]) or 0)>0]
        samples.append({
            "cage":cage,
            "initial_S":num(first["s.density"]),
            "initial_A":num(first["a.density"]),
            "survey_days":[num(r["day"]) for r in rr],
            "first_late_A_instar_observed_day":min(entry) if entry else None,
            "first_A_pupation_observed_day":min(first_pup) if first_pup else None,
            "last_day":num(last["day"]),
            "last_Atrophaneura_total":num(last["total.a"]),
            "last_Atrophaneura_cumulative_pupation":num(last["cumul.ap"]),
            "source_first_A_stage_columns":{k:first[k] for k in ("a1","a2","a3","a4","a5","app","ap","cumul.ap","total.a")},
            "source_last_A_stage_columns":{k:last[k] for k in ("a1","a2","a3","a4","a5","app","ap","cumul.ap","total.a")}
        })
    late_onsets=Counter(str(row["first_late_A_instar_observed_day"])
                        for row in samples if row["first_late_A_instar_observed_day"] is not None)
    pup_onsets=Counter(str(row["first_A_pupation_observed_day"])
                       for row in samples if row["first_A_pupation_observed_day"] is not None)
    return {
        "schema":"chocho_kyoto_stage_demand_source_only_v01",
        "original_DOI":"10.1002/ece3.10164",
        "source_url":URL,
        "source_md5":MD5,
        "n_source_cages":len(cages),
        "n_original_cage_date_rows":len(lines),
        "non_numeric_fields_excluding_date":dict(nonnumeric),
        "instar_support":dict(early_late_support),
        "first_observed_3rd_to_5th_instar_day_distribution":dict(sorted(late_onsets.items(),key=lambda kv:float(kv[0]))),
        "first_observed_pupation_day_distribution":dict(sorted(pup_onsets.items(),key=lambda kv:float(kv[0]))),
        "per_cage_stage_onsets":samples,
        "resource_access_or_quality_observed_at_stage":False,
        "new_butterfly_survival_effect_estimated":False,
        "GEB_PR38_unchanged":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--receipt",required=True,type=Path)
    a=ap.parse_args()
    with urlopen(Request(URL,headers={"User-Agent":"chocho-stages-identity-audit/1.0"}),timeout=35) as response:
        raw=response.read(500000)
    out=analyze(raw)
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(out,indent=2)+"\n")
    # Only summarize source support; avoid dumping potentially misinterpreted stage-by-stage contrasts.
    print(json.dumps({k:v for k,v in out.items() if k!="per_cage_stage_onsets"},indent=2))
    print("Receipt contains stage onset per original cage; no plant-availability measurements in old experiment.")

if __name__=="__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from query_gbif_event_host_first_records import resolve, geometry_by_code, query_years

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--event-csv',type=Path,required=True)
    ap.add_argument('--level3-geojson',type=Path,required=True)
    ap.add_argument('--output-csv',type=Path,required=True)
    ap.add_argument('--output-json',type=Path,required=True)
    a=ap.parse_args()
    geom=geometry_by_code(a.level3_geojson)
    rows=list(csv.DictReader(a.event_csv.open(newline='',encoding='utf-8')))
    tax={}; out=[]
    for i,r in enumerate(rows,1):
        sp=r['species'].strip(); code=r['wgsrpd3_code'].strip()
        if sp not in tax: tax[sp]=resolve(sp)
        m=tax[sp]; error=''; q={'count':0,'years':[]}
        if not m['accepted']: error='REJECTED_TAXON_MATCH'
        elif code not in geom: error='MISSING_WGSRPD3_GEOMETRY'
        else:
            try:q=query_years(m['usage_key'],geom[code])
            except Exception as e:error=str(e)
        years=[y for y,n in q.get('years',[]) if 1750<=y<=2026 and n>0]
        first=min(years) if years else None
        snap=int(r['snapshot_first_year'])
        out.append({**r,'gbif_usage_key':m.get('usage_key'),'gbif_full_coordinate_records_in_region':q.get('count',0),'gbif_full_first_record_year':first,'snapshot_first_minus_full_first_year':None if first is None else snap-first,'has_pre2010_record':None if first is None else int(first<2010),'is_true_post2017_first_gbif_record':None if first is None else int(2018<=first<=2025),'error':error})
        print(i,len(rows),sp,code,'full_first',first,'snapshot',snap,'error',error,flush=True)
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_csv.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    by_group={}
    for group in sorted(set(r['group'] for r in out)):
        rr=[r for r in out if r['group']==group]
        known=[r for r in rr if r['gbif_full_first_record_year'] not in ('',None)]
        by_group[group]={'snapshot_events':len(rr),'full_history_resolved':len(known),'pre2010_records':sum(int(r['has_pre2010_record']) for r in known),'true_post2017_first_records':sum(int(r['is_true_post2017_first_gbif_record']) for r in known),'baseline_2010_2017_first_records':sum(2010<=int(r['gbif_full_first_record_year'])<=2017 for r in known)}
    payload={'schema':'chocho_stageA_event_full_history_audit_v0.1','status':'POSTHOC_FULL_GBIF_HISTORY_EVENT_AUDIT','events':len(out),'by_group':by_group,'claim_boundary':['This audits only cells counted as 2018-2025 events in the capped occurrence snapshot. It does not correct candidate denominators.','Full-history first GBIF record is still a detection date, not colonization date.']}
    a.output_json.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2))
if __name__=='__main__':main()

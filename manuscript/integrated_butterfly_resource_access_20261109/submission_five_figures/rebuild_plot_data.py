"""Rebuild display data from immutable GitHub blobs; no raw ecological refit."""
import argparse,csv,hashlib,io,json,subprocess,urllib.request
from pathlib import Path
P=Path(__file__).resolve().parent
M=json.loads((P/'source_manifest.json').read_text())

def canonical(x):
    if x is None:return 'null'
    if isinstance(x,(int,float)):return f'{x:.9f}'
    if isinstance(x,list):return '['+','.join(map(canonical,x))+']'
    return '{'+','.join(k+':'+canonical(x[k]) for k in sorted(x))+'}'
def fnv(s):
    h=0x811c9dc5
    for ch in s:h=((h^ord(ch))*0x01000193)&0xffffffff
    return f'{h:08x}'
def normalize(d):
    d=dict(d);d.pop('source_commit',None)
    for k in ('scatter','observations'):
        if isinstance(d[k],str):d[k]=[[float(x) for x in l.split(',')] for l in d[k].split(';')]
    return d
def check(d):
    d=normalize(d)
    assert set(d)==set(M['fingerprints'])
    for k,v in d.items():
        got=fnv(canonical(v));assert got==M['fingerprints'][k], (k,got,M['fingerprints'][k])
    print('PASS: 17 source-derived display fields (FNV-1a over 9-decimal canonical numbers)')
    return d

def read_sources(args):
    result={}
    for s in M['sources']:
        if args.git_repo:
            b=subprocess.check_output(['git','-C',args.git_repo,'show',s['ref']+':'+s['path']])
        elif args.source_dir:
            b=(Path(args.source_dir)/s['path']).read_bytes()
        else:
            url='https://raw.githubusercontent.com/'+M['repository']+'/'+s['ref']+'/'+s['path']
            with urllib.request.urlopen(url,timeout=45) as r:b=r.read()
        sha=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        if sha!=s['blob_sha']:raise ValueError(f"Git blob mismatch: {s['key']} {sha}")
        result[s['key']]=b.decode('utf-8')
    return result

def rebuild(R):
    rows=lambda k:list(csv.DictReader(io.StringIO(R[k])))
    m=rows('metrics');g=rows('gaps');i=[r for r in g if r['group']=='island'];main=[r for r in g if r['group']=='mainland'];mi={r['species']:r for r in i}
    s=sorted(m,key=lambda r:(int(r['native_resource_units']),r['species']))
    assert len(s)==len(mi)==239 and set(mi)=={r['species'] for r in s}
    j=json.loads(R['v04'])['figures'];d=json.loads(R['v07']);q=json.loads(R['isolation'])['models'];p=json.loads(R['port'])['report']['models'];c=json.loads(R['community'])['primary'];f=j['fig4'];mc=j['fig2']['marine_circuit']
    total=lambda a,k:sum(int(r[k]) for r in a)
    est=lambda x,ci:[x,*ci]
    groups=[]
    for l,h in [(0,60),(60,120),(120,180),(180,239)]:
        z=[mi[r['species']] for r in s[l:h]];groups.append([len(z),total(z,'added'),total(z,'empty')])
    B,C,E=f['panel_B_224_butterfly_source'],f['panel_C_eight_unseen_worlds'],f['panel_D_south_america']
    out={'source_commit':M['sources'][0]['ref'],
         'scatter':[[int(r['native_resource_units']),int(mi[r['species']]['added'])] for r in s],
         'observations':[r[1:4] for r in d['observations']],
         'fit':{k:d['fits'][k] for k in ['mean_logdistance','sd_logdistance','min_km','max_km','observed','potential']},
         'baseline':[[total(a,k) for k in ['cells','native','added']] for a in [i,main]],
         'port':[est(j['fig2']['port_coefficient'][k]['estimate'],j['fig2']['port_coefficient'][k]['descriptive_HC3_95']) for k in ['introduced','native','introduced_minus_native']],
         'jaccard':[j['fig2']['whole_flora_jaccard'][k+'_median'] for k in ['introduced','native']],
         'count_predictions':[[q['marginal_count_predictions'][f'km{v}'][k] for v in [100,1000,4000]] for k in ['native','contemporary']],
         'shares':[q['adjusted_shares'][f'km{v}'] for v in [100,1000,4000]],
         'hierarchy':[r[1:] for r in d['hierarchy']],
         'paired_shift':est(q['paired_native_current_shift']['distance_beta'],q['paired_native_current_shift']['distance_nine_cluster_95']),
         'top_k':[[int(k),v] for k,v in sorted(j['fig3']['panelA_plants']['top_k_cumulative_share'].items(),key=lambda kv:int(kv[0]))]+[[j['fig3']['panelA_plants']['total_introduced_host_contributors'],1]],
         'gap_groups':groups,
         'port_ratio':[est(p[k]['beta'],p[k]['cluster_t_95']) for k in ['all126_basic','all126_plus_geography','original98_effort','common98_geography_without_effort','common98_geography_plus_effort']],
         'marine':[est(mc[t][k+'_MSE_gain_percent'],mc[t][k+'_descriptive_interval']) for t in ['source_original_one_island_heldout','strict_both_endpoints_unseen'] for k in ['introduced','native']],
         'community':[est(c['H_given_other_MSE_improvement_percent'],[round(x,4) for x in c['bootstrap_percentile_95_unstrat']]),est(B['dataset_days_MSE_improvement_pct_median'],B['dataset_days_MC95']),est(B['record_count_MSE_improvement_pct_median'],B['record_count_MC95'])],
         'transfer':[[C['original_H_after_O_world_equal_pct'],None,None],est(C['dataset_day_median'],C['dataset_day_MC95']),est(C['record_count_median'],C['record_count_MC95'])],
         'south_america':[est(E[k+'_'+t+'_median'],E[k+'_'+t+'_MC95']) for k in ['dataset_day','record_count'] for t in ['within','cross']]}
    for k in ('scatter','observations'):out[k]=';'.join(','.join(str(x) for x in row) for row in out[k])
    return out

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--git-repo');g.add_argument('--source-dir');g.add_argument('--download',action='store_true');p.add_argument('--output',default=str(P/'plot_data.json'));a=p.parse_args()
    out=rebuild(read_sources(a));check(out)
    target=Path(a.output)
    if target.exists():
        previous=json.loads(target.read_text())
        assert normalize(out)==normalize(previous),'Refusing to overwrite different display input'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: pinned public source data and all 17 fingerprints; wrote',target)
if __name__=='__main__':main()

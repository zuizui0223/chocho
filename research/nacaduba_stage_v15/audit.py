"""Source-locked, no-dependencies flowering-coverage reproduction; not insect survival."""
from pathlib import Path
from collections import defaultdict
from math import sqrt
import csv
import json

ROOT=Path(__file__).resolve().parent
YEARS=(2008,2009,2010)
MONTHS=('Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar')
PARTS=('early','mid','late')

def load(filename):
    with (ROOT/'data'/filename).open(newline='',encoding='utf8') as fh:
        return list(csv.DictReader(fh))

def analyze():
    clones=load('clones_kato2012_table4_table5.csv')
    ramets=load('ramets_kato2012_table7.csv')
    assert len(clones)==23 and sum(int(c['ramets']) for c in clones)==39
    assert len(ramets)==36 and len({(r['clone'],r['ramet']) for r in ramets})==36
    expected_means=[65.21739130434783,73.91304347826087,60.869565217391305]
    expected_se=[9.6555205536,8.8183049476,9.4076128343]
    expected_coverage={2008:(9,2,8,19),2009:(6,4,7,17),2010:(8,3,10,21)}
    expected_early={2008:(11,20),2009:(15,21),2010:(12,20)}
    result=[]
    events=0
    for idx,year in enumerate(YEARS):
        vals=[int(c[f'flower_pct_{year}']) for c in clones]
        mu=sum(vals)/len(vals)
        se=sqrt(sum((x-mu)**2 for x in vals)/(len(vals)*(len(vals)-1)))
        assert abs(mu-expected_means[idx])<1e-9
        assert abs(se-expected_se[idx])<1e-8
        by_group={'g':set(),'other':set()}
        starts={'g':[],'other':[]}
        by_clone=defaultdict(int)
        group_flow={'g':0,'other':0}
        durations=[]
        for row in ramets:
            b=row[f'begin_{year}']
            d=row[f'duration_{year}']
            assert bool(b)==bool(d), (row,year)
            if not b:continue
            part,month=b.split()
            start=MONTHS.index(month)*3+PARTS.index(part)
            duration=int(d)
            assert duration>=1 and start+duration<=36
            group='g' if row['seed_source']=='g' else 'other'
            group_flow[group]+=1
            by_clone[row['clone']]+=1
            starts[group].append(start)
            by_group[group].update(range(start,start+duration))
            durations.append(duration)
            events+=1
        for c in clones:
            assert by_clone[c['clone']]*100==int(c['flower_pct_'+str(year)])*int(c['ramets'])
        g=by_group['g'];o=by_group['other']
        counts=(len(g-o),len(o-g),len(g&o),len(g|o))
        assert counts==expected_coverage[year],(year,counts)
        assert (min(starts['g']),min(starts['other']))==expected_early[year]
        assert group_flow['g']==6
        assert group_flow['other']=={2008:19,2009:23,2010:19}[year]
        assert max(durations)=={2008:17,2009:13,2010:17}[year]
        result.append({
            'year':year,'mean_clone_flowering_percent':mu,
            'source_g_flowering_ramets':group_flow['g'],
            'other_flowering_ramets':group_flow['other'],
            'thai_only_10day_bins':counts[0],'other_only_10day_bins':counts[1],
            'shared_10day_bins':counts[2],'combined_10day_bins':counts[3],
            'thai_earliest_10day_bin':min(starts['g']),
            'other_earliest_10day_bin':min(starts['other']),
            'table7_max_duration_10day_bins':max(durations)
        })
    assert events==79
    assert result[1]['table7_max_duration_10day_bins']==13
    return {
        'provenance':'Kato et al. (2012), Silvae Genetica 61(6):228-236, doi:10.1515/sg-2012-0029, Tables 4/5/7',
        'study_site':'Iriomote cultivated clones, not Singapore',
        'n_clones':23,'n_ramets':39,'n_flowering_ramet_year_events':events,
        'years':result,
        'source_discrepancy':{
            '2009_max_in_paper_prose_10day_bins':8,
            '2009_max_in_table7_10day_bins':13,
            'published_table_used_with_transparency':True
        },
        'singapore_daily_buds':None,
        'butterfly_flowering_mismatch_effect':None,
        'biological_colonization_lag':None
    }

if __name__=='__main__':
    print(json.dumps(analyze(),indent=2,ensure_ascii=False))

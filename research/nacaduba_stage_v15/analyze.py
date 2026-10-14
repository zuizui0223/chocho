"""Source-bounded, dependency-free phenology check; not a butterfly synchrony estimator."""
import csv
import json
from pathlib import Path
from statistics import mean,median,stdev
FILE=Path(__file__).with_name('published_monthly_flowering.csv')
def viable_start_days(seq, n):
    if len(seq)<n:return 0
    return sum(all(seq[i:i+n]) for i in range(len(seq)-n+1))
def analyze():
    with FILE.open(newline='') as f:rows=list(csv.DictReader(f))
    assert len(rows)==24
    missing=[x['month'] for x in rows if not x['acacia_auriculiformis_flowering_intensity_percent']]
    assert missing==['1989-02']
    series={x['month']:float(x['acacia_auriculiformis_flowering_intensity_percent']) for x in rows if x['acacia_auriculiformis_flowering_intensity_percent']}
    v=list(series.values())
    assert len(v)==23 and series['1988-01']==65.6 and series['1988-06']==85.6
    assert series['1989-01']==0.6 and series['1989-06']==74.4
    assert min(v)==0.6 and max(v)==85.6
    count={'little':sum(x<33.3 for x in v),'moderate':sum(33.3<=x<=66.6 for x in v),'heavy':sum(x>66.6 for x in v)}
    a=viable_start_days([1]*15+[0]*15,10)
    b=viable_start_days([1,0]*15,10)
    assert (a,b)==(6,0)
    return {'reported_months':len(v),'omitted_month':missing,'mean':mean(v),'median':median(v),
            'max':{'month':'1988-06','intensity_percent':max(v)},'min':{'month':'1989-01','intensity_percent':min(v)},
            'max_min_fold':max(v)/min(v),'cv':stdev(v)/mean(v),
            'month_categories':count,'same_month_january_1988_vs_1989':[series['1988-01'],series['1989-01']],
            'counterexample_constructed_not_observed':{'same_month_binary_bloom_days':15,'larval_days_example':10,
                 'continuous_block_valid_start_days':a,'alternate_day_valid_start_days':b},
            'joint_singapore_butterfly_plant_pheno_effect':None,
            'measurement':'normalized stand flowering crown intensity /180, not % trees flowering, not buds/day',
            'sites':'Peninsular Malaysia plants 1987-1989, Singapore butterfly 2010 unpaired'}
if __name__=='__main__':print(json.dumps(analyze(),indent=2))

"""Offline prior cutoff, observational host placebo and parameter controls."""
import sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_uk_bsbi_historical_host_increment_20261029 import (
    placebo_ids,raw_design_with_era,SEED,TRAIN,TEST,N_PLACEBOS)
from audit_uk_ltmn_butterfly_detection_panel_20261028 import CROSSWALK

def test_temporal_cutoff_fixed_to_after_atlas_1986():
    assert TRAIN==(1987,1995) and TEST==(2005,2021)
    assert N_PLACEBOS==19

def test_host_count_not_no_plant_absence():
    rows=[('S1','Aglais urticae','early',1,[0.,1.],[18.,1.]),
          ('S2','Aglais urticae','late',0,[0.,1.],[19.,1.])]
    count={'S1':300,'S2':500}
    h={'S1':'SP12','S2':'SP13'}
    doc={'SP12':{'h1','other'},'SP13':{'other'}}
    source={'Aglais urticae':{'h1'}}
    x=raw_design_with_era(rows,set(count),{'Aglais urticae'},doc,h,source,count)
    assert x["host"][0][4][-1]==math.log1p(1)
    assert x["host"][1][4][-1]==0.0
    assert x["weather"][0][4][-1]==math.log1p(300)

def test_prevalence_matched_nonhost_not_original_host():
    prev={'h1':2,'h2':1,'n1':2,'n2':2,'n3':1,'n4':1}
    a=placebo_ids({'h1','h2'},prev,0,'Species A')
    assert len(a)==2 and not a&{'h1','h2'}
    assert sorted(prev[x] for x in a)==[1,2]
    assert placebo_ids({'h1','h2'},prev,0,'Species A')==a

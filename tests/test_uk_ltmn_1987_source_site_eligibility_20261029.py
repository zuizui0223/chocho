"""Survey effort-only site eligiblity never opens butterfly outcome."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_uk_ltmn_1987_source_site_eligibility_20261029 import select

def v(s,y,day):
    return dict(UKBMS_SITENO=s,YEAR=y,MONTH=6,DAY=day,START_TIME=.5)

def test_three_early_years_and_late_years_required():
    rows=[]
    for y in (1987,1990,1995,2006,2011,2016):
        for d in (1,2,3,4,5):rows.append(v(1,y,d))
    good,per,years,seen=select(rows,{"1":(123000,240000,("", ""))})
    assert good=={"1"}
    assert seen==30

def test_no_1981_to_1986_training_after_1986_cutoff():
    rows=[]
    for y in (1981,1983,1985,2005,2007,2009):
        for d in range(1,6):rows.append(v(2,y,d))
    sites,*_=select(rows,{"2":(123000,240000,("",""))})
    assert sites==set()

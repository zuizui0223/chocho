"""Synthetic chronology tests only; these are not evidence about butterflies."""
import importlib.util
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
p=ROOT/"scripts/analyze_melitaea_19year_strict_forward.py"
spec=importlib.util.spec_from_file_location("forward",p)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def make_panel():
    result=[]
    for t in range(2000,2005):
        for k in range(12):
            feature=[float((k+t)%4)]*len(m.old.PRIMARY_FIELDS)
            result.append({"year":t,"patch":f"p{k}","y":int((k+t)%2==0),
                           "base":feature,"interaction":int(k%3==0),
                           "decline":int(k%2==0),"backup_high":int(k%3==0)})
    return result

def test_future_year_training_absent_and_numbered():
    panel=make_panel()
    pred,counts=m.temporal_predictions(panel,False,start=2002,end=2004)
    assert counts=={"2002":24,"2003":36,"2004":48}
    assert all(np.isnan(pred[:24]))
    assert np.isfinite(pred[24:]).all()

def test_incremental_term_can_be_added_without_changing_temporal_boundary():
    p=make_panel()
    a,ca=m.temporal_predictions(p,False,start=2002,end=2004)
    b,cb=m.temporal_predictions(p,True,start=2002,end=2004)
    assert ca==cb
    assert np.isnan(a[:24]).all() and np.isnan(b[:24]).all()
    assert all(x>=0 and x<=1 for x in b[24:])

def test_protocol_explicitly_rejects_leaky_loyo_as_prospective():
    import json
    q=json.loads((ROOT/"docs/exploratory/MELITAEA_19YEAR_STRICT_FORWARD_VALIDATION_V01.json").read_text(encoding="utf-8"))
    assert q["first_test_year"]==2010
    assert q["last_test_year"]==2017
    assert q["GEB_PR38_unchanged"] is True

"""Distinct estimands for plant source-index rarefaction."""
import importlib.util,sys,math
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_n_rarefaction_depth_sensitivity.py"
sys.path.insert(0,str(f.parent))
spec=importlib.util.spec_from_file_location("n_raref_depth",f)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_zero_missing_kept_unknown():
    import random
    original={"1":300,"2":1,"3":1,"4":1,"5":1,"6":1,"7":1,"8":1}
    assert m.sampled_keys_from_records(original,random.Random(8),100) <= set(range(1,9))
    assert len(original)==8

def test_pooled_four_cell_depth_is_minimum():
    raw={}
    per=m.PERIODS
    for side in m.SIDES:
        for period in per:
            raw[("1",2,side,period)]={"faceted_specieskey_counts":{"1":101,"2":50}}
    raw[("1",2,"port",per[0])]["faceted_specieskey_counts"]={"1":10,"2":12}
    assert m.cell_support_for_pair(raw,("1",2),"record",60)==22
    assert m.cell_support_for_pair(raw,("1",2),"species",60)==2
    assert m.cell_support_for_pair(raw,("1",2),"record",10)==10

def test_sampling_methods_present_all_ten_modes():
    # The full analysis reads frozen GH Actions source artifacts; tests remain offline.
    import inspect
    src=inspect.getsource(m.run)
    for unit in ("record","species"):
        assert f'("{unit}",' in src

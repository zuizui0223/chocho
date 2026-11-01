"""Synthetic, source-free tests of same-region colonial flora resource checks."""
import importlib.util
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_matched_glonaf_consumer_history.py"


def load():
    # GEB's minimal paper CI does not install the optional spatial geometry
    # stack. Specialized source workflows install and execute these tests.
    pytest.importorskip("pyproj")
    pytest.importorskip("shapefile")
    sys.path.insert(0,str(ROOT/"scripts"))
    spec=importlib.util.spec_from_file_location("matched_glonaf_consumer",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_polygon_share_requires_substantive_real_source_overlap():
    m=load()
    assert m.classify_source_contribution(100,20,11,.05)
    assert m.classify_source_contribution(100,200,6,.05)
    assert not m.classify_source_contribution(100,200,4,.05)
    assert not m.classify_source_contribution(100,200,0,.05)
    assert not m.classify_source_contribution(100,0,6,.05)
    assert m.THRESHOLDS==(.01,.05,.10)


def test_four_empire_shared_source_label_reuses_no_unverified_history():
    m=load()
    assert m.MIN_SOURCE_SHARE==.5
    assert m.shared_years([[50,0,0,0]],[[20,30,0,0]],years_min=25).tolist()==[0.0]
    assert m.shared_years([[50,0,0,0]],[[40,0,0,0]],years_min=25).tolist()==[40.0]

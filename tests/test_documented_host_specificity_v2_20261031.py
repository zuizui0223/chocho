"""Source-only tests: do not claim high-prevalence pseudo-host matching where none exists."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_documented_host_specificity_v2.py"


def load():
    pytest.importorskip("pyproj")
    pytest.importorskip("shapefile")
    sys.path.insert(0,str(ROOT/"scripts"))
    spec=importlib.util.spec_from_file_location("documented_host_specificity_v2",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_entire_unsupported_frequency_class_is_excluded_without_choosing_species():
    m=load()
    assert m.BINS==(1,2,5,10,25,50,100,200,285)
    # Bin index 7 = 101–200 (NOT 201–284).
    freq=np.asarray([150,140,180,155,60,80],int)
    bins={7:[0,1,2,3],6:[4,5]}
    documented=(1<<0)|(1<<1)|(1<<2)|(1<<4)
    selected,support=m.overlap_supported_mask(documented,freq,bins)
    assert selected==(1<<4)
    assert support["eligible_original_taxa"]==4
    assert support["retained_comparable_taxa"]==1
    assert support["excluded_taxa_due_to_no_pseudohost_support"]==3
    assert support["supported_frequency_bins"]==[6]
    assert support["excluded_frequency_bins"]==[7]
    assert support["source_frequency_bins"]["7"]["prevalence_upper_bound"]==200


def test_no_valid_alternative_fails_and_no_species_is_silently_reused():
    m=load()
    freq=np.asarray([150,140,180],int)
    raw=(1<<0)|(1<<1)
    with pytest.raises(RuntimeError,match="No source-supported"):
        m.overlap_supported_mask(raw,freq,{7:[0,1,2]})
    with pytest.raises(RuntimeError,match="No source-active"):
        m.overlap_supported_mask(1<<0,np.asarray([0,10],int),{3:[1]})


def test_method_lock_declares_posthoc_source_support_change():
    m=load()
    import json
    d=json.loads(m.LOCK.read_text(encoding="utf8"))
    assert d["status"].startswith("EXPLORATORY_AMENDMENT")
    assert d["min_supported_region_pairs"]>=20000
    assert "101-200" in d["frequency_bins_source"]

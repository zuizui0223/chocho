import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from run_island17_historical_PLSCIQ_butterfly_LO2 import (
    WPI17,WPI13,RECEIPT,ROOT,PAIRS,ROSTER,TARGETS,MODELS,COMPARISONS,SCOPES,
    original_source,old_port_index_source,designs
)

def test_old_UNCTAD_2026_10_30_exact_edition_and_geoverified_source_population():
    if not WPI13.is_file() or not WPI17.is_file() or not RECEIPT.is_file():
        pytest.skip("Original source-only PortWatch/UNCTAD inputs omitted from anonymous review package")
    lock,q,ids,counts=original_source()
    thirty,twentyfour,twentytwo=old_port_index_source(ids)
    assert len(ids)==51 and len(q)==1275
    assert sum(counts.values())==194
    assert len(thirty)==30 and len(twentyfour)==24 and len(twentytwo)==22
    assert twentytwo.original_island_id.nunique()==17
    assert all(twentytwo.available_2006_2010_quarters.astype(int)>=16)
    assert len(MODELS)==4 and len(TARGETS)==2 and len(COMPARISONS)==5 and len(SCOPES)==2

def test_history_source_target_predictor_order_and_group9_scope():
    if not WPI13.is_file() or not RECEIPT.is_file():
        pytest.skip("Source original 13 port records are missing from anonymous package")
    lock,q,ids,counts=original_source()
    _,_,qualified=old_port_index_source(ids)
    designs_by_target,names,median,receipt,source=designs(q,ids,counts,qualified)
    assert len(median)==17 and "9" in median
    mask=q.island_a.isin(median)&q.island_b.isin(median)
    without9=mask&~q.island_a.eq("9")&~q.island_b.eq("9")
    assert int(mask.sum())==136 and int(without9.sum())==120
    for target in TARGETS:
        assert set(names[target])==set(MODELS)
        assert all("host" not in v.lower() and "ais" not in v.lower()
                   for model in MODELS for v in names[target][model])
        assert "previous_2001_2010_butterfly_Jaccard" in names["new_butterfly_jaccard"]["BASE"]
        assert "previous_2001_2010_butterfly_Jaccard" not in names["earlier_butterfly_jaccard"]["BASE"]
        for model in MODELS:
            x=designs_by_target[target][model]
            assert x.shape[0]==1275 and np.isfinite(x[mask]).all()

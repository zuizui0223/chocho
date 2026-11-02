"""Offline synthetic checks for global urban source estimator."""
import pytest
pytest.importorskip("rasterio")
pytest.importorskip("statsmodels")
import importlib.util,sys,numpy as np,pandas as pd
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_global_hydrowaste_worldcover_adjustment.py"
sys.path.insert(0,str(f.parent))
sp=importlib.util.spec_from_file_location("urban_global",f)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_COG_cache_has_explicit_small_file_handle_limit():
    cache=m.BoundedWorldCoverTiles(max_open=7)
    assert cache.max_open==7 and len(cache.cache)==0
    cache.close()

def test_urban_model_defends_underpowered_support():
    n=90
    df=pd.DataFrame({"port_id":np.arange(n),"in_within_country_comparison":[True]*n,
          "quality":["COMPLETE"]*n,
          "country":np.repeat(["JP","GB","US"],30),
          "WWTP_count_25km":np.tile(np.arange(30)%4,3)})
    result=m.country_model(df,urbancovariates=[])
    assert result["status"]=="INSUFFICIENT_URBAN_SOURCE_SUPPORT"

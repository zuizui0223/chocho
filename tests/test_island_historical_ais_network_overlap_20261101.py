"""Do not confuse country metropole-to-former colony with colony–colony or no route."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/audit_ais_port_colonial_island_overlap.py"


def mod():
    spec=importlib.util.spec_from_file_location("island_hist_ais",SRC)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def p(code,gb=0,es=0,pt=0,nl=0,eligible="1"):
    return {"eligible_for_history_test":eligible,
            "current_WPI_country_codes":code,
            "years_GB":str(gb),"years_ESP":str(es),
            "years_PRT":str(pt),"years_NL":str(nl)}


def test_metropolitan_britain_to_cyprus_not_mislabeled_shared_colony():
    m=mod()
    britain=p("GB")
    cyprus=p("CY",gb=73)
    common,metropole,which=m.classify_pair(britain,cyprus,{})
    assert common is False
    assert metropole is True
    assert which==["GB"]
    # Two separately occupied former colonies share historical GB exposure
    # but neither endpoint is a metropolitan GB country.
    a,b=p("NZ",gb=63),p("CY",gb=73)
    assert m.classify_pair(a,b,{})==(True,False,[])


def test_unsupported_historic_region_is_unknown_not_zero():
    m=mod()
    u=p("FR",eligible="0")
    s=p("CY",gb=73)
    assert m.classify_pair(u,s,{})==(None,None,[])
    assert m.connected_components([("1","2"),("2","3"),("9","10")])==[
        ["1","2","3"],["9","10"]
    ]
    with pytest.raises(ValueError,match="predeclared"):
        m.analyse(30)

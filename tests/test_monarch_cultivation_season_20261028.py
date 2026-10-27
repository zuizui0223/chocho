"""Source-independent scientific controls for host-by-state-year season contrast."""
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_monarch_cultivation_season_20261028 import (
    exact_binomial, canonical_host, table_strata, mh_odds_ratio, summary, group_bootstrap,
)


def row(host, state, year, cult, winter, *, origin="native", confidence="high",
        hex_id="A"):
    return {"host":host,"state":state,"year":year,
            "cultivated":cult,"winter":winter,
            "month":1 if winter else 6,
            "origin":origin,"confidence":confidence,"hex":hex_id}


def test_native_host_same_host_state_year_matching():
    data=(
        [row("Asclepias fascicularis","CA",2020,"yes",True)]*6+
        [row("Asclepias fascicularis","CA",2020,"yes",False)]*4+
        [row("Asclepias fascicularis","CA",2020,"no",True)]*2+
        [row("Asclepias fascicularis","CA",2020,"no",False)]*8+
        [row("Asclepias fascicularis","CA",2021,"no",True)]*3+
        [row("Asclepias curassavica","CA",2020,"yes",True,origin="non-native")]*10
    )
    strata, qc=table_strata(data)
    assert len(strata)==1 and qc["source_qualified_records"]==23
    a=strata[0]
    assert (a["a"],a["b"],a["c"],a["d"]) == (6,4,2,8)
    assert mh_odds_ratio(strata)==pytest.approx(6.0)
    out=summary(strata,keys=("host","state","year"),iterations=99)
    assert out["mh_recording_composition_odds_ratio"]==pytest.approx(6.0)
    assert out["source_information_gate"]=="FAIL_PREDECLARED_SOURCE_SUPPORT"


def test_unknown_and_likely_labels_not_counted_as_strict():
    base=[row("Asclepias incarnata","MA",2019,"likely",True),
          row("Asclepias incarnata","MA",2019,"unlikely",False),
          row("Asclepias incarnata","MA",2019,"maybe",False)]
    strict,_=table_strata(base)
    expanded,_=table_strata(base,label_mode="expanded")
    assert strict==[]
    assert len(expanded)==1
    assert expanded[0]["a"]==1 and expanded[0]["d"]==1


def test_no_pseudocount_or_fake_informative_stratum():
    rows=[row("Asclepias syriaca","OH",2020,"yes",False)]*4 + [
        row("Asclepias syriaca","OH",2020,"no",False)]*6
    strata,_=table_strata(rows)
    assert not strata and mh_odds_ratio(strata) is None
    x=[{"a":2,"b":0,"c":0,"d":5,"n":7,"state":"OH"}]
    assert mh_odds_ratio(x) is None
    assert group_bootstrap([],iterations=10)["valid_draws"]==0


def test_origin_secondary_differs_from_same_species_cultivation():
    source=[
        row("Asclepias syriaca","TX",2020,"no",True),
        row("Asclepias syriaca","TX",2020,"no",False),
        row("Asclepias curassavica","TX",2020,"no",True,origin="non-native"),
        row("Asclepias curassavica","TX",2020,"no",False,origin="non-native"),
    ]
    a,_=table_strata(source,comparison="origin",keys=("state","year"))
    assert len(a)==1 and a[0]["a"]==1 and a[0]["b"]==1
    assert a[0]["c"]==1 and a[0]["d"]==1
    b,_=table_strata(source,origin="native")
    assert b==[]


def test_taxonomic_guard():
    assert exact_binomial("Asclepias incarnata")
    assert not exact_binomial("Asclepias spp.")
    assert not exact_binomial("complex speciosa/syriaca")
    assert not exact_binomial("Cynanchum")
    assert canonical_host("curassavica") == "Asclepias curassavica"
    assert canonical_host("syriaca") == "Asclepias syriaca"
    assert canonical_host("complex speciosa/syriaca") is None
    assert canonical_host("Gomphocarpus sp.") is None
    assert canonical_host("Oxypetalum coeruleum") == "Oxypetalum coeruleum"

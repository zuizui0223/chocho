"""Recompute the published 38-port / 17-island 2001–2010 botanical index.

No presumed historical alien status, observed butterfly establishment or
ship-mediated plant transport is allowed in this ecological source audit.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031"
SRC=D/"audit_island_host_index_uncertainty.py"
spec=importlib.util.spec_from_file_location("historical_38port_index",SRC)
assert spec and spec.loader
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def near(a,b):
    assert math.isclose(a,b,rel_tol=0,abs_tol=2e-11),(a,b)


def test_island_equal_historical_comparison_source_and_bootstrap_replays():
    expected=json.loads((D/"historical_38port_island_uncertainty.json").read_text(encoding="utf8"))
    actual=module.analyze()
    assert actual["schema"]==expected["schema"]
    assert actual["status"]==expected["status"]
    assert actual["sources"]==expected["sources"]
    assert actual["source_ports"]==38
    assert actual["source_islands"]==17
    assert actual["geometric_biological_query_sites"]==114
    assert actual["all_queries_complete"]
    assert actual["bootstrap"]["draws"]==5000
    for cutoff in ("1","5","10","20","50"):
        a=actual["species_key_cutoffs"][cutoff]
        b=expected["species_key_cutoffs"][cutoff]
        assert a["eligible_ports"]==b["eligible_ports"]
        assert a["eligible_islands"]==b["eligible_islands"]
        near(a["island_equal_mean_port_minus_offport_top10_index_fraction"],
             b["island_equal_mean_port_minus_offport_top10_index_fraction"])
        assert a["positive_islands"]==b["positive_islands"]
        assert a["negative_islands"]==b["negative_islands"]
        assert a["zero_islands"]==b["zero_islands"]
        if a["eligible_islands"]<3:
            assert a["island_bootstrap95"] is None
            assert b["island_bootstrap95"] is None
        else:
            for x,y in zip(a["island_bootstrap95"],b["island_bootstrap95"]):
                near(x,y)
            assert a["island_bootstrap95"][0] < 0 < a["island_bootstrap95"][1]
        if a["leave_one_island_out_range"] is None:
            assert b["leave_one_island_out_range"] is None
        else:
            for x,y in zip(a["leave_one_island_out_range"],
                           b["leave_one_island_out_range"]):
                near(x,y)
    assert len(expected["limits"])>=9
    assert "NOT necessarily historically introduced" in " ".join(expected["limits"])


def test_small_n_high_bits_bootstrap_is_not_modulo_cycle_degenerate():
    vals=[-0.0241,0.00153,0.05995,0.06492]
    ci=module.island_bootstrap(vals,module.SEED+10,5000)
    assert ci[0]<0<ci[1],ci
    assert ci[1]-ci[0]>0.01

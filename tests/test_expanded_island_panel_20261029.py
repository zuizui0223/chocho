"""Protect unchanged pilot membership, maximin expansion and source receipts."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

PATH = Path(__file__).resolve().parents[1]/"scripts/acquire_expanded_island_panel_20261029.py"


def module():
    assert PATH.exists(), "Expanded acquisition implementation is missing"
    sys.path.insert(0,str(PATH.parent))
    spec = importlib.util.spec_from_file_location("expanded_island_panel",PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def row(key,lon,rank=None):
    r = dict(island_key=key,representative_lon=lon,representative_lat=0,area_km2=50,
             antarctica_could_be_nearer=0,mainland_coast_gap_lower_km=20,
             mainland_coast_gap_upper_km=21,nearest_wpi_coast_distance_lower_km=10,
             nearest_wpi_coast_distance_upper_km=11)
    if rank:
        r.update(selection_rank=rank,selection_stratum="20_to_100:under100:within50")
    return r


def test_maximin_preserves_pilot_two_then_selects_farthest_nearest_distance():
    m = module()
    candidates = [row("a",0),row("d",180),row("b",90),row("c",-90),row("e",45),row("f",10)]
    pilot = [row("a",0,1),row("d",180,2)]
    chosen = m.expand_selection(candidates,pilot,6)
    assert [r["island_key"] for r in chosen] == ["a","d","b","c","e","f"]
    assert m.expand_selection(list(reversed(candidates)),list(reversed(pilot)),6) == chosen
    assert [r["provenance_origin"] for r in chosen[:3]] == ["ORIGINAL_PILOT_REUSED","ORIGINAL_PILOT_REUSED","NEW_EXPANSION"]


def test_missing_original_island_fails_instead_of_replacing_it():
    m = module()
    with pytest.raises(ValueError):
        m.expand_selection([row("a",0),row("b",90)],[row("a",0,1),row("d",180,2)],8)


def test_original_two_stay_even_when_zeros_and_other_candidates_positive():
    m = module()
    candidates = [row("a",0),row("d",180),row("b",90)]
    for r,n in zip(candidates,[0,0,1000]):
        r["record_count"] = n
    chosen = m.expand_selection(candidates,[row("a",0,1),row("d",180,2)],8)
    assert [r["island_key"] for r in chosen] == ["a","d","b"]


def test_cached_request_and_response_must_match_pinned_hashes(tmp_path):
    m = module()
    req = b'{"geometry":"original"}'
    body = b'{"count":0,"results":[],"facets":[]}'
    (tmp_path/"p_request.json").write_bytes(req)
    (tmp_path/"p_response.json").write_bytes(body)
    meta = dict(request_sha256=m.pilot.sha(req),response_sha256=m.pilot.sha(body),retrieved_utc="2026-10-29T12:00:00Z")
    (tmp_path/"p_receipt.json").write_text(json.dumps(meta))
    result,receipt = m.read_cached(tmp_path,"p",req)
    assert result["count"] == 0 and receipt["retrieved_utc"] == "2026-10-29T12:00:00Z"
    with pytest.raises(ValueError):
        m.read_cached(tmp_path,"p",b"changed geometry")
    (tmp_path/"p_response.json").write_bytes(b'{"count":99}')
    with pytest.raises(ValueError):
        m.read_cached(tmp_path,"p",req)


def test_selection_lock_cannot_be_overwritten_on_resume(tmp_path):
    m = module()
    first = m.lock_selection(tmp_path,{"selected":["a"]})
    assert m.lock_selection(tmp_path,{"selected":["a"]}) == first
    with pytest.raises(ValueError):
        m.lock_selection(tmp_path,{"selected":["b"]})

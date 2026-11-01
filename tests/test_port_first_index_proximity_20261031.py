"""Synthetic unit test only; not a result from the unavailable 14.6M-record file."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PROGRAM = ROOT / "provenance/reviewer_defenses/results/island_observational_robustness_20261031/port_first_index_proximity.py"


def load():
    spec = importlib.util.spec_from_file_location("port_first_index_proximity", PROGRAM)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def occ(species, year, lon, island="island_A"):
    return dict(speciesKey=species, landmass_id=island, year=str(year),
                decimalLatitude="0", decimalLongitude=str(lon))


def test_first_index_near_port_relative_to_same_year_background():
    m = load()
    ports = [dict(port_id="WPI_A", landmass_id="island_A",
                  latitude="0", longitude="0")]
    rows = [occ("B" + str(i), y, 1)
            for y in (2009, 2010, 2011, 2012)
            for i in range(12)]
    rows += [occ("Focal", 2011, 0), occ("Focal", 2012, 2)]
    # An unsupported island cannot be re-labelled as "port-free" or zero.
    rows.append(occ("Ghost", 2011, 1, "unmatched_island"))
    report = m.audit(rows, ports, min_other_records=10)
    focal = [x for x in report["cases"] if x["speciesKey"] == "Focal"]
    assert len(focal) == 1
    x = focal[0]
    assert x["first_index_year"] == 2011
    assert x["prior_opportunity_years"] >= 2
    assert x["first_relative_log1p_nearest_port_km"] < 0
    assert x["difference_first_minus_later5y"] < 0
    assert x["later_eligible_years"] == 1
    assert report["record_audit"]["no_geocoded_port_support"] == 1
    assert report["eligible_islands"] == 1
    assert "NOT_ARRIVAL" in report["status"]


def test_missing_spatial_source_not_a_negative_result():
    m = load()
    with pytest.raises(ValueError, match="empty 2019-WPI"):
        m.audit([], [])
    with pytest.raises(ValueError, match="invalid port coordinates"):
        m.audit([], [dict(port_id="1",landmass_id="a",
                          latitude="200",longitude="0")])

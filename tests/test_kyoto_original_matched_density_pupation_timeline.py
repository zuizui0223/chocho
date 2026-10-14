"""All 30x23 rows below are synthetic parser fixtures, not original butterflies."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
from pathlib import Path

import pytest

PATH=Path(__file__).resolve().parents[1]/"scripts/audit_kyoto_original_matched_density_pupation_timeline.py"
SPEC=importlib.util.spec_from_file_location("source_kyoto_timing", PATH)
assert SPEC and SPEC.loader
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def fake_source():
    groups=[]
    for total in (4,8,12):
        groups.extend([(0,total),(total//4,3*total//4),
                       (total//2,total//2),(3*total//4,total//4),(total,0)])
    rows=[]
    for ix,(s,a) in enumerate(groups):
        for rep in range(2):
            cage=str(2*ix+rep+1)
            for step in range(23):
                day=1+2*step
                decline_factor=1.0 if s <= 2 else .55
                cum_a=min(a,int(a*min(1.,step/12.)*decline_factor))
                cum_s=min(s,int(s*min(1.,step/12.)))
                row={k:"0" for k in m.COLUMNS}
                row.update({
                    "plot.no":cage,"treat.no":cage,
                    "s.density":str(s),"a.density":str(a),
                    "date":"2012-09-01",
                    "day":str(day),"total.a":str(a-cum_a),
                    "total.s":str(s-cum_s),
                    "cumul.ap":str(cum_a),"cumul.sp":str(cum_s),
                    "a1":str(a-cum_a),
                    "s1":str(s-cum_s),
                })
                rows.append(row)
    fh=io.StringIO()
    writer=csv.DictWriter(fh,fieldnames=m.COLUMNS)
    writer.writeheader()
    writer.writerows(rows)
    return fh.getvalue().encode("utf-8")


def audit_fake():
    raw=fake_source()
    return m.analyze_source(raw,checksum=hashlib.md5(raw).hexdigest())


def test_exact_cage_not_690_independent_dose_replicates():
    report=audit_fake()
    assert report["whole_cages"]==30
    assert report["cage_x_visit_original_rows"]==690
    assert report["original_year"]==2012
    assert report["new_butterfly_intervention_or_fitness_observations"] is False
    assert report["GEB_PR38_untouched"] is True
    assert report["planned_days_and_actual_source_days"]["42"]==41


def test_only_exact_native_density_strata_are_compared_without_regression_fabrication():
    x=audit_fake()
    pairs=x["native_density_matched_strata"]
    assert set(pairs)=={"2","3","4","6"}
    assert len(pairs)==4
    assert all(p["cages_per_density_cell"]==2 for p in pairs.values())
    for pair in pairs.values():
        assert pair["lower_S_initial_n"]<pair["higher_S_initial_n"]
        assert len(pair["source_cage_ids_lower_S"])==2
        assert len(pair["source_cage_ids_higher_S"])==2
    assert x["equal_native_stratum_weight_summary_by_day"]["original_final"][
        "number_of_source_matched_A_density_strata"
    ]==4
    assert "p_values" not in x


def test_original_source_identity_fails_closed():
    with pytest.raises(ValueError,match="MD5 mismatch"):
        m.analyze_source(fake_source())
    with pytest.raises(ValueError,match="MD5 mismatch"):
        m.analyze_source(fake_source()+b"x")


def test_monotonic_cumulative_pupation_is_enforced():
    raw=fake_source()
    stream=list(csv.DictReader(io.StringIO(raw.decode())))
    first=stream[0]
    for row in stream:
        if row["plot.no"]==first["plot.no"] and int(row["day"])==45:
            row["cumul.ap"]="3"
        if row["plot.no"]==first["plot.no"] and int(row["day"])==43:
            row["cumul.ap"]="4"
    fh=io.StringIO()
    writer=csv.DictWriter(fh,fieldnames=m.COLUMNS)
    writer.writeheader()
    writer.writerows(stream)
    bad=fh.getvalue().encode()
    with pytest.raises(ValueError,match="cumulative pupation dropped"):
        m.analyze_source(bad,checksum=hashlib.md5(bad).hexdigest())


def test_source_has_no_current_host_biomass_or_adult_flight_claim():
    r=audit_fake()
    warning=" ".join(r["warnings"]).lower()
    assert "not adult eclosion" in warning
    assert "no repeated host edible" in warning
    assert "greenhouse ids" in warning
    assert "day-specific p-values" in warning
    assert len(r["source_cage_records"])==30


def test_original_nonnumeric_day_or_count_cannot_be_interpolated():
    raw=fake_source()
    rows=list(csv.DictReader(io.StringIO(raw.decode())))
    rows[0]["day"]=""
    fh=io.StringIO()
    writer=csv.DictWriter(fh,fieldnames=m.COLUMNS)
    writer.writeheader()
    writer.writerows(rows)
    bad=fh.getvalue().encode()
    with pytest.raises(ValueError,match="missing or invalid source integer"):
        m.analyze_source(bad,checksum=hashlib.md5(bad).hexdigest())

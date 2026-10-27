"""Read-only regression protection for source-verifiable monarch site-information gate."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "provenance/archive/monarch/results/monarch_cultivation_winter_source_MH_20261028.json"


def test_monarch_source_gate_and_rejected_inference_remain_frozen():
    d = json.loads(RESULT.read_text(encoding="utf-8"))
    assert d["schema"] == "chocho_monarch_native_host_cultivation_winter_MH_20261028_v01"
    assert d["source"]["original_zip_sha256"] == "4e4def63669b34738bcbfb8ba32ccf6c61d02cc6992b8f1ebfa8bc350a3987db"
    assert d["source"]["author_thinned_rows"] == 12942
    assert d["source"]["field_eligible_rows_after_rearing_origin_cultivation_geography"] == 10237
    m = d["primary_native_host_cultivation_strict"]
    assert (m["candidate_informative_strata"], m["native_host_taxa_in_informative_strata"], m["state_units"]) == (17, 7, 4)
    assert m["source_information_gate"] == "FAIL_PREDECLARED_SOURCE_SUPPORT"
    assert m["largest_state_record_fraction"] > .96
    assert m["bootstrap"]["valid_draws"] < 950 and m["bootstrap"]["interval"] is None
    assert m["n_informative_source_photo_records"] == 200
    a = d["sensitivity_likely_unlikely_included"]
    assert (a["candidate_informative_strata"], a["native_host_taxa_in_informative_strata"], a["state_units"]) == (41, 10, 10)
    assert a["largest_state_record_fraction"] > .7
    assert a["bootstrap"]["interval"][0] < 1 < a["bootstrap"]["interval"][1]
    h = d["sensitivity_within_100km_host_year"]
    assert h["candidate_informative_strata"] == 7 and h["n_informative_source_photo_records"] == 24
    assert all(v["source_information_gate"] == "FAIL_PREDECLARED_SOURCE_SUPPORT"
               for k,v in d.items() if k.startswith(("primary_", "sensitivity_", "secondary_")))
    assert d["decision"] == "INSUFFICIENT_MATCHED_NATIVE_HOST_CULTIVATION_RECORDING_EVIDENCE"

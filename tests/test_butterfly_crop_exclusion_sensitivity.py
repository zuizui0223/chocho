import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_butterfly_crop_exclusion_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("crop_sensitivity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_crop_classification_exact_and_conservative_genus():
    host_meta = {
        "1": {
            "accepted_name": "Medicago sativa",
            "input_host_name": "Medicago sativa",
            "family": "Fabaceae",
        },
        "2": {
            "accepted_name": "Avena barbata",
            "input_host_name": "Avena barbata",
            "family": "Poaceae",
        },
        "3": {
            "accepted_name": "Poa pratensis",
            "input_host_name": "Poa pratensis",
            "family": "Poaceae",
        },
    }
    crop = {
        "exact_binomials": ["Medicago sativa"],
        "genus_wildcards": ["Avena"],
    }
    strict = mod.classify_hosts(host_meta, crop, conservative_genus=False)
    conservative = mod.classify_hosts(host_meta, crop, conservative_genus=True)
    assert strict == {"1"}
    assert conservative == {"1", "2"}


def test_reconstruction_uses_primary_log1p_expansion_estimand():
    descriptors = {
        "Test butterfly": {
            "species": "Test butterfly",
            "host_family_count": 1.0,
            "resolved_host_species": 1,
            "native_resource_units_descriptor": 2,
        }
    }
    pairs = {"Test butterfly": {"host1"}}
    native = {"host1": {"A", "B"}}
    contemporary = {"host1": {"A", "B", "C", "D"}}
    rows = mod.reconstruct(
        descriptors,
        pairs,
        native,
        contemporary,
        excluded_hosts=set(),
        expected_species=None,
    )
    assert len(rows) == 1
    row = rows[0]
    expected = math.log1p(4) - math.log1p(2)
    assert abs(row["log_resource_expansion"] - expected) < 1e-15


def test_no_non_crop_native_resource_is_none_for_ratio_estimand():
    descriptors = {
        "Test butterfly": {
            "species": "Test butterfly",
            "host_family_count": 1.0,
            "resolved_host_species": 1,
            "native_resource_units_descriptor": 1,
        }
    }
    pairs = {"Test butterfly": {"crop"}}
    native = {"crop": {"A"}}
    contemporary = {"crop": {"A", "B"}}
    rows = mod.reconstruct(
        descriptors,
        pairs,
        native,
        contemporary,
        excluded_hosts={"crop"},
        expected_species=None,
    )
    assert rows[0]["native_resource_units"] == 0
    assert rows[0]["log_resource_expansion"] is None


def test_crop_submission_receipt_omits_structural_portfolio_panel():
    receipt = (
        ROOT / "manuscript" / "butterfly_specialization_supplement_v0.2.md"
    ).read_text(encoding="utf-8")
    renderer = (ROOT / "scripts" / "render_crop_sensitivity_si.py").read_text(
        encoding="utf-8"
    )
    assert "Host-contribution architecture" not in receipt
    assert "**B." not in receipt
    assert "broader diets retained" not in receipt.lower()
    assert "median_effective_contributor_number" not in renderer
    assert "architecture_row" not in renderer
    assert "47.0%" in receipt
    assert "46.0%" in receipt
    assert "Supplementary Table S7. Crop-host exclusion sensitivity" in receipt
    assert "| 0.052 | 38 |" in receipt
    assert "| 0.057 | 37 |" in receipt

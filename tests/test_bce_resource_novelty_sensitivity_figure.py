"""Guard resource-identity and host-inventory sample in supplemental figure."""
from __future__ import annotations

import importlib.util
from copy import deepcopy
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"scripts/render_bce_resource_novelty_sensitivity_figure.py"
spec=importlib.util.spec_from_file_location("bce_source_fig",SOURCE)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def original() -> dict:
    return mod.load_validate(
        ROOT/"provenance/reviewer_defenses/results/butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json"
    )


def test_exact_incidence_partitions_and_source_scope() -> None:
    d=original()
    assert d["sample"] == {
        "focal_butterflies":239,
        "BCE_augmented_butterflies":81,
        "candidate_links":1027,
    }
    assert [x["H"] for x in d["scenarios"]] == [75346,103349]
    assert [x["A"] for x in d["scenarios"]] == [46387,67740]
    assert [x["R"] for x in d["scenarios"]] == [14406,20642]
    assert [x["U"] for x in d["scenarios"]] == [14553,14967]
    assert d["newly_added_regions"] - d["reclassified_as_native"] == d["net_new_regions"] == 414
    assert "not true butterfly utilization" in d["denominator_warning"]


def test_missing_or_altered_source_receipt_is_rejected(tmp_path:Path) -> None:
    import json
    path=ROOT/"provenance/reviewer_defenses/results/butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json"
    o=json.loads(path.read_text())
    o["BCE_Clarke_partial_augmented"]["native_absorbed_incidences"] += 1
    corrupted=tmp_path/"corrupted.json"
    corrupted.write_text(json.dumps(o))
    with pytest.raises(ValueError,match="accounting"):
        mod.load_validate(corrupted)

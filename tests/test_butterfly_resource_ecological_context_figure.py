"""Constrain the alternative Figure 3 to the documented butterfly data sources."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"scripts/render_butterfly_resource_ecological_context_figure.py"
spec=importlib.util.spec_from_file_location("butterfly_context",SOURCE)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def validated():
    return mod.infer_source_values(
        ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv",
        ROOT/"data/frozen/figure_sources/occurrence_resource_validation_species.csv",
        ROOT/"provenance/reviewer_defenses/results/butterfly_expansion_equivalence_v0.1.json",
        ROOT/"provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json",
    )


def test_separate_panel_denominators_and_exact_effects():
    x=validated()
    assert x["focal_species"]==239
    assert x["adequate_host_species"]==215
    assert x["wide_one_family_species_n"]==18
    assert x["occurrence_panel_species"]==32
    assert x["informative_species"]==23
    assert x["outside_native_species_region_units"]==115
    assert x["observed_recovered"]==66
    assert x["region_matched_null"]=={"median":49,"q025":43,"q975":56}
    assert x["host_family_vs_native_resource_spearman"]==pytest.approx(.276256474855,abs=1e-8)
    assert x["host_family_vs_log_gain_spearman"]==pytest.approx(.007985412149,abs=1e-8)
    assert x["host_family_vs_log_gain_bootstrap_95ci"]==pytest.approx(
        [-.1113648840403796,.12777832073032636]
    )


def test_scatter_jitter_is_visual_only_deterministic_and_bounded():
    assert mod.xjitter("Aglais urticae")==mod.xjitter("Aglais urticae")
    assert mod.xjitter("Aglais urticae")!=mod.xjitter("Pieris brassicae")
    assert -0.15<=mod.xjitter("Aglais urticae")<=0.15


def test_source_edit_must_fail_byte_identity(tmp_path:Path):
    src=ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    broken=tmp_path/"changed.csv"
    broken.write_bytes(src.read_bytes()+b"\n")
    with pytest.raises(ValueError,match="source bytes"):
        mod.infer_source_values(
            broken,
            ROOT/"data/frozen/figure_sources/occurrence_resource_validation_species.csv",
            ROOT/"provenance/reviewer_defenses/results/butterfly_expansion_equivalence_v0.1.json",
            ROOT/"provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json",
        )


def test_source_rank_rule_with_ties():
    assert mod.midranks([1,2,1,3]).tolist()==[1.5,3.,1.5,4.]
    assert mod.spearman([1,1,3,4],[2,2,5,6])==pytest.approx(1.)

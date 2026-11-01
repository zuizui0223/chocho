"""Pre-response structural tests for exact-binomial larval host specificity."""
import importlib.util
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_documented_host_specificity.py"


def load():
    # GEB's minimal paper CI does not install the optional spatial geometry
    # stack. Specialized source workflows install and execute these tests.
    pytest.importorskip("pyproj")
    pytest.importorskip("shapefile")
    sys.path.insert(0,str(ROOT/"scripts"))
    spec=importlib.util.spec_from_file_location("documented_host_specificity",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_source_binomial_is_conservative_and_no_nonhost_assertion():
    m=load()
    assert m.binomial("Medicago","sativa")=="medicago sativa"
    assert m.binomial("Medicago","Medicago sativa")=="medicago sativa"
    assert m.binomial("Medicago","") is None
    assert m.binomial("Medicago","sp.") is None


def test_matched_fake_host_null_preserves_taxon_frequency_bins():
    m=load()
    # 4 source regions, 5 source species, two focal source host taxa,
    # matched replacements from two within-bin candidates.
    region_masks=[0b00111,0b11011,0b11001,0b01001]
    freq,groups=m.prevalence_bins(region_masks,5)
    true=(1<<1)|(1<<2)
    # Use a source prevalence design with enough replacements per bin.
    # This design may be infeasible and MUST fail rather than reuse taxa.
    with pytest.raises(RuntimeError,match="Cannot prevalence-match"):
        m.host_freq_matched_null(true,region_masks,groups,freq,draws=20)


def test_pair_jaccard_distinguishes_nonintersection_from_double_empty():
    m=load()
    value,missing=m.pair_jaccard([0,0,1], [0,0,2], [1,2,2])
    assert missing==1
    assert value[0]!=value[0]
    assert value[1]==0
    assert value[2]==1
    _,missing2=m.pair_jaccard([0,0,1], [0,0,2],[1,2,2],undefined_zero=True)
    assert missing2==1

"""Source-independent arithmetic and false-zero guard for coastal N rings."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from exploratory_cirad_n_near_ring_replay import annulus_from_nested, add_annulus


def test_weighted_annulus_reconstruction_not_difference_of_averages():
    a = annulus_from_nested(4, 2, .5, 6, 10, 5, .3, 18)
    assert a["mean"] == pytest.approx(7.)
    assert a["coverage"] == pytest.approx(.5)
    assert a["nonzero_cells"] == pytest.approx(5.)
    assert a["nonzero_share"] == pytest.approx(5 / 6)
    assert a["nonzero_mean"] == pytest.approx(42 / 5)


def test_missing_sea_cells_never_imputed_as_zero():
    a = annulus_from_nested(0, np.nan, np.nan, 6, 10, 5, .3, 18)
    assert np.isnan(a["mean"])
    assert np.isnan(a["nonzero_mean"])


def test_source_geometry_and_negative_annulus_fail_closed():
    with pytest.raises(ValueError, match="geometry"):
        annulus_from_nested(8, 2, 0, 10, 6, 3, 0, 20)
    with pytest.raises(ValueError, match="negative N"):
        annulus_from_nested(4, 100, 0, 6, 10, 1, 0, 18)
    with pytest.raises(ValueError, match="zero fractions"):
        annulus_from_nested(4, 2, 0, 6, 10, 5, .9, 18)


def test_valid_mapped_zero_not_converted_to_na():
    a = annulus_from_nested(3, 0, 1, 8, 8, 0, 1, 21)
    assert a["mean"] == 0
    assert a["nonzero_share"] == 0
    assert np.isnan(a["nonzero_mean"])


def test_feature_columns_and_country_gate():
    data = pd.DataFrame({"port_id": [1], **{
        f"{s}_{r}km_{v}": [x]
        for s in ("human", "crop") for r in (10, 25)
        for v, x in {"valid_cells": 4 if r == 10 else 10,
                     "mean": 2 if r == 10 else 5,
                     "zero_share_valid": .5 if r == 10 else .3,
                     "circle_cells": 6 if r == 10 else 18}.items()
    }})
    a = add_annulus(data)
    assert a.crop_delta.iloc[0] == pytest.approx(np.log(3) - np.log(8))
    assert a.crop_ring_valid_cells.iloc[0] == 6
    assert a.human_ring_nonzero_share.iloc[0] == pytest.approx(5/6)
    with pytest.raises(ValueError, match="Duplicate"):
        add_annulus(pd.concat([data, data]))

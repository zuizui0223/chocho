"""Synthetic HOSTS stage fields audit without attributing source codes to larval ecology."""
import csv
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_hosts_lab_rearing_damage_source_20261027 as source


def test_literal_values_preserved_without_false_boolean_imputation():
    assert source.value("")=="<EMPTY>"
    assert source.value(None)=="<EMPTY>"
    assert source.value("  Yes  ")=="Yes"
    assert source.value("NO")=="NO"
    assert source.value("larvae")=="larvae"


def test_wrong_frozen_original_source_header_rejected(tmp_path):
    f=tmp_path/"hosts.csv"
    with f.open("w",newline="") as out:
        w=csv.writer(out)
        w.writerow(["Lab Rearing","Location","Insect Genus",
                    "Insect Species","Hostplant Genus","Hostplant Species"])
        w.writerow(["Yes","Canada","Test","example","Plant","example"])
    with pytest.raises(ValueError,match="header changed"):
        source.inventory(f,{"Test example"})


def test_prespecified_offpanel_focal_species_list_separate():
    assert len(source.FOCAL)==8
    assert {"Chilades pandava","Cacyreus marshalli","Papilio xuthus"} <= source.FOCAL
    assert "Pieris rapae" in source.FOCAL
    assert "Lab Rearing" in source.REQUIRED
    assert "Damage" in source.REQUIRED

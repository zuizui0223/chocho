"""Synthetic tests of plot as independent unit; not a study of butterflies."""
from __future__ import annotations
import importlib.util
import io
import csv
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/audit_monarch_plot_level_survival_uncertainty.py"
spec=importlib.util.spec_from_file_location("plotcheck",SCRIPT)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def synthetic():
    out=io.StringIO()
    cols=["ID","Plant_ID","PlotNum","Temp","Milkweed","OE_treatment","Surv_adult","Notes"]
    w=csv.DictWriter(out,fieldnames=cols)
    w.writeheader()
    for p in range(1,31):
        temp="ambient" if p<=15 else "elevated"
        for milkweed in mod.HOST:
            for treatment in mod.OE:
                plant=f"{p}-{milkweed}-{treatment}"
                for i in range(2):
                    mid="22b" if (p==1 and milkweed=="tropical" and treatment=="control" and i==0) else f"{plant}-{i}"
                    missing=mid=="22b"
                    w.writerow({"ID":mid,"Plant_ID":plant,"PlotNum":p,"Temp":temp,
                                "Milkweed":milkweed,"OE_treatment":treatment,
                                "Surv_adult":"NA" if missing else ("1" if i==0 else "0"),
                                "Notes":"injured" if missing else ""})
    return out.getvalue().encode()

def test_units_and_plot_bootstrap_are_reproducible():
    fake=synthetic()
    original=mod.git_blob_sha
    try:
        mod.git_blob_sha=lambda raw:mod.BLOB
        result=mod.audit(fake)
    finally:
        mod.git_blob_sha=original
    assert result["plots"]==30
    assert result["plants"]==120
    assert result["known_fates"]==239
    assert result["unknown_fates"][0]["ID"]=="22b"
    assert result["can_prove_heterospecific_pathogen_transmission"] is False
    assert result["bootstrap_plot_stratified_95ci"]==mod.bootstrap_ci(mod.plot_differences(mod_parse(fake)))
    assert result["GEB_submission_PR38_unchanged"]

def mod_parse(fake):
    original=mod.git_blob_sha
    try:
        mod.git_blob_sha=lambda raw:mod.BLOB
        return mod.parse(fake)[0]
    finally:
        mod.git_blob_sha=original

def test_real_blob_hash_rejected_if_fake():
    try:
        mod.parse(synthetic())
    except ValueError as err:
        assert "Git blob SHA mismatch" in str(err)
    else:
        raise AssertionError("fake data were called original")

def test_repeated_temperature_per_plot_does_not_pass():
    fake=synthetic().decode()
    fake=fake.replace("30-swamp-control-1,30-swamp-control,30,elevated,swamp",
                      "30-swamp-control-1,30-swamp-control,30,ambient,swamp")
    original=mod.git_blob_sha
    try:
        mod.git_blob_sha=lambda _:mod.BLOB
        try:
            mod.parse(fake.encode())
        except ValueError as err:
            assert "multiple temperatures" in str(err)
        else:
            raise AssertionError("mixed-temperature plot accepted")
    finally:
        mod.git_blob_sha=original

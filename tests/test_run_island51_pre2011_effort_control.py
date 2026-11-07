import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_island51_pre2011_effort_control import git_blob_sha,early_botanical_features,MODELS,COMPARISONS
def test_git_blob_digest():
    import hashlib
    assert git_blob_sha(b"hello\n")==hashlib.sha1(b"blob 6\0hello\n").hexdigest()
def test_model_roles_frozen():
    assert len(MODELS)==6 and len(COMPARISONS)==7
    assert COMPARISONS[3][:3]==("CRE_to_HCRE","CRE","HCRE")
    assert COMPARISONS[5][:3]==("HRE_to_HCRE","HRE","HCRE")
def test_all_original_ports_early_effort_and_full_island_support():
    from exploratory_island51_circuit_network_holdout import source
    from run_island51_host_circuit_joint import INPUT_PATHS,load_original_port_counts
    import json
    repo=Path(__file__).resolve().parents[1]
    q,ids,_=source(repo/INPUT_PATHS["pairs"],repo/INPUT_PATHS["islands"])
    prep=json.loads((repo/INPUT_PATHS["preparation"]).read_text())
    original=load_original_port_counts(repo/INPUT_PATHS["islands"],ids,prep)
    f,coverage,receipt=early_botanical_features(q,ids,original)
    assert len(coverage)==51
    assert receipt["islands_with_early_plant_records"]==44
    assert receipt["islands_without_early_plant_records"]==7
    assert receipt["zero_record_not_absence"]
    assert all(len(x)==1275 and np.isfinite(x).all() for x in f.values())

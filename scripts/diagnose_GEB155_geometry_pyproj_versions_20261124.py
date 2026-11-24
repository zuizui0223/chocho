#!/usr/bin/env python3
"""Read-only original155 distance provenance: controlled pyproj version comparison.

Recompute the ORIGINAL geographic domain from SHA-verified original sources.
Compare only row counts and numerical error counts to archived distances.
No coordinate, distance, species, observer, original panel or model is exported.
"""
import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import recover_GEB_v05_original155_private_inputs_20261122 as recover
import replay_GEB_v05_original155_partial_panel_20261124 as base

VERSIONS=("3.7.2","3.7.1")
SELECTION_SHA="ae0a874b1d743657fdd95a1e653d7ca0939fc41a56b947f8a077be7785c4332d"

class GeometryGuard(ValueError):
    pass

def require(ok,label):
    if not ok:
        raise GeometryGuard(label)

def worker(root,snapshot,version):
    from pyproj import Transformer, proj_version_str
    import pyproj
    from sklearn.neighbors import BallTree
    require(pyproj.__version__==version, 'pyproj_version_was_not_pinned')
    upstream=root/"upstream"
    source=upstream/"work/temporal"
    forecast=source/"forecast"
    selection=snapshot/"forecast_all_225_disposition.csv"
    require(recover.sha256_file(selection)==SELECTION_SHA,"eligible_taxon_source_sha")
    eligible=pd.read_csv(selection)
    selected=eligible[eligible.stage_disposition.str.startswith("BASELINE_ELIGIBLE")]
    require(len(eligible)==225 and len(selected)==155, "eligible_selection_counts")
    cells=pd.read_parquet(forecast/"candidate_cells.parquet")
    env=pd.read_parquet(forecast/"cell_environment_period.parquet")
    initial=env[env.period==0].set_index("cell")
    cells=cells[(cells.center_region_matches==1)&
                (~cells.region_boundary_flag)&
                cells.cell.map(initial.unclassified).le(.05)].copy()
    geo=pd.read_parquet(forecast/"cell_climate_grid.parquet")
    cells=cells.merge(geo[["cell","climate_id","climate_lat","climate_lon"]],
                      on="cell",validate="one_to_one")
    observed=pd.read_parquet(source/"species_year_cells_10km.parquet",
       columns=["speciesKey","year","cell_x_10km","cell_y_10km"])
    observed=observed[observed.year<=2014].drop_duplicates(
       ["speciesKey","cell_x_10km","cell_y_10km"])
    proj=Transformer.from_crs(6933,4326,always_xy=True)
    test_radians=np.radians(cells[["latitude","longitude"]].values)
    items=[]
    for sp in selected.itertuples(index=False):
        prev=observed[observed.speciesKey==sp.taxon_key]
        require(len(prev)>0,"no_origin_species_cells")
        lon,lat=proj.transform((prev.cell_x_10km.to_numpy()+.5)*10000,
                               (prev.cell_y_10km.to_numpy()+.5)*10000)
        tree=BallTree(np.radians(np.column_stack([lat,lon])),metric="haversine")
        distance=tree.query(test_radians,k=1)[0].ravel()*6371.0088
        chosen=cells.loc[distance<=500,["cell"]].copy()
        chosen["speciesKey"]=sp.taxon_key
        chosen["distance_km"]=distance[distance<=500]
        items.append(chosen)
    new=pd.concat(items,ignore_index=True)
    require(len(new)==280988 and new.speciesKey.nunique()==155,"domain_support")
    original=pd.read_parquet(forecast/"species_cell_period_panel.parquet",
              columns=["speciesKey","cell","period","distance_km"])
    reference=original[original.period==0][["speciesKey","cell","distance_km"]]
    require(len(reference)==len(new),"reference_domain_count")
    merged=reference.merge(new,on=["speciesKey","cell"],
                           suffixes=("_original","_current"),
                           validate="one_to_one",sort=False)
    require(len(merged)==280988,"domain_match")
    old=merged.distance_km_original.to_numpy(dtype=np.float64)
    rebuilt=merged.distance_km_current.to_numpy(dtype=np.float64)
    require(np.isfinite(old).all() and np.isfinite(rebuilt).all(),
            "finite_distance_values")
    diff=np.abs(old-rebuilt)
    a32=old.astype(np.float32)
    b32=rebuilt.astype(np.float32)
    old_log=np.log1p(old).astype(np.float32)
    new_log=np.log1p(rebuilt).astype(np.float32)
    original_source_indices_match=bool(reference[["speciesKey","cell"]].
       reset_index(drop=True).equals(new[["speciesKey","cell"]].
       reset_index(drop=True)))
    metrics={
      "status":"PASS_VERIFIED_SOURCE_TO_GEOMETRY_NUMERICS",
      "pyproj_version":pyproj.__version__,
      "proj_version":proj_version_str,
      "sklearn_version":importlib.metadata.version("scikit-learn"),
      "numpy_version":np.__version__,
      "python_version":platform.python_version(),
      "pairs":len(merged),
      "matching_domain_key_set":True,
      "identical_physical_domain_row_order":original_source_indices_match,
      "original_vs_generated_exact_float64_mismatch_pairs":int(np.count_nonzero(diff!=0)),
      "original_vs_generated_float32_mismatch_pairs":int(np.count_nonzero(a32!=b32)),
      "original_vs_generated_model_feature_log1p_float32_mismatch_pairs":
            int(np.count_nonzero(old_log!=new_log)),
      "maximum_absolute_distance_difference_km":float(np.max(diff)),
      "reference_origin_values_copied_to_generated_data":False,
      "source_records_or_generated_numeric_feature_rows_exported":False
    }
    return metrics

def orchestrate(snapshot,output):
    report={
      "schema":"GEB_original155_geo_numeric_pyproj_version_probe_v1",
      "status":"IN_PROGRESS",
      "original_archive_sha256":recover.ARCHIVE_SHA256,
      "original_panel_sha256":recover.PANEL_SHA256,
      "eligible_species_source_sha256":SELECTION_SHA,
      "controlled_versions":VERSIONS,
      "arms":{},
      "raw_records_exported":False
    }
    with tempfile.TemporaryDirectory(prefix="geb155-geometry-probe-") as tmp:
        root=Path(tmp)
        original_source_files=base.SOURCE_FILES
        verification=base.restore_upstream(root)
        require(verification[recover.PANEL]["sha256"]==recover.PANEL_SHA256,
                "reference_original_panel_sha")
        report["source_inputs_sha256"]={k:v["sha256"] for k,v in verification.items()}
        def launch(python,version):
            cmd=[str(python),str(Path(__file__).resolve()),"worker",
                 "--root",str(root),"--snapshot",str(snapshot),
                 "--version",version]
            process=subprocess.run(cmd,capture_output=True,text=True)
            require(process.returncode==0,"controlled_version_worker_failure_"+version.replace(".","_"))
            return json.loads(process.stdout.strip())
        report["arms"]["pyproj_3_7_2"]=launch(sys.executable,"3.7.2")
        env=root/"pyproj371-venv"
        subprocess.run([sys.executable,"-m","venv","--system-site-packages",str(env)],check=True,
                       capture_output=True)
        python=env/"bin/python"
        install=subprocess.run([str(python),"-m","pip","install",
                       "--disable-pip-version-check","--only-binary=:all:","--no-deps",
                       "pyproj==3.7.1"],capture_output=True,text=True)
        require(install.returncode==0,"alternate_pyproj_wheel_install_failed")
        report["arms"]["pyproj_3_7_1"]=launch(python,"3.7.1")
        report["status"]="COMPLETE"
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("stage",choices=("run","worker"))
    p.add_argument("--snapshot",required=True,type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--root",type=Path)
    p.add_argument("--version",choices=VERSIONS)
    a=p.parse_args()
    if a.stage=="worker":
        require(a.root is not None and a.version in VERSIONS,"worker_args")
        print(json.dumps(worker(a.root.resolve(),a.snapshot.resolve(),a.version),
                         sort_keys=True,allow_nan=False))
        return 0
    require(a.output is not None,"missing_output")
    result={"status":"INCOMPLETE","raw_records_exported":False}
    try:result=orchestrate(a.snapshot.resolve(),a.output.resolve())
    except (GeometryGuard,recover.RecoveryError):
        result["status"]="BOUNDARY_CHECK_FAILED"
    except Exception as e:
        result["status"]="INCOMPLETE_ERROR_"+type(e).__name__
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("ORIGINAL155_GEOMETRY_VERSION_PROBE_"+result["status"],flush=True)
    if result.get("arms"):
        print("GEOMETRY_FLOAT32_DIFF_COUNTS",json.dumps({
            k:{"raw":v["original_vs_generated_float32_mismatch_pairs"],
               "log_model":v["original_vs_generated_model_feature_log1p_float32_mismatch_pairs"]}
            for k,v in result["arms"].items()},sort_keys=True),flush=True)
    return 0 if result["status"]=="COMPLETE" else 1

if __name__=="__main__":
    try: sys.exit(main())
    except Exception as e:
        print("GEOMETRY_PROBE_PRIVATE_WORKER_FAIL_"+type(e).__name__,file=sys.stderr)
        sys.exit(1)

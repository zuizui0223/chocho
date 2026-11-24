#!/usr/bin/env python3
"""Read the original155 private source MANIFEST only; export bounded hash aggregates."""
import argparse, hashlib, json, re, tempfile
from pathlib import Path, PurePosixPath
import recover_GEB_v05_original155_private_inputs_20261122 as source

SHA = re.compile(r"^[0-9a-f]{64}$")
ALLOW = {
    "species_cell_period_panel.parquet": "panel",
    "panel_effort20.parquet": "sensitivity_panel_effort20",
    "panel_radius250.parquet": "sensitivity_panel_radius250",
    "panel_exclude_top_provider.parquet": "sensitivity_panel_exclude_top_provider",
    "panel_precision100.parquet": "sensitivity_panel_precision100",
    "predictions_temporal_effort20.parquet": "sensitivity_prediction_effort20",
    "predictions_temporal_radius250.parquet": "sensitivity_prediction_radius250",
    "predictions_temporal_exclude_top_provider.parquet": "sensitivity_prediction_exclude_top_provider",
    "predictions_temporal_precision100.parquet": "sensitivity_prediction_precision100",

    "predictions_species.parquet": "species_predictions",
    "predictions_temporal.parquet": "temporal_predictions",
    "predictions_spatial.parquet": "spatial_predictions",
    "forecast_species_results.json": "species_scores",
    "forecast_model_lock.json": "model_lock",
}
EXPECTED_SHA = {
    "panel": source.PANEL_SHA256,
    "model_lock": source.LOCK_SHA256,
    "species_scores": "3eaaeeb8b77088314bf42051c5a51c07895c17615bf5f27fc3a33424f25e125a",
}
class ScanError(ValueError): pass

def category(raw):
    if not isinstance(raw, str): return None
    name = PurePosixPath(raw.replace("\\", "/")).name.lower()
    if name in ALLOW: return ALLOW[name]
    if re.fullmatch(r"fit_forecast[a-z0-9_.-]*\.py", name): return "fit_forecast_script_candidates"
    if re.fullmatch(r"model_m[0-5]_[a-z0-9_]+\.joblib", name): return "saved_model_candidates"
    if name in ("requirements.txt", "requirements.lock", "environment.yml",
                "environment.yaml", "pyproject.toml", "uv.lock",
                "poetry.lock", "pip-freeze.txt", "pip_freeze.txt"): return "runtime_manifest_candidates"
    return None

def summarize(records):
    if not isinstance(records, list): raise ScanError("bad_manifest")
    groups = {}
    invalid = 0
    for item in records:
        if not isinstance(item, dict): raise ScanError("invalid_record")
        c = category(item.get("path"))
        if c is None: continue
        value,size=item.get("sha256"),item.get("bytes")
        if not isinstance(value,str) or not SHA.fullmatch(value) or type(size) is not int or size < 0:
            invalid += 1
            continue
        key=(value,size)
        groups.setdefault(c,{}).setdefault(key,0)
        groups[c][key]+=1
    if invalid: raise ScanError("invalid_matching_record")
    summary={}
    for kind,counts in sorted(groups.items()):
        summary[kind] = {
            "entries": sum(counts.values()),
            "distinct_content_hashes": len(counts),
            "hash_groups": [{"sha256": sha,"bytes": sz,"count": n,
               "is_expected_digest": (sha == EXPECTED_SHA[kind]) if kind in EXPECTED_SHA else None}
                for (sha,sz),n in sorted(counts.items())],
        }
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    release=source.gh_json(f"repos/{source.REPO}/releases/{source.RELEASE_ID}")
    if release.get("draft") is not True: raise ScanError("draft_release_not_verified")
    assets={x["id"]:x for x in release.get("assets",[])}
    receipt_id,nbytes,expected=source.RECEIPT
    asset=assets.get(receipt_id,{})
    if asset.get("size")!=nbytes or asset.get("digest")!="sha256:"+expected:
        raise ScanError("source_manifest_asset_binding_failed")
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/"original_manifest.json"
        with path.open("wb") as stream:
            source.stream_asset(receipt_id,stream,nbytes,expected,hashlib.sha256(),Path(temp))
        data=json.loads(path.read_text())
        records=data.get("files")
        if len(records)!=3032: raise ScanError("source_manifest_count_changed")
        groups=summarize(records)
        exact_suffix_presence = {}
        for sensitivity in ('effort20','radius250','exclude_top_provider','precision100'):
            for role, relative in (
                ('panel',f'work/temporal/forecast/panel_{sensitivity}.parquet'),
                ('predictions',f'work/temporal/forecast/predictions_temporal_{sensitivity}.parquet')):
                try:
                    source.find_manifest_record(records,relative)
                    match = True
                except source.RecoveryError:
                    match = False
                exact_suffix_presence[sensitivity+'_'+role] = match

        # Original build_forecast_panel.py SHA-bound source requirements.
        # Report only constant known project-relative path labels, not ZIP paths.
        builder_required = [
            'work/temporal/forecast/candidate_cells.parquet',
            'work/temporal/forecast/cell_environment_period.parquet',
            'work/temporal/forecast/cell_climate_grid.parquet',
            'outputs/forecast_candidate_disposition.csv',
            'work/temporal/forecast/species_recording_days.parquet',
            'work/temporal/forecast/target_group_recording_days.parquet',
            'work/temporal/species_year_cells_10km.parquet',
            'work/temporal/host_region_composition_metrics.parquet',
            'work/remote_sensing/host_region_counts.csv',
            'work/temporal/host_region_family_composition.parquet',
            'outputs/temporal_host_family_support.csv',
            'outputs/forecast_model_lock.json',
            'work/temporal/forecast/species_cell_period_panel.parquet',
        ]
        builder_inputs = {}
        for relative in builder_required:
            try:
                entry = source.find_manifest_record(records, relative)
                builder_inputs[relative] = {
                    'present': True, 'sha256': entry['sha256'],
                    'bytes': entry['bytes'],
                }
            except source.RecoveryError:
                builder_inputs[relative] = {'present': False}
        host_counts_basename_candidates = sum(
            isinstance(entry.get('path'), str)
            and re.fullmatch(r'host[_-]?region.*counts\.csv',
                PurePosixPath(entry['path'].replace('\\\\','/')).name.lower()) is not None
            for entry in records if isinstance(entry, dict)
        )
        occurrence_parts_entries = sum(
            isinstance(entry.get('path'), str)
            and '/work/temporal/occurrence_parts/' in
                ('/' + entry['path'].replace('\\\\', '/'))
            for entry in records if isinstance(entry, dict)
        )

        # 24-Nov later historical-runtime search. Aggregate only, never source paths.
        historical_runtime_categories = {
            'environment_lock_candidate': 0,
            'pip_or_conda_version_snapshot_candidate': 0,
            'original_forecast_execution_log_candidate': 0,
            'distance_intermediate_candidate': 0,
            'geo_library_metadata_candidate': 0,
        }
        for item in records:
            if not isinstance(item, dict) or not isinstance(item.get('path'), str):
                continue
            original_name = item['path'].replace('\\', '/').lower()
            name = PurePosixPath(original_name).name
            if (name in ('requirements.txt', 'environment.yml', 'environment.yaml',
                         'conda-lock.yml', 'conda-lock.yaml', 'pipfile.lock',
                         'poetry.lock', 'uv.lock', '.python-version',
                         'pyproject.toml', 'pixi.lock')
                    or re.fullmatch(r'requirements[-_a-z0-9.]*\.(txt|in)', name)):
                historical_runtime_categories['environment_lock_candidate'] += 1
            if ('pip' in name or 'conda' in name or 'package' in name
                    or 'python_version' in name) and name.endswith(('.json', '.txt', '.log', '.csv')):
                historical_runtime_categories['pip_or_conda_version_snapshot_candidate'] += 1
            if any(t in name for t in ('forecast', 'panel', 'distance')) and name.endswith('.log'):
                historical_runtime_categories['original_forecast_execution_log_candidate'] += 1
            if ('distance' in name or 'nearest' in name) and name.endswith(
                    ('.parquet', '.csv', '.npy', '.npz', '.pkl', '.feather')):
                historical_runtime_categories['distance_intermediate_candidate'] += 1
            if any(t in name for t in ('pyproj', 'proj_version', 'sklearn_version',
                                       'environment_versions', 'runtime_versions')):
                historical_runtime_categories['geo_library_metadata_candidate'] += 1

    result={
        "schema":"geb_original155_manifest_only_code_lineage_v1",
        "release_is_draft":True,
        "manifest_sha256":expected,
        "manifest_bytes":nbytes,
        "source_record_count":3032,
        "original_zip_not_downloaded":True,
        "matching_source_groups":groups,
        "expected_relative_path_suffixes":exact_suffix_presence,
        "original_builder_inputs":builder_inputs,
        "historical_runtime_and_distance_candidates":historical_runtime_categories,
        "original_builder_inputs_present":sum(x["present"] for x in builder_inputs.values()),
        "original_builder_inputs_required":len(builder_inputs),
        "host_region_counts_csv_basename_candidates":host_counts_basename_candidates,
        "occurrence_parts_parquet_source_entries":occurrence_parts_entries,
        "no_private_paths_or_source_rows_exported":True,
        "limits":"Verified manifest paths and registered digests only; existence does not establish historical execution time or code version.",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MANIFEST_ONLY_CLASS_COUNTS",json.dumps({k:v["entries"] for k,v in groups.items()},sort_keys=True))
    print("HISTORICAL_RUNTIME_SOURCE_CANDIDATES",json.dumps(historical_runtime_categories,sort_keys=True))
    print("BUILDER_EXACT_INPUTS",len(builder_inputs),sum(x["present"] for x in builder_inputs.values()),"host_counts_candidates",host_counts_basename_candidates,"occurrence_parts",occurrence_parts_entries)
    print("MANIFEST_ONLY_SOURCE_SCRIPT_HASHES",groups.get("fit_forecast_script_candidates",{}).get("distinct_content_hashes",0))

if __name__=="__main__":
    try: main()
    except (ScanError,source.RecoveryError) as exc: raise SystemExit("MANIFEST_AUDIT_FAILED_"+str(exc))
    except Exception as exc: raise SystemExit("MANIFEST_AUDIT_ERROR_"+type(exc).__name__)

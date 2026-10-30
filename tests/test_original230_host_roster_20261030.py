"""Source-bound host recovery: synthetic identities, counts and release boundaries."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import stat
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/derive_original230_host_roster_20261030.py"
HOST_MEMBER = "inputs/original/insect_host_accepted.csv"
METRIC_MEMBER = "inputs/anthropogenic_species_metrics.csv"
OUTPUTS = {
    "original230_host_edges.csv", "original239_host_count_alignment.csv",
    "plant_taxon_query_candidates.csv", "baseline93_host_coverage.csv",
    "original230_host_roster_manifest.json",
}
PINS = {
    "HOSTS": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
    "WCVP": "65bed76bae9d644ccb6ad200c05f9f5071d89e05",
    "TDWG": "52da7828aba9d461dd133c27b3bd7a4407161f54",
}


def module():
    assert SCRIPT.is_file(), "Original230 source-roster recovery is not implemented"
    spec = importlib.util.spec_from_file_location("original230_roster", SCRIPT)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def sha(data):
    return hashlib.sha256(data).hexdigest()


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def synthetic_inputs():
    names = ["Alpha alba", "Beta beta", "Gamma gamma", "Delta delta"]
    bridge = []
    for name, key in zip(names, ["A", "B", "C", "D"]):
        bridge.append({
            "original_name": name, "current_usage_key": key,
            "current_usage_canonical_name": name, "current_usage_rank": "SPECIES",
            "current_usage_status": "ACCEPTED", "current_accepted_key": key,
            "current_accepted_canonical_name": name, "current_accepted_rank": "SPECIES",
            "synonym": "False", "match_type": "EXACT", "confidence": "99",
            "Papilionoidea_verified": "True", "provisional_eligible": "True",
            "rejection_reasons": "", "current_accepted_key_collision": "False",
            "original_numeric_collision": "False", "conservative_join_eligible": "True",
            "response_sha256": sha(name.encode()), "retrieved_utc": "2026-10-29T00:00:00Z",
        })
    bridge[2].update(current_usage_status="SYNONYM", synonym="True",
                     current_accepted_canonical_name="Gamma novata")
    bridge[3].update(current_usage_rank="GENUS", current_accepted_rank="GENUS",
                     provisional_eligible="False", conservative_join_eligible="False",
                     rejection_reasons="not_species_rank")
    hosts = []
    for butterfly, pid, name, family in [
        ("Alpha alba", "11", "Plantus alba", "Flora"),
        ("Alpha alba", "11", "Plantus alba", "Flora"),  # source alias, same ID
        ("Alpha alba", "12", "Plantus beta", "Flora"),
        ("Beta beta", "12", "Plantus beta", "Flora"),
        ("Gamma gamma", "13", "Planta gamma", "Gammae"),
        ("Delta delta", "99", "Excluded delta", "Deltae"),
        ("Unrelated insect", "900", "Unrelated plant", "Other"),
    ]:
        hosts.append({"input_host_name": "RAW_PROVIDER_ALIAS_DO_NOT_EXPORT",
                      "insect_species": butterfly, "accepted_plant_name_id": pid,
                      "accepted_name": name, "family": family})
    metrics = [{"species": name, "resolved_host_species": str(n)}
               for name, n in zip(names, [2, 1, 1, 1])]
    facets = [
        {"island_key": "island1", "window": "2010_2017", "facet_key": "A", "record_count": "3"},
        {"island_key": "island1", "window": "2018_2025", "facet_key": "B", "record_count": "9"},
        {"island_key": "island2", "window": "2010_2017", "facet_key": "D", "record_count": "2"},
    ]
    expected = {"original_taxa": 4, "eligible_taxa": 3, "excluded_taxa": 1,
                "original_edges": 5, "eligible_edges": 4, "baseline_taxa": 1,
                "baseline_edges": 2}
    return hosts, metrics, bridge, facets, expected


def derive(m, inputs=None):
    hosts, metrics, bridge, facets, expected = inputs or synthetic_inputs()
    return m.derive_roster(hosts, metrics, bridge, facets, expected_counts=expected)


def test_all_conservative_taxa_recovered_without_late_or_ever_recorded_selection():
    m = module()
    result = derive(m)
    edges = result["tables"]["original230_host_edges.csv"]
    assert {(r["original_name"], r["wcvp_accepted_id"]) for r in edges} == {
        ("Alpha alba", "11"), ("Alpha alba", "12"),
        ("Beta beta", "12"), ("Gamma gamma", "13"),
    }
    assert result["summary"]["original_catalog_edges"] == 5
    assert result["summary"]["exported_catalog_edges"] == 4
    assert len(edges) == 4  # Duplicate accepted ID cannot inflate the count.
    gamma = next(r for r in edges if r["original_name"] == "Gamma gamma")
    assert gamma["current_accepted_name"] == "Gamma novata"
    assert gamma["current_usage_status"] == "SYNONYM"
    assert gamma["original_name"] == "Gamma gamma"


def test_baseline_coverage_uses_only_early_records_and_does_not_filter_roster():
    m = module()
    result = derive(m)
    coverage = result["tables"]["baseline93_host_coverage.csv"]
    assert [r["original_name"] for r in coverage] == ["Alpha alba"]
    assert coverage[0]["recovered_host_count"] == 2
    assert result["summary"]["baseline_catalog_edges"] == 2
    assert result["summary"]["eligible_taxa"] == 3


@pytest.mark.parametrize("mutation", ["count", "missing_host", "missing_bridge", "duplicate_metric"])
def test_incomplete_original239_reconciliation_fails_before_export(mutation):
    m = module()
    inputs = synthetic_inputs()
    hosts, metrics, bridge, _, _ = inputs
    if mutation == "count":
        metrics[3]["resolved_host_species"] = "2"  # Even excluded originals must reconcile.
    elif mutation == "missing_host":
        hosts[:] = [r for r in hosts if r["insect_species"] != "Gamma gamma"]
    elif mutation == "missing_bridge":
        bridge.pop()
    else:
        metrics.append(dict(metrics[0]))
    with pytest.raises(ValueError):
        derive(m, inputs)


@pytest.mark.parametrize("bad_id", ["", "11;12", "NA", "-1"])
def test_missing_or_nonatomic_focal_host_ids_are_not_zero_or_invented_links(bad_id):
    m = module()
    inputs = synthetic_inputs()
    inputs[0][4]["accepted_plant_name_id"] = bad_id
    with pytest.raises(ValueError, match="accepted plant ID"):
        derive(m, inputs)


def test_source_id_format_diagnostic_contains_counts_without_provider_values():
    m = module()
    values = ["", "23", "23.0", "23.5", "2.3e+06", "23-1", "23;24", "NA", "PRIVATE_VALUE"]
    rows = [{"insect_species": "PRIVATE_FOCAL", "accepted_plant_name_id": value}
            for value in values]
    rows.append({"insect_species": "PRIVATE_OTHER", "accepted_plant_name_id": "42"})
    audit = m.source_id_format_audit(rows, {"PRIVATE_FOCAL"})
    assert audit["original_taxon_rows"] == {
        "empty": 1, "digits": 1, "decimal_integral": 1, "decimal_fractional": 1,
        "scientific_notation": 1, "hyphenated_numeric": 1, "multiple_id_separator": 1,
        "missing_sentinel": 1, "other": 1,
    }
    assert audit["all_source_rows"] == dict(audit["original_taxon_rows"], digits=2)
    encoded = json.dumps(audit)
    assert "PRIVATE" not in encoded
    assert "23" not in encoded


def test_failing_source_id_preserves_strict_gate_and_attaches_aggregate_audit():
    m = module()
    inputs = synthetic_inputs()
    inputs[0][4]["accepted_plant_name_id"] = "13.0"
    with pytest.raises(ValueError, match="accepted plant ID") as failure:
        derive(m, inputs)
    assert failure.value.format_counts == {
        "original_taxon_rows": {"decimal_integral": 1, "digits": 5},
        "all_source_rows": {"decimal_integral": 1, "digits": 6},
    }
    assert "13.0" not in str(failure.value)


def test_id_diagnostic_is_saved_only_in_incomplete_receipt(tmp_path, monkeypatch):
    m = module()
    workflow, artifact = metadata()
    run_path, artifact_path = tmp_path / "run.json", tmp_path / "artifact.json"
    run_path.write_text(json.dumps(workflow))
    artifact_path.write_text(json.dumps(artifact))
    audit = {"original_taxon_rows": {"scientific_notation": 1},
             "all_source_rows": {"digits": 42, "scientific_notation": 1}}

    def fail_identity(*args, **kwargs):
        raise m.AcceptedPlantIdError(audit)

    monkeypatch.setattr(m, "verify_source_metadata", fail_identity)
    output = tmp_path / "derived"
    status = m.main(["--repo-root", str(ROOT), "--run-metadata", str(run_path),
                     "--artifact-metadata", str(artifact_path), "--output-dir", str(output)])
    assert status == 1
    assert {p.name for p in output.iterdir()} == {"original230_host_roster_manifest.json"}
    receipt = json.loads((output / "original230_host_roster_manifest.json").read_text())
    assert receipt["source_id_format_counts"] == audit
    assert receipt["qualified_host_roster"] is False


def test_conflicting_names_preserved_without_selecting_a_query_name():
    m = module()
    inputs = synthetic_inputs()
    inputs[0][1]["accepted_name"] = "Plantus other"
    result = derive(m, inputs)
    candidate = next(r for r in result["tables"]["plant_taxon_query_candidates.csv"]
                     if r["wcvp_accepted_id"] == "11")
    assert json.loads(candidate["accepted_names_json"]) == ["Plantus alba", "Plantus other"]
    assert "MULTIPLE_ACCEPTED_NAMES" in candidate["ambiguity_flags"]
    assert candidate["query_eligible"] is False
    assert candidate["query_name"] == ""
    assert result["summary"]["exported_catalog_edges"] == 4


def test_same_name_multiple_ids_is_audited_without_adding_inferred_edges():
    m = module()
    inputs = synthetic_inputs()
    for row in inputs[0]:
        if row["accepted_plant_name_id"] == "12":
            row["accepted_name"] = "Plantus alba"
    result = derive(m, inputs)
    candidates = result["tables"]["plant_taxon_query_candidates.csv"]
    for row in candidates:
        if row["wcvp_accepted_id"] in {"11", "12"}:
            assert "NAME_MAPS_TO_MULTIPLE_ACCEPTED_IDS" in row["ambiguity_flags"]
            assert row["query_eligible"] is False
    edges = result["tables"]["original230_host_edges.csv"]
    assert not any(r["original_name"] == "Beta beta" and r["wcvp_accepted_id"] == "11"
                   for r in edges)


@pytest.mark.parametrize("name,flag", [
    ("", "MISSING_ACCEPTED_NAME"),
    ("Planta gamma var. minor", "NOT_EXACT_SPECIES_BINOMIAL"),
])
def test_missing_or_infraspecific_names_remain_unresolved_candidates(name, flag):
    m = module()
    inputs = synthetic_inputs()
    inputs[0][4]["accepted_name"] = name
    candidate = next(r for r in derive(m, inputs)["tables"]["plant_taxon_query_candidates.csv"]
                     if r["wcvp_accepted_id"] == "13")
    assert flag in candidate["ambiguity_flags"]
    assert candidate["query_eligible"] is False


def test_current_concept_collision_cannot_be_marked_conservatively_eligible():
    m = module()
    inputs = synthetic_inputs()
    inputs[2][1]["current_accepted_key"] = "A"
    with pytest.raises(ValueError, match="collision"):
        derive(m, inputs)


def test_nine_exclusions_are_accounted_for_without_exporting_their_host_edges():
    m = module()
    result = derive(m)
    alignment = result["tables"]["original239_host_count_alignment.csv"]
    excluded = next(r for r in alignment if r["original_name"] == "Delta delta")
    assert excluded["conservative_join_eligible"] is False
    assert excluded["expected_host_count"] == excluded["recovered_host_count"] == 1
    assert excluded["exclusion_reasons"] == "not_species_rank"
    encoded = json.dumps(result["tables"], sort_keys=True)
    assert "Excluded delta" not in encoded
    assert "Unrelated plant" not in encoded
    assert "RAW_PROVIDER_ALIAS_DO_NOT_EXPORT" not in encoded


def metadata():
    workflow = {"id": 39006554609, "head_sha": "d53e04fb361701265139ba77b8bca56689da6532",
                "status": "completed", "conclusion": "success",
                "repository": {"id": 1390618992, "full_name": "zuizui0223/chocho"},
                "head_repository": {"id": 1390618992, "full_name": "zuizui0223/chocho"}}
    artifact = {"id": 12090584135,
                "name": "butterfly-author-anonymous-offline-source-replication-v01",
                "size_in_bytes": 3250938, "expired": False,
                "digest": "sha256:a784770f5a215bffb32393c127b66e7320e8ed7ac9b36f921bb9c056f3497a36",
                "created_at": "2026-10-17T05:41:00Z", "expires_at": "2026-11-16T05:40:59Z",
                "workflow_run": {"id": 39006554609, "repository_id": 1390618992,
                                 "head_repository_id": 1390618992,
                                 "head_sha": "d53e04fb361701265139ba77b8bca56689da6532"}}
    return workflow, artifact


def test_verified_metadata_retains_source_identity_without_raw_response_payloads():
    m = module()
    workflow, artifact = metadata()
    workflow["actor"] = {"login": "SOURCE_METADATA_DO_NOT_COPY"}
    receipt = m.verify_source_metadata(workflow, artifact,
                                       now=datetime(2026, 10, 30, tzinfo=timezone.utc))
    assert receipt["artifact_id"] == 12090584135
    assert receipt["run_head_sha"] == "d53e04fb361701265139ba77b8bca56689da6532"
    assert "SOURCE_METADATA_DO_NOT_COPY" not in json.dumps(receipt)


@pytest.mark.parametrize("mutation", ["head", "digest", "run", "expired", "size", "conclusion"])
def test_wrong_source_run_or_artifact_metadata_is_rejected(mutation):
    m = module()
    workflow, artifact = metadata()
    if mutation == "head":
        workflow["head_sha"] = "0" * 40
    elif mutation == "digest":
        artifact["digest"] = "sha256:" + "0" * 64
    elif mutation == "run":
        artifact["workflow_run"]["id"] += 1
    elif mutation == "expired":
        artifact["expired"] = True
    elif mutation == "size":
        artifact["size_in_bytes"] += 1
    else:
        workflow["conclusion"] = "failure"
    with pytest.raises(ValueError):
        m.verify_source_metadata(workflow, artifact,
                                 now=datetime(2026, 10, 30, tzinfo=timezone.utc))


def test_stale_expiry_is_rejected_even_if_api_expired_flag_is_false():
    m = module()
    workflow, artifact = metadata()
    with pytest.raises(ValueError, match="expired"):
        m.verify_source_metadata(workflow, artifact,
                                 now=datetime(2026, 11, 17, tzinfo=timezone.utc))


def test_failed_run_in_precreated_empty_directory_emits_only_incomplete_manifest(tmp_path):
    m = module()
    workflow, artifact = metadata()
    workflow["head_sha"] = "0" * 40
    run_path, artifact_path = tmp_path / "run.json", tmp_path / "artifact.json"
    run_path.write_text(json.dumps(workflow))
    artifact_path.write_text(json.dumps(artifact))
    output = tmp_path / "derived"
    output.mkdir()
    status = m.main(["--repo-root", str(ROOT), "--run-metadata", str(run_path),
                     "--artifact-metadata", str(artifact_path), "--output-dir", str(output)])
    assert status == 1
    assert {p.name for p in output.iterdir()} == {"original230_host_roster_manifest.json"}
    manifest = json.loads((output / "original230_host_roster_manifest.json").read_text())
    assert manifest["status"] == "INCOMPLETE_SOURCE_OR_COHORT_VALIDATION_FAILED"
    assert manifest["qualified_host_roster"] is False


def make_archive(tmp_path, mutation=None):
    hosts, metrics, _, _, _ = synthetic_inputs()
    members = {HOST_MEMBER: csv_bytes(hosts), METRIC_MEMBER: csv_bytes(metrics),
               "code/example.py": b"# Archived code must never be executed.\n"}
    if mutation == "unsafe":
        members["../escape.csv"] = b"bad\n"
    manifest = {"schema": "butterfly_anonymous_offline_resource_package_v0.1",
                "status": "AUTHOR_REVIEW_ANONYMISED_SOURCE_CONDITIONAL_NOT_A_PUBLIC_DATA_DEPOSIT",
                "sha256": {k: sha(v) for k, v in members.items()},
                "upstream_source_pins": PINS,
                "sources_conditioned_on_host_catalogue": True,
                "source_licence_approval_required_before_distribution": True}
    expected_members = {k: sha(members[k]) for k in (HOST_MEMBER, METRIC_MEMBER)}
    if mutation == "member_tamper":
        members["code/example.py"] = b"changed\n"
    if mutation == "missing_manifest_member":
        manifest["sha256"].pop("code/example.py")
    if mutation == "upstream":
        manifest["upstream_source_pins"] = dict(PINS, HOSTS="0" * 40)
    members["SOURCE_MANIFEST.json"] = json.dumps(manifest).encode()
    path = tmp_path / "source.zip"
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, payload in members.items():
            if mutation == "symlink" and name == "code/example.py":
                info = zipfile.ZipInfo(name)
                info.create_system = 3
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
                archive.writestr(info, payload)
            else:
                archive.writestr(name, payload)
    return path, {"expected_sha256": sha(path.read_bytes()),
                  "expected_size": path.stat().st_size, "expected_members": expected_members}


def test_archive_verifies_all_members_and_reads_original_table_without_extraction(tmp_path):
    m = module()
    path, pins = make_archive(tmp_path)
    source = m.read_verified_archive(path, **pins)
    assert len(source["host_rows"]) == 7
    assert len(source["metric_rows"]) == 4
    assert source["provenance"]["manifest_members_verified"] == 3
    assert sorted(p.name for p in tmp_path.iterdir()) == ["source.zip"]


@pytest.mark.parametrize("mutation", ["unsafe", "member_tamper", "missing_manifest_member", "upstream", "symlink"])
def test_archive_path_or_member_drift_prevents_catalogue_use(tmp_path, mutation):
    m = module()
    path, pins = make_archive(tmp_path, mutation)
    with pytest.raises(ValueError):
        m.read_verified_archive(path, **pins)


@pytest.mark.parametrize("field,value", [("expected_sha256", "0" * 64), ("expected_size", 1)])
def test_inner_archive_identity_checked_before_reading(tmp_path, field, value):
    m = module()
    path, pins = make_archive(tmp_path)
    pins[field] = value
    with pytest.raises(ValueError):
        m.read_verified_archive(path, **pins)


def test_export_contains_exact_five_derived_files_and_matching_csv_hashes(tmp_path):
    m = module()
    manifest = m.export_roster(derive(m), tmp_path / "derived", provenance={"synthetic": True})
    assert {p.name for p in (tmp_path / "derived").iterdir()} == OUTPUTS
    assert manifest["status"] == "COMPLETE_ORIGINAL230_SOURCE_CATALOGUE_NOT_LOCAL_HOST_USE"
    for name, digest in manifest["output_sha256"].items():
        assert sha((tmp_path / "derived" / name).read_bytes()) == digest
    assert len(manifest["output_sha256"]) == 4
    assert manifest["no_occurrence_queries_performed"] is True


@pytest.mark.parametrize("mutation", ["extra_table", "raw_column", "preexisting_raw"])
def test_export_fails_closed_on_unallowlisted_table_column_or_preexisting_file(tmp_path, mutation):
    m = module()
    result = derive(m)
    output = tmp_path / "derived"
    if mutation == "extra_table":
        result["tables"]["insect_host_accepted.csv"] = [{"secret": "raw"}]
    elif mutation == "raw_column":
        result["tables"]["original230_host_edges.csv"][0]["decimalLatitude"] = "1.0"
    else:
        output.mkdir()
        (output / "source.zip").write_bytes(b"source archive must never be uploaded")
    with pytest.raises(ValueError):
        m.export_roster(result, output, provenance={})
    assert not (output / "original230_host_edges.csv").exists()

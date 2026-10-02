#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

COPY_PATHS = (
    "pyproject.toml",
    "src",
    "scripts",
    "tests",
    "data/external",
    "data/frozen",
    "docs/exploratory",
    "docs/REPRODUCIBILITY.md",
    "benchmarks/exploratory",
    "provenance/reviewer_defenses",
    "provenance/climate",
    "provenance/crop",
    "provenance/archive/climate/operations",
    "provenance/archive/climate/pilot",
    "manuscript/generated/butterfly_specialization_ecology_blinded_v0.2.md",
    "manuscript/butterfly_specialization_claim_map_v0.2.json",
    "manuscript/butterfly_specialization_supplement_v0.2.md",
)

EXCLUDE_RELATIVE = {
    Path("scripts/build_anonymous_review_bundle.py"),
    Path("scripts/release_preflight.py"),
    Path("scripts/submission_preflight.py"),
    Path("tests/test_release_metadata.py"),
    Path("tests/test_reproducibility_snapshot_integrity.py"),
    Path("tests/test_butterfly_specialization_submission_bundle.py"),
    Path("tests/test_butterfly_specialization_manuscript_bundle.py"),
    Path("tests/test_butterfly_climate_release_postgate_recovery.py"),
    Path("tests/test_butterfly_climate_release_execution_binding.py"),
    Path("tests/test_repository_layout.py"),
    Path("provenance/reviewer_defenses/GEB_V02_DECISION_MEMO.md"),
    Path("provenance/reviewer_defenses/REVIEW_RESPONSE_MAP.md"),
}

IDENTITY_REPLACEMENTS = {
    "zuizui0223": "anonymous-user",
    "ZHANG": "[BLINDED]",
    "Ruiqi": "[BLINDED]",
}

PROVENANCE_KEYS = {
    "workflow",
    "workflow_run_id",
    "run_id",
    "artifact_id",
    "artifact_name",
    "artifact_digest",
    "head_sha",
    "workflow_head_sha",
    "render_workflow_run_id",
    "descriptor_workflow_run_id",
}

TEXT_SUFFIXES = {
    ".py", ".md", ".json", ".toml", ".txt", ".yml", ".yaml", ".cff", ".tsv"
}

RAW_PRESERVE_PREFIXES = (
    "data/external/",
    "data/frozen/",
)

ANON_README = """# Butterfly specialization ecology — anonymous review snapshot

This archive is the anonymized data-and-code review supplement intended to be
uploaded directly with the manuscript during double-anonymous peer review. It is
a review snapshot of the analysis package for a global butterfly-specialization
manuscript.

It contains the scientific code, frozen protocols, de-identified result receipts, the v0.2 supplementary robustness tables,
the exact LepTraits input used by the reconstruction, and byte-exact source inputs
needed to regenerate the three main manuscript figures and the climate supplementary figure. Identifying repository metadata,
Git history, author metadata, cover letters, title pages, release metadata, and
internal hosting/run identifiers are intentionally omitted or redacted for review.

## Quick checks

Python 3.11+ is required; the reference CI uses Python 3.12.

    python -m pip install --upgrade pip
    python -m pip install -e ".[test,figure]"
    pytest -q

The two WCVP/HOSTS sidecar scripts use base R only. Full environment notes are in
docs/REPRODUCIBILITY.md.

## Regenerate the manuscript figures

    python scripts/render_butterfly_specialization_v02_figures.py \
      --anthropogenic-csv data/frozen/figure_sources/anthropogenic_species_metrics.csv \
      --matched-null-json provenance/reviewer_defenses/results/butterfly_resource_expansion_matched_null_v0.2.json \
            --hostbias-null-json provenance/reviewer_defenses/results/butterfly_resource_expansion_hostbias_null_v0.1.json \
      --plant-prominence-json provenance/reviewer_defenses/results/butterfly_host_plant_prominence_expansion_v0.1.json \
      --host-concentration-json provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json \
      --occurrence-csv data/frozen/figure_sources/occurrence_resource_validation_species.csv \
      --occurrence-null-json provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json \
      --occurrence-species-robustness-json provenance/reviewer_defenses/results/butterfly_occurrence_species_robustness_v0.1.json \
      --ceiling-json provenance/reviewer_defenses/results/butterfly_expansion_ceiling_sensitivity_v0.1.json \
      --regional-json provenance/reviewer_defenses/results/butterfly_regional_robustness_v0.1.json \
      --climate-csv data/frozen/figure_sources/climate_distance_sensitivity_species.csv \
      --climate-effect-json provenance/reviewer_defenses/results/butterfly_climate_effect_size_v0.1.json \
      --output-dir results/manuscript-figures

The review snapshot intentionally does not contain author names or a public
repository URL. A public, citable archive is supplied only in the final public
version after double-anonymous review.
"""


def sanitize_text(text: str) -> str:
    for old, new in IDENTITY_REPLACEMENTS.items():
        text = text.replace(old, new)
    text = text.replace("https://github.com/anonymous-user/chocho", "[ANONYMIZED_REPOSITORY]")
    text = text.replace("https://github.com/anonymous-user/TTF", "[ANONYMIZED_PRECURSOR_REPOSITORY]")
    return text


def redact_json(value):
    if isinstance(value, dict):
        out = {}
        for key, item in value.items():
            if key in PROVENANCE_KEYS:
                continue
            out[key] = redact_json(item)
        return out
    if isinstance(value, list):
        return [redact_json(item) for item in value]
    if isinstance(value, str):
        return sanitize_text(value)
    return value


def should_preserve_raw(rel: Path) -> bool:
    s = rel.as_posix()
    return any(s.startswith(prefix) for prefix in RAW_PRESERVE_PREFIXES)


def copy_one(source: Path, destination: Path, rel: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if should_preserve_raw(rel):
        shutil.copyfile(source, destination)
        return

    suffix = source.suffix.lower()
    if suffix == ".json":
        try:
            payload = json.loads(source.read_text(encoding="utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            shutil.copyfile(source, destination)
            return
        destination.write_text(
            json.dumps(redact_json(payload), indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        return

    if suffix in TEXT_SUFFIXES or source.name == "pyproject.toml":
        try:
            text = source.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            shutil.copyfile(source, destination)
            return
        destination.write_text(sanitize_text(text), encoding="utf-8")
        return

    shutil.copyfile(source, destination)


def stage_snapshot(stage: Path) -> None:
    for item in COPY_PATHS:
        source = ROOT / item
        if not source.exists():
            raise FileNotFoundError(source)
        if source.is_file():
            rel = source.relative_to(ROOT)
            if rel in EXCLUDE_RELATIVE:
                continue
            copy_one(source, stage / rel, rel)
            continue

        for source_file in sorted(p for p in source.rglob("*") if p.is_file()):
            rel = source_file.relative_to(ROOT)
            if rel in EXCLUDE_RELATIVE:
                continue
            if "__pycache__" in rel.parts or rel.suffix in {".pyc", ".pyo"}:
                continue
            copy_one(source_file, stage / rel, rel)

    (stage / "README.md").write_text(ANON_README, encoding="utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(stage: Path) -> None:
    rows = []
    for path in sorted(p for p in stage.rglob("*") if p.is_file()):
        rel = path.relative_to(stage).as_posix()
        if rel == "SHA256SUMS":
            continue
        rows.append(f"{sha256(path)}  {rel}")
    (stage / "SHA256SUMS").write_text("\n".join(rows) + "\n", encoding="utf-8")


def scan_identity(stage: Path) -> list[str]:
    failures: list[str] = []
    banned = ("zuizui0223", "ruiqi", "zhang.ruiqi")
    for path in sorted(p for p in stage.rglob("*") if p.is_file()):
        if should_preserve_raw(path.relative_to(stage)):
            continue
        try:
            text = path.read_text(encoding="utf-8").lower()
        except UnicodeDecodeError:
            continue
        for token in banned:
            if token in text:
                failures.append(f"{path.relative_to(stage)} contains {token!r}")
    return failures


def deterministic_zip(stage: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in sorted(p for p in stage.rglob("*") if p.is_file()):
            rel = path.relative_to(stage).as_posix()
            info = zipfile.ZipInfo(rel, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            zf.writestr(info, path.read_bytes())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--keep-stage", type=Path, default=None)
    args = ap.parse_args()

    tmp = None
    if args.keep_stage:
        stage = args.keep_stage
        if stage.exists():
            shutil.rmtree(stage)
        stage.mkdir(parents=True)
    else:
        tmp = tempfile.TemporaryDirectory(prefix="chocho-anonymous-review-")
        stage = Path(tmp.name)

    try:
        stage_snapshot(stage)
        failures = scan_identity(stage)
        if failures:
            raise RuntimeError("anonymous snapshot identity scan failed:\n" + "\n".join(failures))
        write_manifest(stage)
        deterministic_zip(stage, args.output)
        print(json.dumps({
            "output": str(args.output),
            "files": sum(1 for p in stage.rglob("*") if p.is_file()),
            "sha256": sha256(args.output),
        }, indent=2))
    finally:
        if tmp is not None:
            tmp.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

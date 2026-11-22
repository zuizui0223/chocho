#!/usr/bin/env python3
"""Build and immediately replay the internal-only Figure 2A-C source component."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

BASE = Path(__file__).resolve().parent
NAME = "GEB_V05_FIG2ABC_DERIVED_ONLY_NOT_FULL_REVIEW.zip"
AREA = Path("results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv")
ACTIVE = Path("results/three-gate-port-introduced-host-taxa/distinct_introduced_host_taxa_and_port_by_island.csv")
AREA_SHA = "ce1bf5d707fcf040e6ebe7e8dc69ccd5f0d09262f397bd27b39587e43030a3ce"
ACTIVE_SHA = "f1a4e543d479f6524457bb83eaf40c2075917f398161b07daa60794824644efc"
A_SOURCE = "verify_island_port_native_alien_matched_20261118.py"
MARINE_SOURCE = "audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117.py"
STRICT_HELPER = "audit_WCVP126_native_alien_double_unseen_islands_20261117.py"
ONE_SCRIPT = "audit_WCVP126_original_uniform_circuit_crossed_prediction_20261117.py"
STRICT_SCRIPT = "audit_WCVP126_uniform_circuit_double_unseen_20261118.py"
SCRIPT_SHA = {
    A_SOURCE: "5e08f7070a0dd41533813919e34afff5806b27c87c607abee54872963a1d2e78",
    MARINE_SOURCE: "46e738e512daa5f52c6a8b09bea1ae8cd4fe5d340bf6239dbe2cf703394793be",
    STRICT_HELPER: "c102d504a2e38f77ca1097804440ff8d6da2c934f1cd63eb34fe3870808f9b1f",
    ONE_SCRIPT: "aa9ad96aa7a73c48ba8f527bd203ca649b95dd2f5a17d4fda19fe9c9de7a5987",
    STRICT_SCRIPT: "83c1b042488f3eef639a5419e855eb4fd5606231978a26458d34896d49a853bd",
}
B_TABLE = "WCVP_126island_7875_native_alien_flora_similarity.csv"
PLANTS = "WCVP126_flora_similarity_strict_WPI2020_geographic_pairs.csv"
GEOMETRY = "WCVP126_GSHHG_source_coastal_anchors_and_sea_geometry_pairs.csv"
CIRCUIT = "WCVP126_original104_uniform_finvol_electrical_resistance.csv"
CIRCUIT_RECEIPT = "uniform104_WCVP_GSHHG_source_circuit_numerical_audit.json"
SEA_AUDIT = "GSHHG126_sea_geometry_source_coverage_audit.json"
JOIN_AUDIT = "WCVP126_geodesic_port_shipment_source_coverage_audit.json"
ONE_EXPECTED = "WCVP126_native_alien_true_uniform_circuit_vs_port_held_island_validation.json"
STRICT_EXPECTED = "WCVP126_104_true_uniform_circuit_strict_double_unseen.json"
ARTIFACTS = {
    14296292008: (44010011288, "9d3f4277e74f18e6b4935759caa131b93459401046b5c1ac21c7634d68980c87",
                 {B_TABLE: "f0cdefed6f6ff9a1bbd4f07df9efd1410c120abceda9c91e0fe67ff4d9d22010"}),
    14296307100: (44010307551, "c00508ea0384efe73841a82a7ad5e73d3df8a13eed36ecea23b4546238d0d2b4",
                 {PLANTS: "6e8b0d11eb63490aae782e6b8bab381aa07c07fe39123abcc24e9d77d302a799",
                  JOIN_AUDIT: "28937b6c90230e1c24c5ba807d33c457a78b3b7905c9ee59f479f8dc9ee1a221"}),
    14296891019: (44010673580, "36dbee9c52bebf2c8110fd894bd915f11b748c45057c3576dbd67df65bb8a8ff",
                 {GEOMETRY: "de38c3d55489787283c13b59862a7f86b5fcaad1340db51ed65434083d3d527d",
                  SEA_AUDIT: "5d9636fa16811fd5429780ea770d16b7341f519cd6b7e713074ba9a83d34e9f3"}),
    14298480144: (44013658955, "beefe9efb5f3bf8d200842218a1641d660592c8c379123ca60e0707bf967a1af",
                 {CIRCUIT: "e894a82f07eb7f5345669e6bf561f0779a7401a8406a39b5deefd3da587e1d6e",
                  CIRCUIT_RECEIPT: "8c496b7b2300305f8b62d07534bf3ac8bf3a0f0c7e94d02b6f9e57dfa7801d8b"}),
    14298086879: (44013825017, "16a64ea83d4c0ad9fe6f9961cdc303afc6fce79c65ef3c3efb02f1db44224252",
                 {ONE_EXPECTED: "a02264f01b6f0404e8111c7f9066e94fa050a246169d8f82cd16c7779107453f"}),
    14419709321: (44290887759, "cc427d36f673337b7226598d9b5e0f67551d4413f285e159504ac49cd95b9d62",
                 {STRICT_EXPECTED: "08bf98f9963f60dea3fcc4fb9140454ef1adeeccaef6fae0fb3bd61fc8ebb354"}),
}
APPENDIX_COVERAGE = {
    "S1": "Existing source/citation inventory. This add-on checks the distinct 126- and 104-island cohorts.",
    "S2": "Existing Fig1/3BC aggregate replays and separate Fig3A upstream rebuild receipt. Other S2 tests outside this add-on.",
    "S3": "This add-on refits Fig2A-C and the strict shortest-water comparator. Existing Fig2D/MEM10 component separate. "
          "EuDiS, Tonga, alternative metrics, botanical-effort adjustment and Figure S3.1 not run by this add-on.",
    "S4": "Same-source HOSTS×WCVP controls and documentation proxy follow their own pinned workflows.",
    "S5": "Existing Fig4 saved-score component and separate original224 refit/thinning evidence. "
          "155 model retraining and full SI S5 not run by this add-on.",
    "S6": "Endemic checklist and source-identified observations follow the independent qualitative/source validation.",
    "S7": "Historical records, port buffers, nutrient proxies and pathway-control tests outside this add-on.",
    "S8": "Separate rights, provenance and repository-access register. Full SI reproduction remains open.",
}
README = """# Internal Figure 2A–C derived-input replay

Run `python replay_fig2abc.py` after extracting this ZIP. Python 3.12 and NumPy
are sufficient. The runner performs no downloads. Fit results are compared with
archived original fold results and uncertainty intervals.

- A uses original aggregated WCVP native/introduced richness and WPI port
  counts on 126 botanical island regions, retaining the 123 active-region,
  91 port-positive and nine leave-one-Level1-out reruns.
- B recomputes Jaccard from original taxon intersection and union counts.
  All 7,875 native and 7,874 introduced dyads have defined Jaccard. The
  introduced BOU–SSA zero-union dyad remains undefined.
- C uses 5,355 defined botanical dyads on 104 islands, conditionally known
  floristic margins, and derived geographic features. It reruns 104
  one-island-held-out folds, the exact five-fold 1,089-dyad both-islands-unseen
  fit, and the same-cohort strict shortest-water comparator. The source
  circuit is undirected electrical geometry. Its descriptive intervals are
  fold-resampling summaries.

These projected tables contain botanical-region codes, aggregated richness
and port counts, derived dyad features, and no coordinates of WPI ports,
original individual records, GBIF occurrence IDs, or AIS measures. The inputs
do not rebuild the WCVP flora incidence, reconstruct the original WPI conversion,
or solve the upstream GSHHG circuit graph. The study's original converted
WPI GeoJSON edition and original-to-derived rights chain remain unresolved.

Original datasets include WCVP v13 (DOI 10.34885/0yex-xv26, CC BY 3.0),
GSHHG v2.3.7 (DOI 10.5281/zenodo.7007502), and WPI Publication 150.
See MANIFEST.json for the source chain and eight-appendix coverage. The
component is for internal verification. Reviewer anonymity, source
redistribution clearance, whole-SI reproduction and stable reviewer access
remain open.
"""
ENV = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
       "PYTHONDONTWRITEBYTECODE": "1"}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def require_hash(data, expected, name):
    if digest(data) != expected:
        raise ValueError(f"SHA256 mismatch: {name}")
    return data


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_projection(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def json_text(data):
    return json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def collect_sources(artifacts, source_dir):
    for artifact, (_, archive_sha, members) in ARTIFACTS.items():
        path = artifacts / f"{artifact}.zip"
        require_hash(path.read_bytes(), archive_sha, path.name)
        with ZipFile(path) as archive:
            for name, sha in members.items():
                found = [p for p in archive.namelist() if PurePosixPath(p).name == name]
                if len(found) != 1:
                    raise ValueError(f"Expected one {name} in {path.name}")
                data = require_hash(archive.read(found[0]), sha, name)
                (source_dir / name).write_bytes(data)


def projected_a(area_file, active_file, destination, fields):
    active = {row["code"] for row in read_csv(active_file)}
    rows = read_csv(area_file)
    if len(active) != 123 or not active <= {row["code"] for row in rows}:
        raise ValueError("Active source WCVP/HOSTS region identities changed")
    output = []
    for row in rows:
        reduced = {key: row[key] for key in fields if key != "resource_active"}
        reduced["resource_active"] = "1" if row["code"] in active else "0"
        output.append(reduced)
    write_projection(destination, fields, output)


def projected_b(path, destination, fields):
    projection = []
    for row in read_csv(path):
        result = {key: row[key] for key in ("source_region", "destination_region")}
        for target, source in (("native", "native"), ("introduced", "alien")):
            shared = int(row[f"{source}_shared_accepted_taxa"])
            union = int(row[f"{source}_union_accepted_taxa"])
            stored = row[f"{source}_Jaccard"]
            if union:
                if not math.isclose(float(stored), shared / union, rel_tol=1e-13, abs_tol=1e-15):
                    raise ValueError("Stored Jaccard differs from original shared/union counts")
            elif shared or stored not in ("", "NA"):
                raise ValueError("Unexpected value for undefined Jaccard")
            result.update({f"{target}_shared": shared, f"{target}_union": union})
        projection.append(result)
    write_projection(destination, fields, projection)


def projected_c(source, destination, fields, original):
    paths = SimpleNamespace(plants=source / PLANTS, geometry=source / GEOMETRY)
    pairs, undefined = original.input_pair(paths)
    if len(pairs) != 5355 or len(undefined) != 1:
        raise ValueError("Expected original marine support and one undefined pair")
    excluded = undefined[0]
    if (excluded["region_A"], excluded["region_B"], excluded["native_Jaccard_defined"],
        excluded["alien_Jaccard_defined"]) != ("BOU", "SSA", True, False):
        raise ValueError("Marine zero-union pair no longer matches standalone WCVP source")
    circuit = {}
    for row in read_csv(source / CIRCUIT):
        key = tuple(sorted((row["source_region"], row["destination_region"])))
        if key in circuit:
            raise ValueError("Duplicate original resistance dyad")
        circuit[key] = row["log1p_normalized_uniform_circuit"]
    if len(circuit) != 5356 or ("BOU", "SSA") not in circuit:
        raise ValueError("Expected original 5,356 solved resistance dyads")
    output = []
    for row in pairs:
        reduced = {key: row[key] for key in fields if key != "circuit"}
        reduced["circuit"] = circuit[tuple(sorted((row["island1"], row["island2"])))]
        output.append(reduced)
    write_projection(destination, fields, output)


def run_original(source, one):
    key = ONE_SCRIPT if one else STRICT_SCRIPT
    out = source / ("one_fit" if one else "strict_fit")
    args = ["--plants", str(source / PLANTS), "--uniform-circuit", str(source / CIRCUIT),
            "--island-area", str(source / AREA), "--out-dir", str(out)]
    if one:
        args += ["--shortest-water", str(source / GEOMETRY),
                 "--uniform-circuit-receipt", str(source / CIRCUIT_RECEIPT),
                 "--original-sea-audit", str(source / SEA_AUDIT),
                 "--original-shipping-audit", str(source / JOIN_AUDIT)]
    else:
        args += ["--geometry", str(source / GEOMETRY), "--circuit-receipt", str(source / CIRCUIT_RECEIPT)]
    subprocess.run([sys.executable, str(source / "scripts" / key), *args], env=ENV,
                   check=True, capture_output=True, text=True, timeout=120)
    result = ONE_EXPECTED if one else STRICT_EXPECTED
    return json.loads((out / result).read_bytes()), json.loads((source / result).read_bytes())


def choose(document, keys):
    return {key: document[key] for key in keys}


def build(args):
    with tempfile.TemporaryDirectory(prefix="fig2abc-build-") as work:
        work = Path(work)
        source = work / "source"
        component = work / "component"
        (source / "scripts").mkdir(parents=True)
        (component / "data").mkdir(parents=True)
        for name, base in ((AREA, AREA_SHA), (ACTIVE, ACTIVE_SHA)):
            destination = source / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(require_hash((args.audit_root / name).read_bytes(), base, str(name)))
        collect_sources(args.artifacts, source)
        for name, sha in SCRIPT_SHA.items():
            path = args.one_held_script if name == ONE_SCRIPT else args.audit_root / "scripts" / name
            data = require_hash(path.read_bytes(), sha, name)
            (source / "scripts" / name).write_bytes(data)
            if name != ONE_SCRIPT:
                (component / name).write_bytes(data)
        sys.path.insert(0, str(source / "scripts"))
        sys.dont_write_bytecode = True
        import audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 as original
        import replay_fig2abc as replay

        projected_a(source / AREA, source / ACTIVE, component / "data/port_flora_regions.csv", replay.PORT_FIELDS)
        projected_b(source / B_TABLE, component / "data/flora_jaccard_counts.csv", replay.FLORA_FIELDS)
        projected_c(source, component / "data/marine_fit_inputs.csv", replay.MARINE_FIELDS, original)

        fitted_a = subprocess.run([sys.executable, str(source / "scripts" / A_SOURCE)], env=ENV,
                                  check=True, capture_output=True, text=True, timeout=30)
        source_a = json.loads(fitted_a.stdout)
        expected_a = choose(source_a, ("model_scopes", "geographic_sensitivity"))
        del expected_a["geographic_sensitivity"]["interpretation"]
        expected_b = {
            "source_island_pairs": 7875,
            "defined_pairs": {"native": 7875, "introduced": 7874},
            "undefined_introduced_pair": ["BOU", "SSA"],
            "median_Jaccard": {"native": 0.00448654037886341, "introduced": 0.0460304551145545},
        }
        actual_one, archived_one = run_original(source, one=True)
        actual_strict, archived_strict = run_original(source, one=False)
        one_keys = ("native_flora", "introduced_flora")
        strict_keys = (
            "source_common_regions", "source_common_pairs", "source_n_test_double_unseen_pairs",
            "n_source_folds", "original_shortest_double_unseen_verified_before_electrical_fit",
            "native_uniform_circuit", "alien_uniform_circuit", "origin_paired_comparison",
        )
        # Check complete analytical blocks, including actual heldout losses and intervals.
        n_one = replay.compare(choose(actual_one, one_keys), choose(archived_one, one_keys), "one-held original")
        n_strict = replay.compare(choose(actual_strict, strict_keys), choose(archived_strict, strict_keys), "strict original")
        expected = {"A": expected_a, "B": expected_b,
                    "C": {"one_held": choose(archived_one, one_keys),
                          "strict": choose(archived_strict, strict_keys)}}
        (component / "EXPECTED_FIT_RESULTS.json").write_text(json_text(expected), encoding="utf-8")
        (component / "README.md").write_text(README, encoding="utf-8")
        shutil.copyfile(BASE / "replay_fig2abc.py", component / "replay_fig2abc.py")
        manifest = {
            "schema": "GEB_v05_Fig2ABC_derived_replay_v1",
            "status": "INTERNAL_ONLY",
            "release_ready": False,
            "source_licenses_cleared": False,
            "indirect_anonymity_verified": False,
            "stable_reviewer_access_available": False,
            "full_SI_replayed": False,
            "appendix_coverage": APPENDIX_COVERAGE,
            "source_git_snapshots": {
                "one_held": "8f9b37ac9ef9f18727fefcee7a1fc61360b4d503",
                "strict": "9a718997bce21f5b99329a91844bd31c375c6d16",
                "A_aggregate_source": "be3bf40a59023b2336da3876d50fd2c9e8ca76ae",
            },
            "source_code_sha256": SCRIPT_SHA,
            "source_aggregate_sha256": {str(AREA): AREA_SHA, str(ACTIVE): ACTIVE_SHA},
            "source_artifacts": {
                str(artifact): {"run": run, "archive_sha256": sha, "member_sha256": members}
                for artifact, (run, sha, members) in ARTIFACTS.items()
            },
            "original_shipped_result_fields_checked": {"one_held": n_one, "strict": n_strict},
            "archive_file_sha256": {
                path.relative_to(component).as_posix(): digest(path.read_bytes())
                for path in sorted(component.rglob("*")) if path.is_file()
            },
        }
        (component / "MANIFEST.json").write_text(json_text(manifest), encoding="utf-8")
        # Capture actual leaf-level results from the ZIP's standalone offline runner.
        archive_path = work / NAME
        forbidden = re.compile(
            rb"zuizui0223|zhang[ _.-]?(?:ruiqi|rachel)|[a-z0-9_.+-]+@[a-z0-9.-]+\.[a-z]{2,}"
            rb"|/(?:home|Users)/[a-z0-9_.-]+/", re.I
        )
        with ZipFile(archive_path, "w") as archive:
            for path in sorted(component.rglob("*")):
                if not path.is_file():
                    continue
                name = path.relative_to(component).as_posix()
                data = path.read_bytes()
                if forbidden.search(name.encode()) or forbidden.search(data):
                    raise ValueError(f"Direct source-identification token in {name}")
                item = ZipInfo(name, (2026, 11, 22, 0, 0, 0))
                item.external_attr = 0o644 << 16
                item.compress_type = ZIP_DEFLATED
                archive.writestr(item, data)
        extracted = work / "extracted"
        with ZipFile(archive_path) as archive:
            if archive.testzip() is not None:
                raise ValueError("Component ZIP CRC failed")
            archive.extractall(extracted)
        replay_run = subprocess.run([sys.executable, str(extracted / "replay_fig2abc.py")],
                                    cwd=extracted, env=ENV, check=True, capture_output=True, text=True, timeout=120)
        result = json.loads(replay_run.stdout)
        result["archive_sha256"] = digest(archive_path.read_bytes())
        result["archive_bytes"] = archive_path.stat().st_size
        result["archive_files"] = len(manifest["archive_file_sha256"]) + 1
        args.out.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(archive_path, args.out / NAME)
        (args.out / "FIG2ABC_COMPONENT_QA.json").write_text(json_text(result), encoding="utf-8")
        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", type=Path, required=True, help="Directory of original numeric-ID ZIPs")
    parser.add_argument("--one-held-script", type=Path, required=True,
                        help="Original one-island-fit script at 8f9b37ac9ef9f18727fefcee7a1fc61360b4d503")
    parser.add_argument("--audit-root", type=Path, required=True, help="Checkout containing the original aggregate CSVs")
    parser.add_argument("--out", type=Path, required=True)
    print(json_text(build(parser.parse_args())), end="")

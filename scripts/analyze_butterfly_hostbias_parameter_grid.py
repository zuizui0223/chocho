#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-protocol-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=299)
    args = ap.parse_args()

    base = json.loads(args.base_protocol_json.read_text(encoding="utf-8"))
    bandwidths = [0.1, 0.2, 0.4]
    exponents = [0.5, 1.0, 2.0]
    rows = []

    target_script = Path(__file__).with_name(
        "analyze_butterfly_resource_expansion_hostbias_null.py"
    )
    with tempfile.TemporaryDirectory(prefix="hostbias-grid-") as tmpdir:
        tmp = Path(tmpdir)
        for bandwidth in bandwidths:
            for exponent in exponents:
                protocol = dict(base)
                protocol["native_log_bandwidth"] = bandwidth
                protocol["usage_exponent"] = exponent
                protocol["permutations"] = args.permutations
                protocol["seed_tag"] = (
                    f"chocho-geb-hostbias-grid-v01-bw{bandwidth}-exp{exponent}"
                )
                protocol_path = tmp / f"protocol_bw{bandwidth}_exp{exponent}.json"
                result_path = tmp / f"result_bw{bandwidth}_exp{exponent}.json"
                protocol_path.write_text(
                    json.dumps(protocol, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                subprocess.run(
                    [
                        sys.executable,
                        str(target_script),
                        "--protocol-json",
                        str(protocol_path),
                        "--descriptors-csv",
                        str(args.descriptors_csv),
                        "--insect-host-csv",
                        str(args.insect_host_csv),
                        "--native-distribution-csv",
                        str(args.native_distribution_csv),
                        "--contemporary-distribution-csv",
                        str(args.contemporary_distribution_csv),
                        "--output-json",
                        str(result_path),
                    ],
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
                result = json.loads(result_path.read_text(encoding="utf-8"))
                model = result["models"]["native_range_usage"]
                total = model["total_introduced_added_units"]
                mean = model["mean_log_expansion"]
                median = model["median_log_expansion"]
                diag = model["matching_diagnostics"]
                observed = result["observed"]
                rows.append(
                    {
                        "native_log_bandwidth": bandwidth,
                        "usage_exponent": exponent,
                        "species": result["species"],
                        "observed_total_added": observed[
                            "total_introduced_added_units"
                        ],
                        "null_total_q025": total["q025"],
                        "null_total_median": total["median"],
                        "null_total_q975": total["q975"],
                        "p_total_high": total["p_high"],
                        "observed_mean_log_expansion": observed[
                            "mean_log_expansion"
                        ],
                        "null_mean_log_median": mean["median"],
                        "p_mean_high": mean["p_high"],
                        "observed_median_log_expansion": observed[
                            "median_log_expansion"
                        ],
                        "null_median_log_median": median["median"],
                        "p_median_high": median["p_high"],
                        "median_native_match_difference": diag[
                            "median_absolute_log1p_native_breadth_difference_across_replicates"
                        ],
                        "mean_sampled_other_users": diag[
                            "mean_other_lepidoptera_users_per_sampled_host_median"
                        ],
                    }
                )

    payload = {
        "schema": "chocho_butterfly_hostbias_parameter_grid_v0.1",
        "status": "POSTHOC_PARAMETER_ROBUSTNESS",
        "permutations_per_cell": args.permutations,
        "bandwidths": bandwidths,
        "usage_exponents": exponents,
        "grid": rows,
        "summary": {
            "cells": len(rows),
            "cells_total_p_lt_0_05": sum(r["p_total_high"] < 0.05 for r in rows),
            "cells_mean_p_lt_0_05": sum(r["p_mean_high"] < 0.05 for r in rows),
            "cells_median_p_lt_0_05": sum(
                r["p_median_high"] < 0.05 for r in rows
            ),
            "cells_observed_total_above_null_q975": sum(
                r["observed_total_added"] > r["null_total_q975"] for r in rows
            ),
            "cells_observed_total_within_null_95": sum(
                r["null_total_q025"]
                <= r["observed_total_added"]
                <= r["null_total_q975"]
                for r in rows
            ),
            "minimum_total_p_high": min(r["p_total_high"] for r in rows),
            "maximum_total_p_high": max(r["p_total_high"] for r in rows),
        },
        "interpretation_rule": (
            "The attenuation of the simple host-identity excess is parameter-robust "
            "if observed expansion remains inside the combined-null distribution "
            "across reasonable native-range bandwidths and other-Lepidoptera-use "
            "weight exponents."
        ),
        "claim_boundary": (
            "This grid probes reasonable tuning choices but does not make HOSTS "
            "consumer degree a pure measure of biological host commonness."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

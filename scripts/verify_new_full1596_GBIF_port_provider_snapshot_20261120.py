#!/usr/bin/env python3
"""Legacy alias for archived NEW GBIF 1596 WPI port attribution QA.

The original first exported v1 JSON was superseded by source-equivalent v2
full-provider archive with SHA256 e79c3ad4..., produced only from the four
previously successful 399-port source shards. Keep this old entry point
read-only and forward to the canonical v2 verification script.

Neither version retroactively recovers original 2026 GBIF datasetKey sources,
reproduces the historic WKT port-buffer exact occurrence roster, or mints DOI.
"""
from verify_GEB_v05_new_port1596_attribution_20261120 import run

if __name__ == "__main__":
    run()
    print("LEGACY_GBIF1596_VERIFY_ALIAS_TO_CANONICAL_CURRENT_SOURCE_V2_PASS")

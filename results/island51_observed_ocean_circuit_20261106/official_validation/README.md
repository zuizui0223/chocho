# Reproduce the official Circuitscape numerical checks

The primary Python run computes three complete 51-island distributed-port circuit matrices. The independent official check uses **Circuitscape 6.1.0, Julia 1.13.3, double precision and CHOLMOD** on two fixed pairs (`9–15`, `715–1577`) for each complete exported marine network. It does not compute a second full 51-island matrix.

The original successful execution is preserved byte for byte under `original_execution/`. A portable entry point at [`scripts/verification/validate_island51_observed_ocean_circuits.jl`](../../../../scripts/verification/validate_island51_observed_ocean_circuits.jl) regenerates source, ground and Circuitscape configuration files from the exported node weights. It accepts paths at runtime and requires no machine-specific source path.

## Run the six checks

From the repository root, set `JULIA_BINARY` to a Julia **1.13.3** executable and `CS_CHECKOUT` to an installed Circuitscape **6.1.0** checkout. The installation section below provides the exact source pin and dependency snapshot.

If the Python results and network exports already exist, this is the only computation required:

```bash
"$JULIA_BINARY" --project="$CS_CHECKOUT" \
  scripts/verification/validate_island51_observed_ocean_circuits.jl \
  /tmp/chocho-island51-observed-ocean-networks \
  /tmp/chocho-island51-observed-ocean-results \
  /tmp/chocho-island51-official-checks
```

Arguments are `NETWORK_DIR RESULTS_DIR OUT_DIR`. `NETWORK_DIR` must contain `uniform_conductances.txt`, `ais_conductances.txt`, `wind_conductances.txt` and `island_node_weights.csv`. `RESULTS_DIR` must contain the corresponding `island51_observed_circuit_pairs.csv`. The graph exports contain each undirected edge once, with compact one-based node IDs and conductance in the third column.

To generate those inputs from the frozen repository bundle first:

```bash
python -u scripts/run_island51_observed_ocean_circuits.py \
  --input-root data/exploratory/island51_observed_ocean_20261106 \
  --coast-dir data/exploratory/island51_observed_ocean_20261106/coast_0p0625 \
  --marine-ports-only \
  --out-dir /tmp/chocho-island51-observed-ocean-results \
  --network-dir /tmp/chocho-island51-observed-ocean-networks
```

The portable Julia check writes six rows to `official_checks.csv`, source/ground/configuration files, and a `verification_summary.toml` only after all checks pass. Large voltage or current maps are disabled. It checks positive finite resistance, agreement with the reference pair table, current conservation at every graph node, and equality of branch dissipation with `b'v`. The independent numerical tolerance is `1e-7`; the Circuitscape solver option is `residual_tolerance=1e-9`.

The real rerun using this portable entry point is retained under `portable_reexecution/` with its log. These are small numerical verification outputs, not new biological fits.

## Install the pinned official environment

The official package checkout is pinned to:

- Repository: <https://github.com/Circuitscape/Circuitscape.jl>
- Tag: `v6.1.0`
- Commit: `253a4189296a01158b8e38f48b9fcb1745815161`
- Julia: `1.13.3`
- Official Linux x86-64 archive: <https://julialang-s3.julialang.org/bin/linux/x64/1.13/julia-1.13.3-linux-x86_64.tar.gz>
- Julia archive SHA256: `5c61cb4b40332cb5a3504e7779e2f5dac093caa4a1116ded19f01413f85c63e4`

`environment/Project.toml` and `environment/Manifest.toml` are the exact successful package-project snapshots. The manifest identifies Circuitscape itself with `path = "."`; consequently it belongs **inside the pinned source checkout**. The archived environment directory alone is not a standalone installed Circuitscape package.

```bash
JULIA_BINARY=/path/to/julia-1.13.3/bin/julia
CS_CHECKOUT=/tmp/chocho-circuitscape-v6.1.0
export JULIA_DEPOT_PATH=/tmp/chocho-circuitscape-depot
export JULIA_NUM_THREADS=2

git clone --branch v6.1.0 --depth 1 \
  https://github.com/Circuitscape/Circuitscape.jl "$CS_CHECKOUT"
test "$(git -C "$CS_CHECKOUT" rev-parse HEAD)" = \
  253a4189296a01158b8e38f48b9fcb1745815161
cp results/island51_observed_ocean_circuit_20261106/official_validation/environment/Project.toml \
  "$CS_CHECKOUT/Project.toml"
cp results/island51_observed_ocean_circuit_20261106/official_validation/environment/Manifest.toml \
  "$CS_CHECKOUT/Manifest.toml"
"$JULIA_BINARY" --project="$CS_CHECKOUT" -e 'using Pkg; Pkg.instantiate()'
```

The exact Julia manifest SHA256 is `0be996454077b0de717e4e3ee8d6934aff2c5c896c04fd3d9354bed7d231f1a2`.

## Distributed current and ground convention

For each pair, aggregate the two islands' port weights to a signed node injection `b=w_a-w_b`. Duplicate node rows add. Each island's positive weights must sum to one. No terminal ports are short-circuited.

Circuitscape 6.1.0 advanced mode skips a component when its explicit source strengths sum exactly to zero. The script therefore chooses one negative-injection node as a direct ground and omits only that node's source row. Its omitted sink is supplied by current conservation at the ground; all other negative and positive source rows remain. The output measure is `b'v`, with zero ground voltage. Independent all-node Kirchhoff and branch-energy checks verify the representation.

The package-generated configuration includes a legacy `[Version] version = 5.0.0` boilerplate value. Actual runtime versions come from `VERSION` and `pkgversion(Circuitscape)` and are recorded in the receipts and portable summary. The script requires the pinned versions explicitly.

## Tests and preserved evidence

The portable boundary tests invoke real Circuitscape on a four-node chain. They check six hand-derived values, duplicate node-weight aggregation and rejection of an intentionally wrong reference matrix value:

```bash
JULIA_BINARY="$JULIA_BINARY" CIRCUITSCAPE_PROJECT="$CS_CHECKOUT" \
  python scripts/verification/test_validate_island51_observed_ocean_circuits.py
```

The original series, ring and distributed-chain fixture script, all seven assertion results, successful six-case script, full injection definitions, configurations, raw checks and logs are retained under `original_execution/`. Original machine paths remain in that immutable evidence; use the portable entry point for reexecution. `evidence_manifest.json` inventories the small evidence files and current portable scripts by SHA256.

These checks establish numerical agreement between two solver implementations for the represented networks and terminal distributions. They do not validate source AIS units, demonstrate ship movements, reconstruct directional winds or establish causal/historical biological transport.

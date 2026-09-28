# Script map

Scripts are grouped by role so the repository root of `scripts/` stays navigable.

## `paper/` — manuscript, review and release

- `render_butterfly_specialization_manuscript_figures.py` — regenerate the five manuscript figures from frozen figure-source inputs.
- `render_butterfly_specialization_blinded_manuscript.py` — render the blinded Markdown source.
- `build_blinded_review_docx.py` — build the editable double-anonymous review DOCX.
- `build_anonymous_review_bundle.py` — build the de-identified reviewer code bundle.
- `release_preflight.py` — report remaining administrative blockers for v1.

## `resource/` — host-resource reconstruction and ecological analyses

- `build_butterfly_resource_envelope_pilot.py` — reconstruct the frozen S1 resource descriptors.
- `build_butterfly_contemporary_resource_envelope.py` — reconstruct contemporary host-resource geography.
- `analyze_butterfly_anthropogenic_resource_expansion.py` — native-to-contemporary resource expansion.
- `analyze_butterfly_resource_expansion_mechanism.py` — host-contribution decomposition.
- `analyze_butterfly_host_breadth_geography.py` — taxonomic versus geographic specialization.
- `analyze_butterfly_host_specialization_hierarchy.py` — within-family portfolio hierarchy.
- `build_wcvp_hosts_sidecar.R` and `build_wcvp_hosts_contemporary_sidecar.R` — base-R WCVP/HOSTS sidecar builders.

## `climate/` — independent climate route

- `freeze_butterfly_climate_release_panel.py`
- `acquire_butterfly_resource_envelope_occurrences.py`
- `apply_butterfly_climate_release_quality_gate.py`
- `analyze_butterfly_resource_envelope_climate.py`
- `test_butterfly_climate_release_hypothesis.py`

These preserve the frozen independent panel, quality gates, cross-fit and primary-test boundary.

## `diagnostics/` — mapping and historical diagnostics

- `map_butterfly_resource_envelope_wgsrpd3.py`
- `compare_butterfly_native_contemporary_occurrence_overlap.py`
- `diagnose_butterfly_resource_sampling_effort.py`
- `apply_butterfly_resource_envelope_climate_gate.py`

Diagnostics are retained for auditability. They are not additional manuscript analyses and must not be used for response-driven model search.

For claim boundaries use `manuscript/butterfly_specialization_claim_map_v0.1.json`.

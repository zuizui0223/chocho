# Script map

The repository keeps executable history, but only a small set of scripts are normal entry points for the current paper.

## Submission and release

- `render_butterfly_specialization_v02_figures.py` — regenerate the current three main figures plus Supplementary Figure S1 from vendored figure sources and frozen result receipts.
- `build_blinded_review_docx.py` — build the editable double-anonymous review DOCX.
- `build_anonymous_review_bundle.py` — build the de-identified reviewer code bundle.
- `release_preflight.py` — report the remaining administrative blockers for the v1 release.

## Resource reconstruction and ecological analyses

- `build_butterfly_resource_envelope_pilot.py` — reconstruct the frozen S1 butterfly resource descriptors.
- `build_butterfly_contemporary_resource_envelope.py` — reconstruct contemporary host-resource geography.
- `analyze_butterfly_anthropogenic_resource_expansion.py` — native-to-contemporary resource expansion.
- `analyze_butterfly_resource_expansion_mechanism.py` — historical within-butterfly host-contribution decomposition.\n- `analyze_butterfly_host_contribution_concentration.py` — current unit-preserving across-plant species/genus contribution concentration.\n- `analyze_butterfly_expansion_equivalence.py` — current bootstrap precision and strict post-hoc diet-breadth equivalence diagnostic.
- `analyze_butterfly_poaceae_sensitivity.py` — exclude all Poaceae users or one-family Poaceae specialists to test whether redistributed grasses create the diet-breadth result.
- `analyze_butterfly_host_breadth_geography.py` — taxonomic versus geographic specialization.
- `analyze_butterfly_host_specialization_hierarchy.py` — within-family portfolio hierarchy.
- `build_wcvp_hosts_sidecar.R` and `build_wcvp_hosts_contemporary_sidecar.R` — base-R WCVP/HOSTS sidecar builders.

## Independent climate route

- `freeze_butterfly_climate_release_panel.py`
- `acquire_butterfly_resource_envelope_occurrences.py`
- `apply_butterfly_climate_release_quality_gate.py`
- `analyze_butterfly_resource_envelope_climate.py`
- `test_butterfly_climate_release_hypothesis.py`

The climate route preserves the frozen panel, quality gates and independent-test boundary.

## Diagnostics and historical execution helpers

The remaining mapping, overlap, sampling-effort and pilot-gate scripts are retained because they document the route to the frozen result. They are not additional manuscript analyses and should not be used for response-driven model search.

For the exact current manuscript claim boundaries, use `manuscript/butterfly_specialization_claim_map_v0.2.json`.

# GEB original155 — complete source-input-to-model-panel reproduction (24 November 2026)

**PASS in the checked scope.** The original archived `build_forecast_panel.py` was executed **unchanged** (source SHA-256 `99af0755252f2bb81c53f1b07a504e34d15092ccbd8a10a0bc16dde0e0b57ec6`) from the recovered 13 originally expected source-input paths. The first original private Draft-release archive had 10/13 panel-builder inputs; the missing candidate selection, family support and origin host-count CSVs were separately source-byte and receipt-verified in the newer 24-Nov Draft-release addendum (not silently treated as original 21-Nov archive members).

[Executed source-builder CI #45318106188](https://github.com/zuizui0223/chocho/actions/runs/45318106188) **success**, [source receipt artifact #14886548204](https://github.com/zuizui0223/chocho/actions/runs/45318106188/artifacts/14886548204). This is a real source-code execution and all-column comparison, **not** a score-only replay.

| Validation level | Outcome |
|---|---|
| Reconstructed keyed `speciesKey × cell × period` | **842,964 / 842,964**, all keys matched, zero duplicates |
| Columns in original and rebuilt panels | **96 / 96**, all numerical/categorical values matched within stated tolerance |
| `native_host_count` | **0 differences**, maximum absolute difference **0.0** |
| `introduced_host_count` | **0 differences**, maximum absolute difference **0.0** |
| `introduced_share` | **0 differences**, maximum absolute difference **0.0** |
| Full-panel all-column mismatches | **0** |
| Preset comparison tolerance | **RTOL 1e−6; ATOL 1e−8** |
| Original saved panel SHA-256 | `5ae5b48b752fad05067537de5fb95d4ee74f17f3ba1114443d337bd7fd3695a1` |
| Fresh output panel SHA-256 | `52e4126fe8df941011726b7f6d9ba28fe45cf180330a0b9929c4a4b5ef698a97` |

**Critical distinction:** the regenerated Parquet's *bytes* differ from the original, while all compared values match (not a byte-identical recreation). The late-added `host_region_counts.csv` covers 87,235 species×regions and contains source-missing introduced shares, which are handled by the original script's `fillna(0)` rule after the join. It was not reverse-filled from the final panel by this run. However, recovering an input on 24 November and finding it produces matching output **does not independently prove when or from which original upstream data that CSV itself was created**.

No raw GBIF occurrence partitions were reprocessed in this job; saved species/provider-day counts, past occurrence cells, environmental source intermediates and host-composition tables were used. Consequently the claim **raw original GBIF→day preprocessing→all 96-column panel** is still open. The earlier 93-column partial replay and 10-of-13 inventory remain authentic historical checkpoints, superseded for the *current ability to run the full panel builder* only.

The separate original155 temporal+spatial model refits succeeded; unseen-species fold 3 numerical reproduction and two additional original temporal sensitivity refits remain mismatched. This step neither fixes nor conceals them. Nothing was uploaded with private observation rows or author details. Neither a permanent reviewer-accessible repository nor a redistribution licence for the full source bundle has been established. Stage A, Draft PRs #179/#189, no coauthor requests or journal submission.

[Machine-readable QA record](../../../results/geb_v05_private_source_recovery_20261120/ORIGINAL155_ALL96_PANEL_BUILDER_REPRO_VERIFIED_20261124.json).

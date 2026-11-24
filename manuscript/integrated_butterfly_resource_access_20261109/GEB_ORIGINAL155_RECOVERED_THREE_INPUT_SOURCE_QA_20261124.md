# 24 November: source verification of the newly found original155 input files

**Verified actual recovery, not merely file names.** Two assets were added on 24 November 2026 to the previously existing **Draft** release `temporal-analysis-recovery-20261121` (ID 453667492): ZIP **726,857 bytes** and JSON source receipt **1,302 bytes**. The original four-part 5.55-GB release remains independently pinned and unchanged.

[New three-input source verification CI #45317873271](https://github.com/zuizui0223/chocho/actions/runs/45317873271) **PASS**: the private 2026-11-24 ZIP and accompanying JSON were each stream-verified against GitHub's original asset size and SHA-256; all required CSV members were opened, their full bytes hashed and each hash located in the separate JSON receipt. Row counts, expected required columns, selected taxa/families and basic origin-count/share value ranges were validated. No source rows, direct project-local private paths or CSV bodies were exported.

| Newly recovered input | Verified rows | Special checks | SHA-256 |
|---|---:|---|---|
| `outputs/forecast_candidate_disposition.csv` | **225** | **155** unique eligible butterfly taxon IDs | `f275181f49a7d8da99cae815717c2f2ce608eb14ca3be5609482d82f307a7003` |
| `outputs/temporal_host_family_support.csv` | **109** | Thresholds select **18** host-plant families | `4644a2ed2f233052089b8880c25f22072158421f5f0396d668398899414cbb0e` |
| `work/remote_sensing/host_region_counts.csv` | **87,235** | All 87,235 species×region keys unique; 42,418 source `introduced_share` values missing | `5c074721f4ebb3e2e413d6e4926320b97f5cb1b78f7004b0138ceeecd33d2b8d` |

The original panel-builder fills missing `introduced_share` with zero **after the source join**, so 42,418 missing source values are not alone a failure. The added ZIP has **four non-directory members**: the three required tables plus one additional member; the independent check examines these three specified CSVs, not a claim that the archive contains only three.

**Identified provenance boundary:** The original 21-Nov 3,032-file archive did **not** include these exact three inputs. Their addition on 24 November and agreement with final derived data cannot, by itself, prove the 2026-10-29 original source timestamps or independence from reconstruction using the previously frozen panel. We therefore keep the historic source-lineage issue distinct from input content integrity.

The frozen original code `analysis_source/build_forecast_panel.py` has SHA-256 `99af0755252f2bb81c53f1b07a504e34d15092ccbd8a10a0bc16dde0e0b57ec6`. A separate [full source-conditional panel regeneration CI #45318106188](https://github.com/zuizui0223/chocho/actions/runs/45318106188) was initiated with **all 13** source inputs and exact original code. The new three-source proof above does **not** claim the complete 842,964-row original panel has been regenerated at the time of this receipt.

This finding does not update the four-figure manuscript statistics or settle the separate unseen-species fold 3 and two temporal sensitivity numerical discrepancies; third-party redistribution and anonymous reviewer access remain separate open Stage A gates. The recovered source release remains **Draft**.

# GEB v0.5 Figure 2A–C and four-component replay verification

**2026-11-22 · Internal source audit for Draft PR #189**

Figure 2A–C now runs alongside the existing Figure 1, Figure 2D, Figure 3B–C and Figure 4 saved-score checks. GitHub Actions rebuilt the 39-file package and ran all four components after extraction at audit commit `81b869e37a9ca2ec711c99d775c6d1d34d9caf5e`. An independent local build produced the same inner ZIP SHA-256. The run ended with `GEB_V05_LIMITED_FOUR_COMPONENT_REPLAY_PASS_NOT_FULL_PAPER`. [CI run 44970859551](https://github.com/zuizui0223/chocho/actions/runs/44970859551) · [Job 137050100625](https://github.com/zuizui0223/chocho/actions/runs/44970859551/job/137050100625)

| Archive | Bytes | SHA-256 |
| --- | ---: | --- |
| Actions delivery archive, artifact 14733606979 | 563,382 | `3a4beb33c72093aaf146ecca8da28650e2264fa221f7f9f8ab12bf8c9abbd09d` |
| Inner 39-file composite, `GEB_V05_INTERNAL_COMPOSITE_REPRO_NOT_FOR_JOURNAL_REVIEW.zip` | 560,519 | `40369ba5820f0e40ed56758f64165f951c82c6d3d64e5613175533d3c2055c95` |
| Standalone Figure 2A–C build input, 11 members | 448,525 | `0d34e188fb6ee2184a00f54853b61d49dfa085008679e6c537bae2e9386357d3` |

Artifact retention ends **2027-02-20 01:18:47 UTC**. The internal link uses repository access and exposes the Actions delivery archive. [Artifact 14733606979](https://github.com/zuizui0223/chocho/actions/runs/44970859551/artifacts/14733606979)

## Figure 2A–C evidence

The builder checks six original Actions archives and each consumed member against pinned SHA-256 values. It also checks five original scripts and two aggregate input tables. It projects only the columns needed for the three replays. The extracted Figure 2A–C runner verifies its packaged members and works offline with Python 3.12 and NumPy. The full composite runner also needs internet for dependencies and pinned TDWG geography.

The replay reports **1,947 scalar comparisons**. Numerical comparisons use relative tolerance `1e-8` and absolute tolerance `1e-11`. Other values and structure are checked for an exact match. The comparisons use three kinds of evidence.

| Panel | Fields | Comparison and retained result |
| --- | ---: | --- |
| A | 240 | Fresh outputs of the hash-pinned original script and aggregate inputs. Checks include 126-island HC3 fits, 123 resource-active and 91 port-positive sensitivities, and nine leave-one-Level1-out reruns. The introduced coefficient is +0.2854, native +0.0954, and the paired difference +0.1900 (HC3 95% interval −0.0187 to +0.3987). |
| B | 7 | Expected scope counts, the undefined pair’s region codes and both medians derived from the pinned source counts. Among the 7,875 pairs from 126 islands, 7,875 native and 7,874 introduced Jaccards are defined. The introduced BOU–SSA pair has zero union and remains undefined. Median Jaccard is 0.00448654 native and 0.04603046 introduced. |
| C | 1,700 | Archived original result fields. These comprise 1,308 one-island-held-out and 392 strict both-islands-unseen fields, including fold-level losses, descriptive intervals and the strict shortest-water comparator. Original code was rerun against the archived blocks before testing the projected inputs against those same blocks. |

Panel C uses **104 islands and 5,355 pair inputs**. The 104 one-island-held-out folds gave relative MSE reductions of **−0.110% native** and **+4.239% introduced**. The five strict folds evaluate 1,089 pairs with both islands unseen and yielded **−0.044% native** and **+2.735% introduced**. Descriptive fold-resampling intervals for the strict native (−0.984% to +0.431%), introduced (−1.649% to +5.508%) and paired introduced-minus-native (−1.693 to +5.400 percentage points) results all include zero. The electrical-resistance exposure is an undirected geometric measure. Figure 2D is a separate 126-island, fixed-239-butterfly potential-host analysis and retains the nine-Level1 paired slope contrast `p=0.0677`.

## Composite contents and evidence limits

A bytewise member comparison with the prior 28-file composite found 25 unchanged component or receipt members, 11 added Figure 2A–C members and three updated wrappers (README, manifest and runner). The new package has 31 component files, five provenance receipts and three wrappers. CRC, file names, source byte equality and direct identity-token screening passed on the actual ZIP. The projected tables contain botanical region codes, aggregate counts and derived pair features. Upstream archives, original occurrence records, WPI port coordinates and raw AIS inputs are excluded.

The Figure 1 and Figure 3B–C component replays aggregate results for the fixed 239-butterfly cohort. The appended Figure 3A receipt records a separate HOSTS/WCVP upstream rebuild. The Figure 2D component reruns the fixed-239 isolation and MEM10 calculations from processed inputs. The Figure 4 component recomputes the original 155 and 224 saved-score summaries. The distinct original 224 source aggregation and 676 model-fold refit have their own receipts. The included `source_qa_receipts/GBIF155_original_source_and_missing_RF_panel_gate.json` is a historical pre-recovery record. Original 155 recovery and model refitting are tracked by the separate private-source workflow. The 39-file archive contains neither the recovered 155 RF training panel nor a 155 model refit. [Original 155 recovery workflow](../../.github/workflows/audit-GEB-v05-original155-recovered-archive-20261122.yml)

The Figure 2A–C manifest carries an appendix-by-appendix S1–S8 coverage ledger. Within S3, the WCVP/WPI aggregate fits and circuit models have replay coverage. EuDiS, Tonga, botanical-effort adjustment, alternative metrics and Figure S3.1 require separate evidence. The replay leaves original WCVP incidence construction, the original WPI GeoJSON edition and conversion chain, and the upstream GSHHG circuit solver outside its scope.

**Release status remains internal.** Full SI replay, third-party redistribution clearance, indirect anonymity verification, stable reviewer access and a permanent archive DOI remain open. Direct token screening is a limited check. This audit records successful replay of the packaged computations. [Machine-readable verification receipt](../../results/geb_v05_private_source_recovery_20261120/FOUR_COMPONENT_INTERNAL_REPLAY_VERIFIED_20261122.json)

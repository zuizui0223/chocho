# Final scientific integration review

## Disposition and scope

**Ready within the reviewed scope, with the limitations already stated in the report. No blocking scientific wording, denominator or interpretation defect found.** This was a brief read-only integration review. No numerical suite, sampling run or source query was repeated. Only this review file was written.

Reviewed the integrated result and both new data-folder READMEs against the completed diagnostic/summary files and the retained independent source, taxonomy, repair/join and numerical reviews. Repository: `/workspace/scratch/10d2ad986b7b/chocho`.

| Reviewed prose | SHA256 at review |
| --- | --- |
| `provenance/archive/host_hindcast/ISLAND_DATED_HOST_LATITUDE_RESULTS_20261030.md` | `e1f2e8931eb94493a0dac87ea0ff59fee9bf7ec4950185fc27124c381852c0d1` |
| `data/exploratory/island_dated_host_20261030/README.md` | `42bb53869f3777ddaa529cae71442aa60e936be42fc92df4a1a1eaf3e6065dd3` |
| `data/exploratory/island_latitude_null_20261030/README.md` | `33669889052c86fb22b4283e13fd8ae70547ac58620adc6b3f04b6265d20ab70` |

The first prose pass used report hash `bf7361016bc922dd6a4bc57d513036c3c7cb676c7ecd2795946f1f6c8e5d9011`. After the descriptive host conditional-support subsection was added, its final text and additional group-variation counts were read and checked, yielding the reviewed report hash above. The worker's completed `/tmp/chocho-dated-workflow/host_conditional_support_review.md` was also read. Its copy at the report's relative data-folder link was installed during final assembly and its presence was verified. Final delivery-manifest assembly, any subsequent prose changes and publication remain root's integration responsibility.

## Verified scientific statements

| Topic | Verification and assessment |
| --- | --- |
| Distance endpoint and scale | Correctly defined as mean `log(1 + upper coast gap km)` over later-positive cells. The displayed ratio is `exp(observed T - exact conditional mean T)` on the geometric scale of **1 + gap**, and is explicitly distinguished from arithmetic mean kilometres. |
| Fixed denominators | The 2,385/66, 2,197/62 and 2,385/35 scenario denominators match the preserved diagnostic. Forced positives remain in the full positive denominators. |
| Ratios | Longitude-only 0.488059/0.549816/0.587710 and latitude-conditioned 0.777125/0.807778/0.712248 match the reviewed outputs. Changes are not described as causal percentages explained. |
| Tail probability interpretation | Hits 33/416/24, each from exactly 100,000 draws, and simultaneous upper bounds 0.000475536/0.004615193/0.000367814 match. These are conditional table-probability bounds, not ecological effect confidence intervals. Independent draws refer to tables, not ecological observations. |
| Conditional support | Variable positives 36/34/18, forced positives 30/28/17, variable species 26/24/15 and variable islands 11/11/8 match. Loss of randomization support is not portrayed as added replication. |
| Signed concentration | Primary 95.8902% and sensitivities 106.7219%/99.0569% are correctly described as signed arithmetic shares of the net observed-minus-reference statistic. The explanation for the share above 100% is correct. No causal or probability interpretation is attached. |
| Taxonomy coverage | 93 butterflies, 1,440 links, 1,062 source plant IDs; 1,040 successful primary queries plus 29 successful accepted-name checks; 22 unqueried source names. Final 695 concepts cover 865 links; 345 queried IDs plus 22 unqueried IDs and 575 links remain unresolved. All agree with completed outputs. |
| Butterfly catalogue completeness | The host summary confirms 84/93 butterflies have at least one resolved concept and only 4/93 have all source IDs resolved. The other nine are not recoded as having no hosts. |
| Recorded-host denominators | All 2,385 candidate queries are complete. The counts 336/2,385 = 14.09%, 20/66 = 30.30%, 1,961 no-match/unresolved, 88 no-match/fully resolved and 0 missing sum and agree with the summary and CSVs. |
| Descriptive fractions | The report explicitly notes that the 20/66 and 336/2,385 fractions share observations and compositions differ. Neither their contrast nor the 46 later-positive no-match cells is treated as an adjusted effect or biological host absence. |
| Plant evidence and repair | The resolved manifest confirms 28 complete targets, 0 missing, 0 complete-zero targets, 728,040 early plant records and 16 botanical-support descriptors. The original 25 complete/3 missing first-run manifest is retained and clearly distinguished from the authoritative resolved input. |
| Date and location limits | Correctly distinguishes 2010–2017 query membership from species-specific dates, same-event co-observation and the earlier +/-3-year point rule. Exact-polygon query membership is not presented as complete coordinate-uncertainty containment. |
| Biological interpretation | Catalogue compatibility plus a dated plant record is not described as feeding, larval development, resource abundance, origin status, wild rather than cultivated availability, or a causal host mechanism. |
| Port status | Remains unestimated. The original 23-island sample/support failures, 8/1/10/4 comparison cells, stopped 120-island baseline route and inability of fixed island margins to identify additive island-constant port effects are correctly retained. WPI proximity is not historical shipping. |
| Study status | The text discloses that butterfly outcomes and prior concentration were already known; this is exploration rather than preregistration, independent replication or a validated global mechanism. |

The READMEs point readers to the correct resolved plant manifest, retain the distinction between query SUCCESS and taxonomy eligibility, state the endpoint and date limitations, and identify the preserved first-run numerical outputs. They are consistent with the integrated result.

## Added host conditional-support diagnostic: independently verified

Joined `dated_host_primary_risk_cells.csv` to the primary rows of `latitude_exact_cell_occupancy_counts.csv.gz` by current butterfly key and target island. Both sides have exactly 2,385 unique matching keys, and every later-positive label agrees. This is a direct table reconciliation, not another null run or effect test.

| Conditional status | Recorded known-host concept | No match with unresolved catalogue | No match with fully resolved catalogue | Total |
| --- | ---: | ---: | ---: | ---: |
| Variable allowed cells | 17 | 72 | 0 | 89 |
| Variable observed positives | 11 | 25 | 0 | 36 |
| Forced positive cells | 9 | 21 | 0 | 30 |
| Never-positive allowed cells | 310 | 1,868 | 88 | 2,266 |

Consequently all 88 fully resolved no-match cells are never positive **under the current conditional fiber**. This is a constraint induced by the preserved margins, not evidence that butterflies cannot occur without hosts. The 17 recorded-host versus 72 unresolved no-match variable cells still provide variation in a recording proxy; they do not provide a resolved biological host-present versus host-absent comparison. It would therefore be inaccurate to say that every statistical recorded-host association is mathematically unidentifiable. The defensible next-step conclusion is that biological host-effect interpretation needs better taxonomic and local-recording support, with a fixed comparison and explicit treatment of recording intensity.

The four fully resolved source catalogues are Biblis hyperia, Hamadryas februa, Junonia genoveva and Lasiommata megera. Even their no-match cells are not ecological absence observations. These details support narrowing the next task toward unresolved concept interpretation and local evidence; they do not justify changing the current risk set or selecting a favorable null.

The added paragraph's grouping counts were independently recomputed from the same 89 variable cells: both recorded-host indicator values occur in **12 of 32 refined rows**, **10 of 26 species** and **7 of 11 target islands**. The text correctly states that such mixed labels do not alone prove that a complete-table score can vary under all joint margins. The added subsection neither computes nor implies a new host effect, null tail or biological absence comparison. Its next-step wording is supported.

## Integration recommendation

The current scientific account is appropriately bounded: latitude conditioning weakens the recorded-distance contrast, which remains concentrated in the western tropical group; actual early catalogue-compatible plant recording is now available, but the host mechanism and the port effect remain unestimated. Keep the explicit distinction between successful acquisition, complete taxonomy and recorded plant presence alongside the headline fractions.

No additional analysis or source retrieval is needed to deliver the reviewed findings. The added support subsection preserves the conditional meaning of never-positive cells and the distinction between recorded-host variation and biological absence. No required scientific correction is outstanding in the reviewed prose. The supporting review is present at its relative link.

## Final freeze addendum — 2026-10-30

Read the now-present host conditional-support subsection and its final wording, including the precise denominator **"32 species-by-longitude-by-zone rows containing variable cells"**. The current integrated report SHA256 is **`daa8b97f0ad37b9f804dcaee098e0808363ff96f922c6984f49e6ca342f1b489`**. Earlier hashes above retain the review history.

The independently confirmed mixed recorded-host indicator counts are **12/32 refined rows, 10/26 species and 7/11 target islands**, all scoped to cells with variable exact occupancy. The prose correctly avoids treating this as proof of a nonconstant complete-table host score or a biological host-present/host-absent contrast. **No scientific blocker remains.**

The supporting review is present at `data/exploratory/island_dated_host_20261030/host_conditional_support_review.md` and is byte-identical to the reviewed worker report, SHA256 **`75af7e6e56f11d49ae0ca8f5c1ef95db9b5a2d794c6b39835802eb8e72edd0eb`**.

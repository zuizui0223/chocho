# Original individual monarch larval development: host-resource conversion source gate

**2026-10-26, posthoc original-data audit.** Eight survival rows had already been inspected, so this is NOT a preregistered hypothesis or an independent discovery. Separate exploratory branch; frozen GEB main and PR #38/#70 untouched.

## Goal and data
Global butterfly host geography measures opportunity. The prior broad caterpillar field survey provided only 15 exact butterfly larval photos from one species, so could not identify butterfly recruitment. Fisher et al. (2026) Insects 17:215, DOI 10.3390/insects17020215, openly share butterfly **Danaus plexippus** laboratory experiments using upper/younger (Top) vs bottom/older (Bottom) Asclepias syriaca leaves. Original GitHub repo kelseyefisher/MilkweedLeafAgeMonarchSuccess at pinned tree acbb93c75e1f4810355c6ce03e6076e2303b2a98. It is not a nonnative plant-invasion experiment. Published authors have already studied these leaf-age effects. Our job is to **audit actual life-stage conversion, replication, and inference boundaries**, not to claim their treatment effect as new.

## Fixed posthoc analysis rules, before computing new trial-stratified summaries
1. Input exact SurviveToPupation.csv source Git blob 33f269e78f235cd1039394162a6871370db5e020; four source YearTrial blocks, two leaf treatments, exactly 60 assigned caterpillars per block-treatment. Outcome success = survived to **pupation**, not adult eclosion. Verify all eight rows, year-trial labels and aggregate source denominators.
2. Compute **common Mantel–Haenszel Top vs Bottom odds ratio** stratified by 4 YearTrial blocks. Report pooled top and bottom success proportion, risk difference, four separate block effects, and leave-one-block-out common-OR range. Under each block's fixed margins, generate exact hypergeometric distribution and convolution across blocks of total top survivors, with two-sided probability-order tail; this is an exact conditional test under **equal treatment effect null**, not robust generalization beyond the four experiment trials.
3. 10,000 bootstrap draws of **four whole YearTrial blocks sampled with replacement**, seed 20261026; show 95% percentile of common OR and paired mean risk difference, and proportion of resamples with Top better. Do not treat 480 larvae as 480 independent source environments.
4. Read day-6 larval mass (2324_GD_6DayMass.csv) and pupal development days (2324_GD_PupalDays.csv) using source-locked Git blobs, same YearTrial × Top/Bottom cells. Report available sample counts and median/mean as **conditional on measurement** and on reaching pupation. No imputation; do not infer successful metamorphosis, adult abundance, or individual stage transition from stage-specific nonmatched IDs.
5. Source leaf age/phenotype indicators (LeafCharacteristics.csv and PercentNitrogen.csv) describe manipulated leaves, NOT experimentally measured causal mediators for particular surviving individual butterflies. Do not conflate leaf characteristics with survival at same plant.
6. If total source counts, exact Git SHA, treatment names or sampling blocks differ, emit NOT EVALUABLE and archive source failure, never repair outcome by changing threshold or group identities. No additional outcome-directed hypothesis selection.

## Decision
This independent original-source experimental audit may demonstrate that **host leaf quality and juvenile developmental performance** are separate biological filters between mapped candidate host resources and realized butterfly reproduction. It cannot prove host relocation expands butterfly abundance, macro ecological competition, or fitness superiority of invasive milkweed. Published experiment remains original authors' result. Preserve a separately reproducible code and aggregate receipt; no raw individual animal records uploaded to chocho.


## Source-level response discordance: fixed two-scenario audit (after CI failure, before re-estimation)

The original source SHA-pinned trial summary SurviveToPupation.csv for block 20241 lists Top=30/60 and Bottom=36/60, whereas BOTH independent 262-row PupalDays and PupalMass files list Top=36 and Bottom=30 pupal outcome records for exactly this block. The other 3 blocks match perfectly, and original total pupae remain 262 regardless. This is a SOURCE INTEGRITY CONTRADICTION affecting treatment assignment, not evidence of biological inconsistency. The first CI correctly failed closed rather than silently repairing data.

Do NOT choose the truth by effect direction. Calculate and archive BOTH source-trust scenarios:

- S0 table-first: preserve the original 8-row pupation summary as is; original 20241 Top30, Bottom36.
- S1 pupal-record-first: change ONLY 20241 treatment-specific successful pupae to 36 Top, 30 Bottom as indicated by both pupal individual files, with the same 60 assigned in each. Explicit hypothetical sensitivity, NOT a correction to the published data.
- For both, compute the 4-block MH common OR, 10,000 four-trial cluster bootstrap, fixed-margin hypergeometric conditional probability-order tail, and leave-one-trial-out OR range.
- Individual growth/day data are reported separately as recorded in the original source with 20241 treatment assignment ambiguity flagged. No claimed survival-to-pupa individual linkage for 20241 is permitted. No renaming of source Top/Bottom.
- If interpretation changes across scenarios, label SOURCE_TREATMENT_LABEL_INDETERMINATE. Even if both agree, original source discrepancy remains unresolved. Do not infer adult emergence in either scenario.

This is a data-provenance audit after inspecting raw record disagreements, not an outcome-driven new ecological hypothesis.

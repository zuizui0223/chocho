# Nature Communications target: temporal resource × access evidence gate
**Version:** 2026-11-08 v0.1 | **Status:** frozen study-design and source-support screen; no new ecological model fitted

## Central question
When does the human-mediated redistribution of candidate larval-host plants become associated with subsequent butterfly community redistribution across sea barriers, beyond geography, climate, and recording effort?

This is an **observational, testable question**. Current archived plant–butterfly relations are not equivalent to established butterfly breeding or verified plant transport; do not claim otherwise. The target journal is *Nature Communications*, not a statement of journal acceptance/readiness.

### Distinct contribution against contemporary work
Gross, Kawahara & Daru (2026), `Climate and regional plant richness drive diet specialization in butterfly caterpillars` (*Nature Communications*, DOI 10.1038/s41467-026-73236-4) addresses climate/plant-family richness and diet specialization. Our proposed novelty is **historically ordered redistribution of potential resources interacting with access to predict new island butterfly indexing**. It must be supported by truly time-prior host identities, independent island/species/geographic transfer, and evidence insensitive to observation bias—not by re-demonstrating an existing static resource association.

## Fixed existing source populations (NEVER pool denominators)
- **Resource geography (separate existing JBI/GEB input):** 239 butterflies and 126 WGSRPD3 insular / 237 mainland botanical regions; +118.9% vs +40.9% potential resource expansion, matched +10.7 pp share; this is a **contextual baseline, not proof of colonisation**.
- **Separate existing butterfly prediction validation (context, not the intended mechanism test):** 155 butterflies × 208,482 held-out species×cell×period rows. Archived temporal and spatial species-equal MSE improvements from adding static host identities are +4.1202% and +2.9952%. These results predict standardized record probabilities and do not document host arrival dates or independent butterfly establishment. Do not merge these denominators with the 225- or 239-butterfly source populations.
- **Island/port assemblage source:** a different **225-butterfly candidate panel**, 51 original landmass groups, 194 original botanical port buffers, 182 marine-qualified electrical terminals, and 1,275 non-independent island pairs. This is the source for a prospective *record-index* risk-set design, not 239.
- **Archived 2001–2010 prior host identification:** one pinned per-butterfly `island_dated_host_match_all225.csv` and independently verified `island_dated_host_match_validation.json`. The common newly indexed 2011–2022 sample has 82 butterflies, **1,163 new butterfly×island record pairs**, only **5** with an explicitly matched *recorded* prior host (five distinct butterfly species). This is measured prior-record overlap, not known resource absence elsewhere.
- **Archived 2001–2010 general botanical indexing:** 194/194 original ports re-queried and audited; 149 indexed ports and 45/51 islands with at least one historical GBIF Plantae speciesKey. General Plantae speciesKeys are NOT verified larval hosts.
- **Modern marine:** uniform geography-only circuit plus AIS density (2015–2021) and mean wind speed. Modern AIS and wind are **not historically pre-2011 predictors**. Uniform geometry can be used as an a priori access proxy; historical port connectivity (2006–2010 UNCTAD PLSCIQ index) is available only on 17/51 source groups, and is an index, not voyages.

## Publication-defining prospective test (not yet run)
**Unit:** butterfly species `b` × island `i` × validation window `t`, based on a fully enumerated source-indexed risk set. Never use only positive butterfly records as a full risk set. Distinguish: not indexed; indexed before baseline; newly indexed later; failed/no-effort source. Do not interpret non-indexed as ecological absence.

**Exposure:** `PriorHost_bi(<=2010)`, a species-specific count/indicator of named/accepted larval host plant taxa with independently indexed local evidence from 2001–2010 or earlier, *not* the 2001–2024 pooled plant inventory. The original static accepted host relations can define identity but cannot document local use or date of introduction. Keep a binary indicator of source effort and disaggregated coverage; distinguish a true record 0 from missing indexing.

**Access proxy:** source-defined marine effective resistance from uniform ocean field. For potential source islands use only **pre-2011 butterfly-indexed** source observations in the risk-set reconstruction. Distances, island area, climate, and pre-2011 botanical/butterfly source effort are controls; coastline resistance is not a directed path or documented stepping-stone relay. Modern AIS is sensitivity/context only. If climate for an island and year is unavailable, do not claim climate-controlled estimates.

**Primary hypothesis:** the *interaction* between measured prior host suitability/availability and source-indexed access improves out-of-time/both-islands-unseen prediction of later butterfly indexing beyond additively including both, geography and effort. A generic potential-resource increase is insufficient. Report resource-only, access-only, additive, interaction models on an identical held-out sample and loss metric.

**Competing explanations / falsifiers:** all-plant richness and observation effort; non-host plant identities with frequency-matched sampling; placebo host–butterfly links preserving degree and plant family; geographical distance; climate; local AIS graph-node label placebo; independent taxonomic/spatial validation; earlier butterfly community baseline. No post hoc model selection based on a desired positive sign.

**Evaluation:** record-index endpoints not colonisation; calibration/log-loss/Brier/PR-AUC where class support permits, and island+species-blocked evaluation with island family-wise resampling. Lock all primary metrics and folds before model outputs. Report coverage/selection across *all* 225 species and 51 islands, and any species overlap between test and training.

## Present hard gate and interpretation
The archived early-host **case-only** corroboration is **5 of 1,163** later butterfly record pairs (0.43%). This alone does **NOT** establish that prior hosts are absent from the remaining 1,158; it does show that a model claiming an observed *species-specific pre-2011 host trigger among subsequent positives* lacks observed support with current sources. The earlier family-preserving null (44/199 draws at least as large for pair weighting) does not establish exceptional prior-host matching. It is not a valid independent p-value.

**Gate:** Do not fit an interaction and call the result Nature Communications-ready using modern 2001–2024 candidate hosts or the 194-port pre-2011 *general* flora as though those were identified early hosts. First recover the actual time-indexed accepted host-plant identities in a reproducible public source and the negative/effort-defined species×island at-risk population; quantify source coverage per island/species, not just global counts. A conservative **exploratory feasibility benchmark** is at least 50 independently source-supported later *case* matches across at least 10 butterfly species and 10 island groups; meeting this threshold is only a screening rule, **not** a power calculation or sufficient proof. If it cannot be met, analyse only record-based correspondence and do not claim real distribution response.

## Nature Communications editorial and replication constraints
- Article's significance must advance understanding rather than repackage descriptive results. Use a clearly new *dynamic/interaction* finding only if independent temporal tests actually support it.
- The journal does not accept pre-submission enquiries and expects the full manuscript; no journal submission is authorised here.
- The minimum data used to interpret, replicate and extend conclusions must be available to readers/reviewers under disclosed terms; source-locked output summaries do NOT replace unavailable original HOSTS×port identity rows. Scientific manuscripts that significantly overlap other under-review papers are not independently publishable without review/disclosure.
- Preserve the original GEB global +54.9% paper and the existing JBI Draft PR #151 without modifying either, pending explicit manuscript-overlap strategy with coauthors.

## Source receipts and literature
- `provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_all225.csv`
- `provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_validation.json`
- `provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_host_family_null_dated_audit.json`
- `provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_audit.json`
- `results/island51_pre2011_public_general_flora_20261107/source/GBIF_original51_early_general_flora_islands.csv`
- `results/island51_pre2011_public_general_flora_20261107/source/GBIF_original194_port_2001_2010_coverage.csv`
- [Nature Communications aims](https://www.nature.com/ncomms/ncomms/aims), [peer-review criteria](https://www.nature.com/ncomms/editorial-policies/peer-review), [Guide to authors](https://www.nature.com/ncomms/submit/guide-to-authors), [data/code availability](https://www.nature.com/ncomms/editorial-policies/reporting-standards), [Gross et al. 2026](https://doi.org/10.1038/s41467-026-73236-4)

**Inference freeze:** Existing positive geography and port-flora associations remain unchanged. No ecological claim of ship-mediated plant transport, successful larval feeding, reproduction, source island travel or colonisation has been earned by this gate.

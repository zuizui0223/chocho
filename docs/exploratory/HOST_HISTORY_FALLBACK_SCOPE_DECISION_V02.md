# Ecological hypothesis narrowing after the failed sequential-choice pilot

**Date:** 2026-10-10. **Status:** independent prior-art and source-feasibility review; **not a new biological finding**. GEB submission PR #38 unchanged.

## What was newly checked in this continuation

1. Verified previous *Euphydryas editha* serial host-choice source: 145 trials, 29 maternal identities, 93 sequential eligible trials from 27 females. Added lagged *Plantago* egg share to date/trial-order model: held-out log-loss improvement **+0.00448** with maternal-cluster 95% interval **−0.01817 to +0.02566**; the predeclared predictive gate failed. Do **not** search alternate lags or categories after this result. Original source / receipt: https://github.com/zuizui0223/chocho/actions/runs/38012350290.
2. Revisited Jones & Agrawal (2019) *Oikos*, https://doi.org/10.1111/oik.06001, Dryad https://doi.org/10.5061/dryad.8hd6764. They already **randomized female oviposition experience** and showed a preference shift, also measured larval performance separately. Choice hosts were two **native** *Asclepias incarnata* subspecies, not an experimental nonnative-host introduction. Their paper itself reports approximately 22% pupation on initially preferred *A. i. pulchra* versus 43% on alternative *A. i. incarnata*, with no significant four-host survival difference (Fisher p≈0.17); these are **published group summaries, not reanalysed data or adult recruitment**.
3. Original Dryad `oviposition_experience.xlsx` and `caterpillar_performance.xlsx` were both attempted through fixed individual-file URLs on GitHub Actions. **Both returned HTTP 403**; the 3 source-integrity/unit tests passed, but `source_bytes_complete=false`, `row_level_schema_audited=false`, `effect_estimated=false`. Audit run: https://github.com/zuizui0223/chocho/actions/runs/38012834051. Do not cite CI success as biological confirmation.
4. Gowri & Monteiro (2024), *Annals of the New York Academy of Sciences*, https://doi.org/10.1111/nyas.15090, already experimentally tracked an **acquired food-odor preference** for five generations and its disappearance after withdrawal. Learned preference, multigenerational persistence and withdrawal per se are therefore not undiscovered biology. Importantly, this manipulated an odor added to leaves, **not** replacement of real native/exotic host-plant species.
5. Singer & Parmesan (2019), *Evolutionary Applications*, https://doi.org/10.1111/eva.12775, already documented host-adoption reversals and experimentally measured life-history consequences in *E. editha*. Host reversion itself should not be a headline novelty. Exposing a hysteresis pattern in a single separate population is not sufficient to supersede this long line of work.

## Decision: no reanalysis-only second ecology paper from these sources

**Do not** promote the E. editha exploratory result, Jones & Agrawal's published learning experiment, isolated native/exotic performance contrasts or static chocho host-envelope numbers to independent confirmation of ecological hysteresis. Current resources lack a single combined cohort where randomized **host-use history**, standardized contemporary **host availability and removal**, and **viable offspring production** all occur together.

A defensible new empirical target is:

> Does experimentally imposed past reliance on a redistributed host change offspring production after that host disappears, **when present edible native-host biomass is equal**? Is any loss of recruitment due instead simply to the decrease in remaining biomass?

The key is not just whether host choice remembers the past. It is whether past host exposure changes a butterfly's **functional ability to recover** and whether that effect can be distinguished from contemporary resource scarcity.

## Revised factorial experiment — three current-resource arms

Randomize independent maternal lineages (nested within multiple populations and at least two feasible butterfly species) to two past-resource regimes over experimental generations, then randomly assign these lineages to a current resource condition:

| Historical experimental supply | Present mixed hosts, sum biomass B | Present native-only, compensated to B | Present native-only, uncompensated B/2 |
| --- | --- | --- | --- |
| native-only experience | H0 × R | H0 × C | H0 × S |
| native+introduced experience | H1 × R | H1 × C | H1 × S |

- **Primary estimand:** `(H1_C - H1_R) - (H0_C - H0_R)`. This tests a *history × host-composition* interaction **at equal total edible biomass**, with randomized assignments to history and current treatments.
- **Secondary estimand:** `(H1_S - H1_C) - (H0_S - H0_C)`. This isolates the added exposure to *resource shortage* beyond host replacement at matched final identity.
- **Primary measured outcome:** **flight-capable adult offspring produced per initially assigned mated female**, retaining all zeroes and failed lineages; not egg count, GBIF observation, coarse WGSRPD3 opportunity, or inferred carrying capacity.
- **Mechanistic secondaries:** female acceptance of ancestral host, eggs on native hosts, tracked neonate-to-adult survival, plant stage and usable biomass, host chemistry, enemy exposure. Behavioral choice alone does not imply a fitness difference.
- **Validity conditions:** admitted host species must support development to flight-capable adults; matched plant developmental stages and measured *available edible biomass through time*; individual family/lineage is replication unit; maintain multiple populations and population-level variation; report treatment attrition and zero reproduction as randomized.
- **Interpretation:** a nonzero primary interaction can identify a causal effect of *experimentally assigned recent history under controlled laboratory conditions*. It is **not** equivalent to a causal effect of centuries-long botanical invasion or genetic evolution. Deliberate multigeneration/common-garden contrasts and replicated host histories are required to attribute mechanism.

**Positive, null and opposite predictions all matter.** Positive history-dependent loss supports impaired fallback, null supports functional reversibility under matched supplies, and the opposite sign supports enhanced fallback after broad-host history. If effect only appears in the uncompensated arm, then short-term *resource quantity* rather than host-use memory may explain vulnerability.

## Scope and administrative boundaries

- Frozen design: `docs/exploratory/RANDOMIZED_HOST_HISTORY_FALLBACK_PROTOCOL_V01.json` (v02 revision before any new experimental outcome).
- There is **no** newly estimated biological effect for this six-cell factorial design, and no claim that experiments have been performed.
- Existing PR #38 remains the ready-to-submit descriptive **potential host-resource geography** paper, not altered by this exploration.
- Stop acquiring more static butterfly/host geographic Jaccards to manufacture significance; prioritize an experiment or an independently collected dataset containing every necessary exposure and offspring-fate field.

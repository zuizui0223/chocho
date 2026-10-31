# Submission decision — GEB integrated versus split-paper tracks

**Date:** 2026-10-31  
**Repository:** `zuizui0223/chocho`  
**Decision:** **Draft an integrated GEB Research Article as the first-choice alternative, while preserving the already frozen GEB v0.2 manuscript as a fully viable independent submission.** Do not submit or overwrite either version without final coauthor approval.

## Why GEB first

The primary ecological hypothesis is *partner redistribution can change consumer resource geography without any shift in recorded dietary breadth, and that geography contains predictive information about where the consumer is recorded*. This is a general, macroecological question. The measured outcomes are global butterfly × WGSRPD3 resource opportunity (239 species), held-out GBIF recording probability (155 species), and exploratory island-pair assemblage correspondence (51 islands). GEB specifically seeks broad-scale ecological patterns and globally relevant conceptual advances, including anthropogenic impacts. This is a tighter fit than using 51 selected islands to headline an invasion-transport mechanism that the data did not identify.

GEB official scope: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/ProductInformation.html  
GEB author guide (checked October 31, 2026): https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html  
**Category:** Research Article. Typical main text ≈5,000 words, typically 6–8 display pieces, structured abstract ≤300 words, double-anonymous review, free-format initial submission. Its official submission portal is https://authors.wiley.com/journal/GEB.

## Journal ranking conditional on scientific framing

| Priority | Journal | When it is suitable | Main editorial/reviewer risk |
|---|---|---|---|
| 1 | **Global Ecology and Biogeography** | General partner-redistribution biogeography; 239-species resource reconstruction PLUS out-of-time/spatial host predictive validation | Three sample frames must read as one question; an uncontrolled harbor correlation must not be mistaken for confirmation of colonization |
| 2 | **Ecography** | When an integrated manuscript foregrounds time-/space-heldout species-distribution predictions, spatial modelling and the host-aware information gain | Need a clearly transferable methodological/ecological contribution, not just an after-the-fact boost in GBIF index prediction |
| 3 | **Journal of Biogeography** | When stronger independent island or historical biogeographic contrasts become primary, or GEB editors seek a more geographical framing | Pooled post-outcome port plants, nonindependent island pairs, old assemblage correlation exceeding new, and geography-preserving nulls are serious weaknesses for a gateway mechanism |
| Conditional | **Global Change Biology** | If an independent biological response of butterfly ranges to global environmental/plant change is documented with stronger temporal support | Current gains are *resource potential and recording prediction*, not an identified biological effect of host introductions on colonization |
| Conditional | **Diversity and Distributions** | If conservation implications are independently substantiated and central to the paper | At present there is no measured restoration benefit or plant-removal risk; added conservation relevance would be speculative |

Verified journal scope:  
GEB: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/ProductInformation.html  
Ecography: https://onlinelibrary.wiley.com/journal/16000587  
JBI: https://onlinelibrary.wiley.com/page/journal/13652699/homepage/productinformation.html  
GCB: https://onlinelibrary.wiley.com/page/journal/13652486/homepage/productinformation.html  
D&D: https://onlinelibrary.wiley.com/page/journal/14724642/homepage/productinformation.html

## The single manuscript claim: four measurement levels, not four independent confirmations

1. **Opportunity:** fixed HOSTS identities + WCVP regional statuses; 239 butterflies, 206 expanded, 26,530 native butterfly × region units + 14,553 added = **+54.9%** (crop-host exclusion +46–47%). 38 host species contribute half of added units. Host-family breadth vs proportional gain: Spearman ρ = 0.008, CI crosses zero; **do not claim equivalence** (post-hoc TOST p = 0.078).
2. **Occurrence geography:** 32-species independent climate-stratified GBIF validation panel, **66/115** outside-native species × region observations inside introduced-host contemporary resource envelopes versus region-stratified matched null. Separate source and selected panel; neither observed feeding nor time-indexed host introduction.
3. **Predictive information:** 155 already-supported butterflies, **208,482 held-out evaluation pairs**; host-resource and host-family composition on top of background + climate + land cover lowers equal-species-weighted MSE **4.12%** temporal, **3.00%** spatial and **2.74%** at >100 km from fitting locations (post-hoc subset). Four sensitivity tests are positive; new-species transfer is *not* improved. These are *retrospective GBIF recording-probability evaluations* using static host information, NOT a prospective dispersal/colonization forecast.
4. **Island-port spatial correspondence:** 51 islands, 1,275 dependent pairs: pooled port candidate introduced-host Jaccard vs new butterfly Jaccard **ρ = +0.300**, distance + prior-butterfly partial rank **+0.215**; prior butterfly assemblage correlation is higher (**+0.418**). A 225-species comparison gives GEB gain vs island-port candidate-host count **ρ = +0.652**, but both share host catalogue inputs. In an original 17-island pre-2011 plant series, geographic- and plant-prevalence-matched nulls **do not support an independent host-index effect** on later first-index butterfly records. Keep pooled port data as spatial *association*, not temporal mediation or ship route.

**Main abstract and conclusion:** levels 1–3 provide the general GEB contribution; level 4 is exploratory evidence of the geographic interface with islands, and a conspicuous *negative control* against causal gateways. The island signal alone must not be described as validation of +54.9% realized resource opportunity.

## Recommended integrated article layout (approximately five main figures)

- **Figure 1 — Where opportunities change:** a conceptual/map panel plus 239-species native and contemporary resource-unit distributions (26,530 → 41,083).
- **Figure 2 — Which hosts make the gain:** 670 contributors, 38 hosts generating half, host-family breadth vs proportional gain (ρ = 0.008), crop exclusion as SI.
- **Figure 3 — The observed butterfly link:** 66/115 outside-native butterfly × region observation recovery compared with matched geometric/null expectation.
- **Figure 4 — Hosts improve held-out recording predictions (NEW)**: M3 vs M5 temporal and spatial MSE, ≥100-km sensitivity; 4 observation sensitivities in SI; explicitly flag unseen-species failure and recording-probability target.
- **Figure 5 — Island/port spatial correlation (NEW)**: 51-island pair Jaccard association with old-butterfly and distance-controlled comparisons, plus 17-island geography-matched historical null as a clear negative control. Use island jackknife for dependence and describe pooled 2001–2024 plant data.

Existing v0.2 already has **3 main figures and 1 supplement**. New Figures 4 and 5 must be generated from existing fixed source data and audited; merely referencing them is not enough. Figure numbering is provisional and existing figures may need consolidation.

## Submission gate (work still needed)

1. **Manuscript integration:** create a distinct complete alternative `integrated_GEB_research_article_draft_v0.1.md` with separate Methods, Results, and Discussion for prediction/islands. The existing `manuscript/butterfly_specialization_ecology_v0.2.md`, its figure manifests, frozen claim map and anonymous reviewer bundle remain byte-identical.
2. **Word and focus edit:** enforce abstract ≤300 English words and main body approximately 5,000 words. Remove duplicate climate/host-null details from main body to SI before any submission. Keep main ecological discovery clear in title, first abstract sentences, figures and cover letter.
3. **Science and source:** report sample size and observational outcome for each module; audit time leakage from contemporary WCVP and pooled 2001–2024 port plants. The 155-species forecast improvement must not be converted into a benefit of *introduced* hosts or a causal port effect; no independent species-generalization improvement.
4. **Author and figure checks:** produce/validate new Figure 4/5 PDF/vector + captions, add new module data-source and SI text, regenerate blind version, update methods/data availability, anonymous code/records and distinct versioned cover letter.
5. **Final submission approval:** corresponding author details, author affiliations/ORCIDs, author contributions, funding, declarations, data permissions, all-coauthor approval and portal submission are still author-only steps; **no journal submission has occurred**.

## Stop/split rule

**Choose integrated GEB** only if the prediction cohort and island exploratory association can be explained clearly inside ~5,000 words **without upgrading their causality** and if new figures/records can be reviewed. Otherwise **submit the current frozen GEB v0.2 intact** first and develop a second, genuinely distinct island-port paper only after the same-island historic port versus nonport host detection and independently dated alien-status/recording controls yield a novel result. Do not submit overlapping manuscripts claiming the same 54.9%/host-source dataset as two independent discoveries.

**Important provenance:** The older draft cover letter claims authors will approve the submitted version; that is a template, not evidence of actual approval. No author metadata or journal account access has been inferred.

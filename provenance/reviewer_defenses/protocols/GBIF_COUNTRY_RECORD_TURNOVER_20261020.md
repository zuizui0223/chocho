# Frozen follow-up: country-level butterfly detection turnover, 2010–2017 vs 2018–2025

**Protocol fixed before examining paired country-level transition outcomes.** This
is an exploratory source-conditional diagnostic of GBIF detection history, not
a preregistered test, a dated host-plant introduction experiment or a test of
actual butterfly colonization. Do not change the frozen GEB submission.

## Inputs and fixed target

Reuse the completed **924 country facet queries** from 20 October 2026 and the
same frozen 229 non-colliding butterfly names. The 8 unresolved names and
both names in the Oeneis accepted-key collision remain excluded. Use only
2010–2017 and 2018–2025; do not include the partial 2026 calendar year.
Read the frozen 255 species × botanical Level-1 resource-contrast strata,
including their country-by-species candidate cells, *without* refitting or
retuning plant ranges, geography, native hosts or climate filters. The
baseline candidate set comprises 3,188 species × country cells; only 48 were
positive in the pooled accepted-basis comparison.

A cell with no indexed record is called **not recorded**, never absent.
An indexed record is not verified wild presence, breeding or establishment.
Record dates do not tell when the host plant arrived. No host-use inference is
permitted from country overlap.

## Frozen reporting background and outcomes

For each basis scope separately, count records by species × country × window,
using the existing frozen CSV. The accepted scope combines the five accepted
basis-of-record categories; the core scope has HUMAN_OBSERVATION and
PRESERVED_SPECIMEN only. They overlap; never sum them.

For every focal species × country cell, other-species effort is the number of
**other** non-colliding, completed butterfly taxa with indexed records in that
country in **each** window. Exclude the focal species from both efforts, and
do not use its record to qualify itself. Retain only countries meeting the
chosen minimum *in both windows*, so newly recorded focal butterflies cannot
solely create their own effort reference. This measures taxon-recording
opportunity, not calibrated field survey effort.

Record the four binary states: 00 (neither window), 01 (2018–2025 only),
10 (2010–2017 only), 11 (both windows). The primary comparison is **01
among early-unrecorded 00+01 cells**; it is a contrast in GBIF discovery,
not real colonization. Hold the early-unrecorded pool fixed within each
species × Level-1 stratum; retain 0<A<N exposure contrasts. With N candidate
countries, A introduced-added-host countries, K 01 countries and O 01
added-host countries, use expected overlap AK/N and excess
100 × sum(O−AK/N)/sum(K) pp.

As prespecified descriptive falsification, separately evaluate **10 among
late-unrecorded 00+10 cells** using the same statistic, and report 00/01/10/11
state counts by exposure. Opposite temporal asymmetries are not a test of
extirpation or demographic loss.

## Fixed diagnostics and information gate

Primary: accepted bases, focal record threshold >=1 in each window, >=3
other-focal taxa recorded in **both** windows. Sensitivities:
- basis: accepted vs core;
- minimum other-taxon recording opportunity per window: 1, 3, 5;
- focal recorded threshold in a positive window: >=1 vs >=5.

These give 2 × 3 × 2 = **12** transparent specifications, applied to both
directions, with no after-the-fact search for a best result. Negative/zero
candidate cells do not prove species absences, and a transition based on >=5
means <5 counts as unrecorded for that diagnostic only.

Compute an exact conditional hypergeometric reference tail or validated
Monte Carlo approximation, genus-cluster bootstrap 95% intervals, and the
number of strata with 0<A<N and 0<K<N. As in the prior protocol, **50**
informative strata, **25** taxa and **25** positive target cells are the
minimum for an exploratory precision claim. Below these thresholds, mark the
inference unresolved regardless of point estimate. Never claim causality
from conditional-tail values, positive intervals or the GBIF source size.

## Audit and reporting

Save full source hashes, 12 complete combination tables, primary state table,
which countries and genera carry positive transitions, information gates, and
a repeatable code/test path. Report date-window differences in GBIF observation
effort using other-species richness; report geographic concentration and
sensitivity to Caribbean/Central America as posthoc checks only. In particular,
do not describe late-only 01 as expansion, lagged population recruitment,
establishment, successful larval development or altered butterfly diversity.
Keep these exploratory files out of the frozen GEB anonymous review snapshot.

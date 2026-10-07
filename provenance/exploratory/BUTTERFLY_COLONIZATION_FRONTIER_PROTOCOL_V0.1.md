# Prospective butterfly colonization-frontier validation v0.1

**Frozen:** 2026-10-07  
**Status:** prospective prediction set; do not update candidate membership using future occurrence outcomes.

## Purpose

This protocol freezes a falsifiable prediction set derived from the current resource-rewiring analysis before any future occurrence records are used for validation.

The ecological prediction is:

> Butterfly occurrence should be more likely to appear in introduced-only host-resource regions that are already climatically similar to held-out observed regions and have adequate background sampling effort.

This is a prediction of **future detection / realized colonization opportunity**, not guaranteed colonization.

## Frozen candidate definition

A species × WGSRPD3 region entered the candidate pool only if, at the time of freezing:

1. the region lay outside the species' native-host resource envelope;
2. the region entered the contemporary resource envelope through introduced host distributions;
3. the butterfly had not been observed in the current occurrence panel in that region;
4. background sampling effort from the independently assembled butterfly panel was sufficient;
5. climate mismatch was estimable from the frozen cross-fit climate model.

The full eligible pool contained **806 species × region units across 24 species**.

The high-confidence prospective set frozen here contains the **44 units with climate-compatibility rank >= 0.80**. The exact list is stored in:

`provenance/exploratory/butterfly_colonization_frontier_high_confidence_v0.1.csv`

## Primary future test

At a future prespecified census date, query new butterfly occurrences **after the freeze date** without changing the candidate set.

Primary comparison:

- detection incidence among the 44 high-confidence candidates;
- versus eligible frontier units with climate-compatibility rank < 0.50 from the same frozen 806-unit pool.

Species identity should be treated as a cluster in inference.

## Secondary tests

1. Time-to-first-detection as a function of frozen climate compatibility.
2. Effect of background sampling effort.
3. Whether candidates connected to existing resource geography by introduced-host stepping stones are detected earlier.
4. Whether high exact-host co-user exposure predicts slower realization, as expected if biotic competition offsets resource availability.

## Interpretation boundaries

A future new record is not proof that the introduced host caused colonization or that the butterfly locally uses that host.

Failure to detect a butterfly does not prove ecological absence.

Local host use, host quality, abundance, phenology, dispersal barriers, enemies and habitat availability remain unmeasured filters.

The prediction ledger must not be retrospectively edited after future occurrence outcomes are inspected.

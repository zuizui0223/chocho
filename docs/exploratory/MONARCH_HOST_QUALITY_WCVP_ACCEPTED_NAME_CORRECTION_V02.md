# WCVP accepted-name priority correction for monarch quality panel v0.2

**2026-10-08; corrective taxonomic audit, decided BEFORE seeing corrected outcomes.**

The initial 127-species WCVP matching run 37744576643 successfully extracted Greenstein et al. (2022) classifications, but treated 16 exact botanical strings as ambiguous because both an **accepted botanical name record** and one or more **synonym/homonym records** resolved to multiple accepted IDs. Thus the original Asclepias primary panel had only 44 accepted taxa and **wrongly excluded** `Asclepias curassavica`, a known independently high-performance H3 host with 41 native and 97 introduced WCVP regions. The original pilot's null p=0.5561 is retained for audit only and **MUST NOT be quoted as the final test of the planned taxonomic population**.

## Correct matching rule before outcome inspection

For each of the **originally frozen 127** scientific names:

1. Search WCVP records with exact case-sensitive `taxon_name` and `taxon_rank=species`.
2. If there is **exactly one unique** accepted record ID among rows satisfying `taxon_status=accepted` and `plant_name_id=accepted_plant_name_id`, select that accepted ID **even if a homonymous or synonym spelling also exists**. This is an exact-accepted-name preference, not genus matching or a search for high-performance taxa.
3. Otherwise accept a spelling as a synonym only if all matching species-name rows unambiguously resolve to one accepted plant ID.
4. If multiple accepted IDs survive after accepted-name priority, exclude and report as ambiguous. Never guess identity.
5. Preserve Greenstein 2022 high/low classes exactly, original WCVP distribution flags, native-range control, quartile permutation seed, sample construction and hypothesis. Rerun all predeclared primary/secondary/sensitivity tests; **do not choose among v0.1 and v0.2 based on p-value**.

Source frozen botanical reference: `matildabrown/rWCVPdata` at commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`. Expected positive control: accepted `Asclepias curassavica` ID `500848`, 41 native / 97 introduced regions. If the unambiguous accepted ID is not recovered, abort corrected analysis.

The corrected test remains **exploratory**, not an independent confirmatory analysis. Strongly evidence-ranked H3/L3 versus weakly evidenced H1/H2 etc. are not equivalent measurements of quantitative fitness; WCVP range is not local plant biomass or butterfly occupancy. Do not change GEB main manuscript.

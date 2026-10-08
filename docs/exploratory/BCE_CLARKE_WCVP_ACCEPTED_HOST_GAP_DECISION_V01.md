# Evidence-tier-filtered European butterfly host links after WCVP taxonomic reconciliation

**Date:** 2026-10-08. **Status:** ACCEPTED-ID CROSSWALK COMPLETE / POTENTIAL SOURCE-COVERAGE GAP / ECOLOGICAL AND GEOGRAPHIC INFERENCE STILL OPEN. No changes to original manuscript/PR #38.

## Explicit scope and provenance

The source is the [Butterfly Conservation Europe (BCE)](https://www.bc-europe.eu/taxonomy.php) species checklist **derived** from Clarke (2024), *Ecology and Evolution*, DOI [10.1002/ece3.10834](https://doi.org/10.1002/ece3.10834), which filters host records to source evidence ranks **1–3**. The original relational Dryad records [10.5061/dryad.1vhhmgr12](https://doi.org/10.5061/dryad.1vhhmgr12) could not be retrieved from this environment (401/403 and supplementary timeout); BCE suppresses individual link references/evidence grades and may share underlying references with original global HOSTS. **Independently curated does not mean original biological observations are independent.**

The full 92 exact butterfly binomial overlap with chocho's 239 species was source-verified in [run 37777242991](https://github.com/zuizui0223/chocho/actions/runs/37777242991), with 83 interpretable species-level BCE foodplant pages. BCE–HOSTS original lexical comparison on those 83: 1,390 BCE and 1,068 HOSTS associations; **338** plant scientific-binomial spellings match, **1,052** BCE exact spellings not present on the frozen HOSTS butterfly-host list, and 730 original HOSTS host spellings not present on BCE.

## WCVP accepted plant identity resolution

[Run 37777881774](https://github.com/zuizui0223/chocho/actions/runs/37777881774) successfully retrieved the fixed source `rWCVPdata` commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`, preserved frozen HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`, generated accepted-name sidecars, and reconciled all **1,052 BCE-only botanical binomial spellings** with unique accepted WCVP species identities. These 1,052 links involved **689 distinct BCE botanical binomial strings**: **685** uniquely resolved through exact accepted species status and **4** through synonym-to-accepted mapping (under the fixed protocol).

When checked against the original *same butterfly × accepted plant ID* frozen HOSTS pairs, the complete 1,052 source-level candidate pairs split into:

| Outcome after unique WCVP botanical taxon ID reconciliation | Butterfly–plant pairs |
| --- | ---: |
| Accepted ID already exists as a frozen butterfly–host edge (nomenclatural duplication) | **25** |
| Accepted ID not represented in frozen HOSTS for the same butterfly | **1,027** |
| Unmatched or ambiguous WCVP accepted species ID among these 1,052 | **0** |

**Interpretation:** Exact botanical synonym/name inconsistencies explain only 25/1,052 lexical differences. A large **source-coverage divergence remains after accepted-ID matching**. This is more serious than a simple spelling problem, and directly challenges any assumption that the fixed HOSTS inventory exhaustively represents documented larval foodplants for these European species. Nevertheless, the remaining 1,027 should be called **candidate additional published host edges**, NOT 1,027 independently witnessed feeding events or universally viable resources.

## What remains unresolved

1. **Source provenance and biological quality:** Although BCE says the list retains Clarke's evidence grades 1–3, its pages do not identify which grade or original paper supports each link. Grade 3 includes eggs/oviposition and does not establish development to adulthood. HOSTS and BCE may include overlapping primary citations.
2. **Geography:** Clarke/BCE provides host use in Europe; the fixed chocho model extrapolates every recorded host association to every listed WCVP WGSRPD3 plant region worldwide. Extending 1,027 European candidates in that same manner is deliberately a **one-direction upper-bound geographical sensitivity**, not actual global feeding. Conversely, excluding them does not demonstrate true biological absence.
3. **Butterfly taxonomy and reporting scale:** Only 83 species-level host pages of the 92 exact-European butterfly overlaps (out of chocho's 239 species) enter these counts, selected by the existing dataset's exact nomenclature. No inference to other geographies or butterflies is justified.
4. **Opposite differences:** 730 original HOSTS names lacked exact BCE counterparts. Neither checklist is a complete truth set, and different inclusion criteria mean adding the BCE-only ones is one-directional rather than full reference-set reconciliation.

## Prespecified follow-on sensitivity

The independent protocol `BCE_CLARKE_ONE_DIRECTION_HOST_AUGMENTATION_PROTOCOL_V01.json` was committed **before the WCVP accepted-ID results appeared**. It adds ALL 1,027 accepted-ID-absent candidate pairs for applicable butterflies to a separate model copy and holds other focal butterflies unchanged. The code first verifies all 239 original species×region native/contemporary/added counts and the global 14,553 / 206 totals against frozen primary figure data, then computes changes in native and contemporary resource units. If the frozen 239 baseline fails, the script stops rather than silently revising source data. The result must remain an exploratory, geographically biased sensitivity; it cannot be promoted to a revised corrected global GEB effect without source-specific adjudication, broader coverage and independent inference.

**Main GEB manuscript, frozen figures and PR #38 unchanged.**

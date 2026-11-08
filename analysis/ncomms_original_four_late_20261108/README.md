# Nature Communications: four genuine 2018–2022 first-indexed butterfly events recovered from the original-year archive

**8 November 2026. SOURCE-VERIFIED PARTIAL CHRONOLOGY; not the full original Parquet or a population-level ecological mechanism test.**

## What was recovered

The 83,573-row original species×island×year source (SHA256 `8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f`) was used by the author on 2026-10-30, but was intentionally not committed into GitHub. Its contemporaneously archived **five preselected prior-host links** retain the precise indexed first year and all recorded years for 10 original species×island endpoints.

The specific [source JSON](../../provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/full_original_island_chain_years_observed_20261030.json) is pinned by Git blob `68de667c1f93c4e7d99f25ffba54b341d7f1f697`. Its six 2011–2017 source endpoints agree exactly with the independently stored [5,665 original butterfly first-index/follow-up rows](../../provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv); the other four are source-verified **additional 2018–2022 positive event identities**.

| Butterfly species | Original island/landmass group ID | Indexed first year | Original candidate prior host in selected source |
|---|---:|---:|---|
| *Colias croceus* | 169 | **2019** | *Medicago sativa* |
| *Colias hyale* | 169 | **2019** | *Medicago sativa* |
| *Pieris rapae* | 951 | **2019** | *Lepidium didymum* |
| *Nymphalis antiopa* | 715 | **2021** | *Alnus incana* |

These are four **record-index** identities, not confirmed insect immigration/breeding or host plant consumption. The original five links were chosen *because* they had pre-index candidate host correspondence, so neither their resource association nor their chronological order estimates a population-level host effect or independent stepping-stone pathway. Two *Colias* links even share the same island pair and plant taxon.

## Source visibility after recovery

| Original 225 butterflies × 51 island source | Source-count |
|---|---:|
| Species×island possible index cells | 11,475 |
| Pre-2011 indexed butterfly cells, original marginal total (individual ids missing) | **524** |
| 2011–2017 indexed butterfly first-records with original species/island/year ids | **302** |
| New 2018–2022 first-indexed butterfly cells NOW independently reidentified from five-link archive | **4** |
| Total source-verified 2011–2022 new indexed identities | **306** |
| Original 2018–2022 later first-indexed marginal count | **128** |
| **2018–2022 later case identities still unrecovered** | **124** |
| Total 2011–2022 first indexed cells in original source (from both independent margins) | **430** |

The regenerated 11,475-row source-visibility CSV sets these four additional cells to `CONFIRMED_FIRST_INDEX_2018_2022_ORIGINAL_YEAR_SOURCE` and stores their exact years. It intentionally does **not** label any of the remaining unknown combinations as an ecological absence, and does not assume historical HOSTS presence/absence for other butterflies. The five preselected source hosts are stored only in a prominently selection-qualified column.

## Verification and reproducibility

- [Archived original source row chronology](../../provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/full_original_island_chain_years_observed_20261030.json) was originally validated against the privately provided Parquet in PR #104. Our current workflow re-verifies the **saved selected records and SHA**, not the inaccessible raw bytes.
- [Original 225×51 visibility builder](../../scripts/build_ncomms_island_risk_visibility.py), [four-case overlay](../../scripts/recover_ncomms_four_original_late_year_cases.py), [exact/corruption unit tests](../../tests/test_ncomms_four_original_late_years.py) and [dedicated workflow](../../.github/workflows/chocho-ncomms-four-original-year-recovery.yml).
- Source-checked GitHub workflow: [full reconstruction and independent 11,475-row audit](https://github.com/zuizui0223/chocho/actions/runs/42456536409). The artifact holds the enlarged visibility CSV, exact four recovered case CSV and pinned receipt.
- Other GEB and JBI research sources/manuscripts are unchanged; NC source-dependent tests are excluded from the independent GEB anonymized export.

## Ready-to-run full-original importer (not executed on the missing real bytes)

- [SHA-locked original raw Parquet importer](../../scripts/import_ncomms_exact_original_year_parquet.py) and [eight synthetic-only regression tests](../../tests/test_ncomms_exact_original_year_importer.py) were added to this Draft branch.
- The importer requires an **authorized local path** and checks the original SHA-256 before loading anything. It then checks the exact six original columns, **83,573 unique rows, 224 original species, 3,384 island IDs, 2001–2024**, original 14,425,536 precise recorded individuals, both positive source-year cohorts and the original 225×51 species/island margins. It never uploads the raw Parquet.
- Its output, when the original becomes readable, will be a derived, explicitly index-based 225×51 table of historical first-index years and verified recording-window status—not true occupancy or colonization.
- [Synthetic-only SHA/reconciliation pipeline](https://github.com/zuizui0223/chocho/actions/runs/42456772018) passed; this CI **does not claim to have read the original file**. It tests fail-closed behavior on deliberately synthetic rows.
- Run only after real bytes become accessible: `python scripts/import_ncomms_exact_original_year_parquet.py --original-year-parquet /authorized/original/global_original_species_island_year.parquet --output-dir /private/results`.

## Original file status and next step

ChatGPT Library currently indexes two copies of the matching original filename and 411,552-byte declared size under the project's butterfly and frog folders. Attempts to materialize or read either copy in this conversation were **blocked by file authorization**. Therefore neither copy has been byte-loaded or newly hash-checked and no claim of actual complete-Parquet recovery is made.

The originally known **83,573 rows; 224 species; 3,384 island groups; years 2001–2024; six columns; SHA256 8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f** remain the exact validation lock for a future authorized byte copy. Until the original can be attached to an accessible current conversation, the 524 pre-2011 and remaining 124 2018–2022 individual case identities cannot be exactly and safely regenerated from marginal counts alone.

The Nature Communications hypothesis `time-prior species-specific host resource × marine accessibility -> later butterfly record response` remains untested because the public record-index risk set and complete dated HOSTS identities cannot be reconstructed. Historic **positive** island-resource, port-flora, and host-forecast results remain valid at their separate original populations and are not repudiated.

**No journal submission, branch merge or false-negative model was performed.**

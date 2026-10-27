# Independent larval-host evidence audit for the 13 butterfly introduction events

**2026-10-27, frozen before compiling the new independent primary-source stage matrix.** Exploratory follow-up to chocho #78, independent of the original HOSTS Location×WCVP host pair-overlap null. This protocol does NOT retroactively validate the 145 proposed newly adopted hosts, infer an invasion rate from case reports, or modify the frozen GEB manuscript.

## Fixed target cohort and evidence source contract

- Fix all **13** prior literature-curated butterfly×receiving-area events from `data/exploratory/invasive_butterfly_verified_entry_events_20261027.csv`, representing **8 butterflies**, of which **5** are in the original 239 and **3** are external targeted additions. Do not add taxa because they give favorable results or merge off-panel 3 into the original 239-source summary.
- Find *independent* local source accounts that document an identified **consumer butterfly**, its **receiving locality**, a host taxon to the highest warranted rank and **life stage**. Only explicitly local evidence can count toward event-compatible larval use. A Europe-wide review cannot be silently attached to the Mallorca event; an experiment in Italy cannot be silently promoted to a wild Spanish observation.
- Record original research author/date/DOI/URL (not just the GloBI/HOSTS intermediate provider). Missing or conflicting citations remain unresolved, not zero. Never consider a duplicate of the original HOSTS resource.csv to be independent validation without comparing original source citation provenance.
- Botanical origin is a separate source, ideally exact accepted Kew/WCVP species and original recipient WGSRPD region. Use `NATIVE_AT_RECIPIENT`, `ALIEN_NATURALIZED_RECORDED`, `ALIEN_CULTIVATED_NOT_ESTABLISHED`, `PLANT_STATUS_UNRESOLVED`. Horticultural presence alone is **NOT** WCVP naturalization. For genus-level `Pelargonium spp.`, do not claim an exact species-level WCVP introduction.

## Explicit independent evidence tiers

- `FIELD_LARVA_HOST_REARED`: locally collected larva on a named host, reared to pupa or adult on the same documented taxon, with field→captive transition transparent.
- `FIELD_EGG_AND_LARVA_HOST`: locally observed eggs and larvae on named host without documented offspring completion; not proof of viable adult recruitment.
- `FIELD_LARVAL_DAMAGE_OR_OVIPOSITION`: directly documented larval damage or eggs on host at recipient, no feeding/rearing outcome independently established.
- `CONTROLLED_ONLY_LARVAL_COMPLETION`: controlled no-choice/choice or lab host feeding to adult, **NOT** observed wild host use.
- `REGIONAL_HOST_LIFE_CYCLE_REVIEW`: official regional host/life cycle description with no cited focal field individual; supportive context only.
- `ADULT_OR_ESTABLISHMENT_ONLY`: direct adult occurrence/genetics or invaded pest status without a specific local juvenile-host connection.
- `NO_SOURCE_VERIFIED_LOCAL_HOST_EVIDENCE`: insufficient literature located (not evidence no interaction exists).

The strongest tier per original event is reported **separately** from stage records at different places. Multiple evidence records from the same paper/citation are NOT independent replications.

## Primary decision metrics (frozen before retrieval)

1. Count original 13 receiving-area events with direct, independently source-verified juvenile×host field evidence (`FIELD_LARVA_HOST_REARED`, `FIELD_EGG_AND_LARVA_HOST`, `FIELD_LARVAL_DAMAGE_OR_OVIPOSITION`) and distinct species. A literature panel must never be extrapolated to 239 butterflies or interpreted as representative.
2. For those with documented host taxon botanical origin, report `NATIVE_AT_RECIPIENT` vs `ALIEN_NATURALIZED_RECORDED` vs `ALIEN_CULTIVATED_NOT_ESTABLISHED`; for uncertain genus-level plants keep unresolved. These are qualitatively different routes; do not call use of a recipient-native plant de novo evolved host adaptation.
3. Report independently verified controlled host capacity unsupported by local wild observations as a separate outcome. A particularly informative counterexample is experimentally surviving on native European `Geranium` but predominantly using ornamental `Pelargonium` in wild European `Cacyreus marshalli` source reports.
4. **Chronology conflict check:** `Papilio demoleus` Dominican Republic first scientifically documented specimen **29 March 2004**, but EPPO's later distribution page says **March 2001**; do not overwrite the historical source as one confirmed arrival year. For `Pieris rapae` NZ, original 1966 reference labels successful arrival 1930 while a modern NZ butterfly organization reports first local observation 1929; keep first sighting vs establishment separately and explicitly unresolved.
5. Source publication year is distinct from original field observation year, first invasive detection, first host observation, experimental collection year or laboratory replicate. Unknown year remains NULL.
6. Maintain a native-host study comparator where available, but NO native-source absence test from an unsearched paper corpus.

## Stop conditions and outputs

An exact field host name, locality, developmental stage, full source citation and plant-origin source are required for a high-grade mechanistic route. A location-ambiguous review can support discussion but never count as proof of field host use for a named receiving event.

The result may identify contrasting confirmed **mechanisms of interaction realization** but CANNOT estimate their global frequencies or explain new climatic niche expansion. Keep data quality/source overlap auditable, preserve zero/unknown distinction, and do not enter GEB/main/PR #38 or PR #78 result tables.

Do not copy copyrighted article text or raw interaction matrices; capture only short, independently written evidence descriptions and source URLs.

# Pre-outcome protocol: source-positive country cells to butterfly-point climate

**Frozen 2026-10-27 before opening NEW point outcomes for climate-novel or climate-analogue cases.** This is a source-audited extension of chocho #72 and #73, not a pre-registered prospective hypothesis. No modifications to main or GEB submission (#38, #70).

## Central question

The earlier #72 country-level result includes 73 **species×country positive recording cells**, 32 located in countries whose introduced-host resource *area* is climate-novel-dominant and 41 in climate-analogue-dominant countries. These are NOT 73 individual insects. Can the actual indexed butterfly observation points be located in climatic environments that independently corroborate that coarse *host-resource opportunity class*? Can observed butterfly points actually fall within the introduced-host geographic footprint, rather than merely sharing a country?

## Stage 1: derive and freeze exact source-cohort identities, before touching new point outcomes

Read exactly the five original Git-SHA-frozen files from #72: the 235×140 source-conditional species-country fraction table, the original country design, original GBIF country-count CSV, 239-taxon manifest and source GBIF query receipts. Call #72 `load_areas`, `load_observations`, `opportunity_class`, `PRIMARY` and reproduce its selection exactly:
- complete noncolliding accepted GBIF species, no WCVP native-host guard, country added-host area >=5%, >=80% of added area either native-host BIO1/BIO12 5–95% analogue or outside it;
- >=3 *other* recorded butterfly taxa in the same country (no own-target effort leakage);
- positive if the frozen 2010–2025 accepted-bases species×country GBIF index count >=1; only butterfly×WGSRPD Level1 strata with BOTH host-resource opportunity classes survive.
- Expected fixed source positives exactly **73 = 32 novel + 41 analogue**, with exactly **20 species**. Retain each case species, original accepted GBIF `taxonKey`, country, original frozen pooled index count, native-host region climate source class, analogue/novel host-resource area fractions, original matched Level1 group. No location or individual source metadata in the frozen roster.
- Summarise frequency of long-tailed positive counts and assess before retrieval whether finite source-specific query caps can cover all 73 cases; do not choose positive-rich species, favor novel or select cases on live points. New live GBIF counts can drift and must be audited rather than silently reconciled.

## Stage 2: georeferenced GBIF point recording diagnostic (NOT colonization)

Use only the 73 frozen cases. Request accepted-source taxon, ISO country, year 2010–2025, original 5 GBIF `basisOfRecord`, `hasCoordinate=true`, `hasGeospatialIssue=false`, `occurrenceStatus=PRESENT`. Page every eligible query without sampling (300 rows/page), cap 10,000 rows per case and 40,000 records total with explicit unscorable overflow; limit API to ~1 request/sec, bounded retries and timeout. Keep original vs live index drift and potential duplicates; validate accepted speciesKey/taxonKey, country, year, finite coordinates, reported uncertainty, life-stage metadata. Count per case, class and butterfly species using both equal-case and record-weighted summaries, while explicitly not treating individual records as biological replicates.

Stage 2 does **not** alone calculate climate, because coordinates without same-source matched high-resolution BIO1/BIO12 layers would not prove anything about the insect climatic niche. A geographic country-level source class must never be assigned to all points in that country. Do not present Stage 2 as proof of realized niche entry.

## Stage 3: independent consistent climate audit

Before calling a point climate shift, use a *single coherent* climate dataset for both observed butterfly points and the comparator/baseline. First evaluate availability of source-aligned 1981–2010 CHELSA v2.1 BIO1/BIO12 grids plus verified WGSRPD polygons; if unavailable, use a clearly labelled independent WorldClim v2.1 10-arcminute BIO1/BIO12 **at both observation points and native-host polygon representative points** (not CHELSA host thresholds applied directly to WorldClim point values). Predefine 5–95% BIO1/BIO12 rectangle, elevation limitations, spatial autocorrelation, and explicit estimation unit. Never use contemporary climate at one location against mismatched-era or mismatched-model thresholds without sensitivity/matching.

## Stage 4: scientific gate

At minimum require all 73 source cohort identities verified (never substitute missing taxa with zeros), all selected point queries complete with original-to-live drift audit, and a measurable consistent climate source. Even then, without independent wild reproduction, documented host use, introduction chronology, survey effort and repeated detections, call results **documented point-climate recording patterns** or *potential host-resource niche opportunity*, not butterfly adaptation, colonization, population growth, or novel interaction competition.

Do not merge the result into the original GEB manuscript or release author-only raw botanical geometry, source coordinates, GBIF observer identities or exact localities in the code archive. Derived species×country quantities and checksums are permissible.

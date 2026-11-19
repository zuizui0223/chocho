# 104-island host-masked alien flora negative control — verified 19 Nov 2026

[Protocol](PROTOCOL.md) · [Verified condensed score summary](VERIFIED_HOST_MASKED_SUMMARY_20261119.json) · [GitHub Actions run #44445904218](https://github.com/zuizui0223/chocho/actions/runs/44445904218) · [Original full JSON artifact #14493153837](https://github.com/zuizui0223/chocho/actions/runs/44445904218/artifacts/14493153837)

**Source audit passed:** 5,355 original WCVP coastal whole-flora dyads exactly regenerated (including two response-undefined dyads). Models evaluate 5,353 resource-defined dyads and 1,088 strict both-islands-unseen held-out pairs. Raw inputs are SHA-locked to archived HOSTS×WCVP / PR118. 7,887 alien accepted plant taxa: 668 known host IDs of the frozen 239-butterfly panel (H), 7,219 other unlisted taxa (O).

| Extra introduced plant-composition predictor | E held-dyad MSE improvement | S retrospective improvement |
|---|---:|---:|
| W all alien flora | +26.363% | +33.768% |
| **H only registered HOSTS hosts** | **+33.821%** | **+41.804%** |
| O only other introduced plants | +19.539% | +26.112% |
| H+O | +33.707% | +41.584% |

- Add H after O: **+17.609% E** [descriptive fivefold 95% +8.539,+28.224], **+20.939% S** [+11.121,+31.503].
- Add O after H: **−0.173% E** [−0.481,+0.058], **−0.379% S** [−0.862,+0.133].
- Add H after W: +9.740% E, +11.672% S; add O after W: +7.583% E, +6.926% S. Several strongly shared predictors can add complementary Jaccard-weight information.
- H positive in all five held-out folds under E and S. All flags for empty H/O unions were zero in this matched cohort.

**Ecological interpretation:** host-flora composition is more informative about **potentially supported** butterflies than composition of taxa not catalogued as hosts. Yet the response is mathematically constructed from the same WCVP plant geography and HOSTS links used to define H, so this is **source structural validation, not independent host-specific feeding evidence**. O is a catalogue-exclusion negative control, not biologically demonstrated non-hosts; it is larger and has different taxon/richness profiles. Independent butterfly assemblage gains were much smaller, and no population or colonization effect has been measured.

**Integrity:** full original artifact SHA256 `3d9531806d25c0d08022ec63b99b7c7224d89303383139abe964b31aa85ceb62`; full JSON is preserved in the linked Actions artifact. The shorter JSON committed here retains the key numerical estimates and source linkage. The two successful runs (push and PR) are **duplicate executions of the same input**, not independent biological replications.

# Phenological host-window raw source access: feasibility decision

**2026-10-10** · Source availability audit complete, ecological analysis **not** run.

### Executed audit

Workflow: https://github.com/zuizui0223/chocho/actions/runs/38011851108  
GitHub Actions tests: **3 passed**. All audit steps completed and saved the original-source access receipt. This is **a software/audit success only**, not a positive scientific result.

Source: Toftegaard et al. (2018) https://doi.org/10.5061/dryad.jp328r6

| Original file | Attempted primary URL | Response | Gate |
| --- | --- | --- | --- |
| Plant phenology in the field.txt | https://datadryad.org/downloads/file_stream/79537 | HTTP **403** | SOURCE_ACCESS_BLOCKED |
| Sampling site areas.txt | https://datadryad.org/downloads/file_stream/79538 | HTTP **403** | SOURCE_ACCESS_BLOCKED |

Neither file's raw bytes or schema was verified. The repeated blocked result agrees with independent direct-web and container download attempts. The journal's public source-description establishes 2010–2013, three regions, six hosts and host phenology with egg/larva annotations, but it does **not** provide a verified original row sample or denominators. `source_bytes_complete=false`; `row_level_ecological_effect_estimated=false`.

**Decision:** STOP the local phenological complementarity pilot until the original files can be obtained legally through an available source or author-provided archive. Do not invent a result or select a proxy endpoint from the published abstract. Only then can the protocol's egg-negative denominator, phenology stage, flight timing and site-year validation be inspected.

**Scientific priority regardless of gate:** history-dependent host fallback as set out in `HOST_REDISRIBUTION_HYSTERESIS_HYPOTHESIS_V01.md`. That requires matched historical exposure, replicated native/introduced host choice, tracked adult emergence and host-removal contrasts. No present existing network metric supplies those variables.

**Isolation:** main GEB PR #38 was not modified by this branch.

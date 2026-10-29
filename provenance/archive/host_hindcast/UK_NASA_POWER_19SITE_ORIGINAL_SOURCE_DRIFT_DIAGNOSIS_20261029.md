# Independent NASA POWER original source drift diagnosis

2026-10-29. After the first UKBMS historical BSBI HOSTS model [Action #40806287237](https://github.com/zuizui0223/chocho/actions/runs/40806287237) **failed closed** on a full-raw-JSON SHA mismatch for the same original NASA POWER 1981–2021 data, before any late-test host model score was generated.

## Why and what is known
The previous NASA POWER 19-source history Action #40802452296 succeeded with 19 original daily 1981–2021 responses, one raw-body SHA per station; those SHA values protect both meteorological numeric values and dynamically generated response metadata. Re-querying the same 19 sites on 29 October changed the aggregate raw hash. This could be a provider-generated HTTP-JSON response wrapper, a source revision, changes in the original locality coordinates, or another original-system error. **No explanation is assumed** without source comparison.

## Diagnostic frozen before outcomes
- Retrieve the original historical NASA POWER source receipt artefact #40802452296, exactly 19 unlinked raw response hashes, original sites/UKBMS visit selection; refuse altered taxonomy or different year range.
- Retrieve each 1981–2021 station daily T2M/PRECTOTCORR using **identical source method** but do not fit a butterfly model. Compare how many original response SHA match, compute the canonical SHA of the 14,975 **numeric day values only** (remove source JSON wrapper and site coordinates); record missingness.
- Re-fetch exactly **one original station twice inside the same source run** to check whether identical daily values can have different HTTP response SHA even without a source date change. Never export original station identity, coordinates, daily temperature/precipitation values, biological observation rows, or per-location digests.
- Distinguish (1) raw-source body changed from (2) same-session numeric data changed. The earlier run did not persist numeric-only source hashes, so it is **not possible** to claim original numeric contents are identical unless another independently valid legacy source is available.
- If raw body differs but current canonical daily values are internally stable and physically valid, a subsequent **new independently source-versioned exploratory hindcast** may be run, fitting *all* climate and host comparators against exactly the same current daily source snapshot. The earlier run must not be claimed reproduced. Keep the 1987–95 training/2005–21 test dates, pre1987 plant-only predictors and biological source gates unchanged.
- If current numeric values are unstable within the same run, SOURCE UNAVAILABLE and **do not report ecological model results**. Rerunning blindly is not permitted.

Source NASA POWER original daily API, https://power.larc.nasa.gov/docs/services/api/temporal/daily/; original UKBMS and BSBI files SHA-pinned from PR #97 and this branch.

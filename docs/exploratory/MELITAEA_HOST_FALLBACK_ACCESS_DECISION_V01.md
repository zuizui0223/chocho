# Åland butterfly host fallback: source and novelty gate (2026-10-10)

**Status:** SOURCE ACCESS BLOCKED; INDEPENDENT TEST NOT PERFORMED. No new ecological effect.

## Concrete source

DiLeo et al. (2024), Dryad DOI https://doi.org/10.5061/dryad.905qfttrg, latest version dated 2024-07-25. Public README lists **36,704 patch × annual survey rows from 2004–2013** and categorical abundance (0–3) for **Plantago lanceolata** and **Veronica spicata**, butterfly patch occupancy and nest count. The two plant species are already used by *Melitaea cinxia* in the Åland Islands, **not** an alien-versus-native plant treatment.

These unusually valuable observations could, conditional on sufficient genuine plant declines, test a response of local butterfly occupancy after one larval host declines in the presence of another. The latter is an ecological endpoint independent of a mechanically computed resource envelope, but **not a randomized causal removal** and **not historical adoption of an invasive host**.

## Executed fixed feasibility attempt

- Protocol: `MELITAEA_TWO_HOST_FALLBACK_PROTOCOL_V01.json`, committed before inspecting original event counts.
- Audit: `scripts/audit_melitaea_host_fallback_panel.py`; latest validated CI https://github.com/zuizui0223/chocho/actions/runs/38014340894.
- **4 software tests passed** in the latest [exact-head audit](https://github.com/zuizui0223/chocho/actions/runs/38014340894), including source ZIP-path validation, duplicate patch-year rejection, exposure/outcome ordering and fail-closed HTTP access. An earlier audit exposed and led to correction of a duplicate-year key bug, and an intermediate synthetic fixture used incorrectly escaped newlines; both were repaired before this green receipt. These are software integrity checks, **not biological results**.
- Official Dryad v2 `/download`: **HTTP 401**; official version 3 direct file stream: **HTTP 403**. Neither original archive nor any of its 36,704 rows was obtained in the runner.
- As a consequence, `source_verified=false`, `effect_estimated=false`. Number of host declines, distinct patches, independent networks, year balance and conditional next-year extinction differences **remain unknown**.
- The audit script now additionally accepts an actual unchanged published version 3 ZIP as a local `--archive` argument, hashes it, requires the exact `empirical_models/data/RAWDATA/fall_survey_2004_2013.csv` path and validates time order and support before effect fitting. An automatically generated or synthetic CSV is not a substitute.

## Published results already covering generic host abundance

- Hanski & Singer (2001), https://doi.org/10.1086/321985: host composition and local host preference already shown to shape butterfly colonization.
- Opedal et al. (2020), https://doi.org/10.1002/ecy.3186: **19-year** Åland metacommunity dataset; colonization and extinction already modelled against abundances of both *Plantago* and *Veronica*, connectivity, and climate. Generic abundance–extinction associations are not novel. Existing evidence suggests *Veronica* can influence parasitoid interactions differently than *Plantago*, so functional backup is not guaranteed by the label “alternative host.”
- DiLeo et al. (2024) itself evaluated butterfly extinction risk with host-environment variables and population genetics. A repackaging of current host abundance effects is not new.

**Remaining question**, conditional on source access and adequate actual declines: does *recent change* in one host alter the other host's predictive relationship with subsequent **observed local persistence**, at matched pre-decline state and geographic sampling support? The answer can be positive, null or negative. A single-species observational time series cannot, by itself, identify historical eco-evolutionary hysteresis; that requires an experimental history × resource-loss design with adult offspring output.

## Decision

Do not fit the fixed model or claim plant-loss buffering while the source gate is blocked. Do not lower preregistered event/independence criteria to rescue a sparse sample. The experiment-targeted randomized resource-history `RANDOMIZED_HOST_HISTORY_FALLBACK_PROTOCOL_V01.json` remains a separate, potentially causal route. The GEB submission PR #38 is untouched.

Exact optional local follow-up, only **if** authentic original ZIP is available:

```bash
python scripts/audit_melitaea_host_fallback_panel.py \\
    --archive DiLeo_et_al_DRYAD_ver3.zip \\
    --receipt results/melitaea_host_fallback_source_gate_v01.json
```

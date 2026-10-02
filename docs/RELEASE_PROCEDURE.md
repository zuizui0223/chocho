# Public archival release procedure — v1.0.0-butterfly

This procedure is administrative only. The current v0.2 scientific claim map, species panels, thresholds, result receipts and manuscript statistics are the release boundary.

## Why this release uses a manual Zenodo draft

The public DOI should appear inside the exact tagged source snapshot. A manual Zenodo draft allows the DOI to be reserved before publication, so the reserved DOI can be written into the exact release commit before the Git tag is created.

## Inputs already complete

- Package and CITATION version: `1.0.0`
- Intended Git tag: `v1.0.0-butterfly`
- Current v0.2 claim map, manuscript, supplement and review-package builders
- Current paper CI validating v0.2 tests and four-display figure regeneration
- Release preflight: `python scripts/release_preflight.py`

## Final sequence

1. **Choose the software license**
   - Add `LICENSE`.
   - Use the chosen SPDX identifier consistently in archive metadata.

2. **Finalize identifying metadata**
   - Final author list/order.
   - Affiliations and corresponding author.
   - ORCIDs where available.
   - CRediT roles.
   - Acknowledgements.
   - Funding.
   - Conflict-of-interest statement.

3. **Create a Zenodo draft manually**
   - Select resource type: software.
   - Enter the final title and creators.
   - Do not publish yet.
   - Reserve a DOI in the draft.
   - Keep the draft; deleting it loses the reserved DOI.

4. **Insert the reserved DOI into the release candidate**
   - Add `doi: <reserved DOI>` to `CITATION.cff`.
   - Insert the same DOI in the separate v0.2 title page.
   - Insert the same DOI in the final public Data and Code Availability statement.
   - Do not put the public DOI in the double-anonymous review manuscript.

5. **Keep peer-review access separate from the public release**
   - Initial GEB peer review uses the de-identified reviewer ZIP uploaded directly as supplementary review material.
   - No external anonymous reviewer URL is required when that ZIP is attached to the submission.
   - Do not expose the public author-identifying archive in the blinded manuscript.

6. **Run release preflight**
   - `python scripts/release_preflight.py`
   - It must report `READY`.

7. **Run final CI on the exact release commit**
   - Current v0.2 offline paper test suite.
   - Current three-main-plus-one-supplement figure regeneration.
   - Blinded v0.2 DOCX build/anonymity checks.
   - Anonymous reviewer bundle build/identity scan/tests/current-figure regeneration.

8. **Create the Git tag**
   - Tag the exact green release commit as `v1.0.0-butterfly`.
   - Do not move or reuse the tag later.

9. **Create an archive from that exact tag**
   - Export the exact tagged source archive.
   - Confirm that it contains the current v0.2 manuscript/claim map, frozen inputs/result receipts and reserved DOI metadata.

10. **Upload the exact tag archive to the existing Zenodo draft**
    - Upload the archive to the draft that already owns the reserved DOI.
    - Check version = `1.0.0`.
    - Check creators/order, license, title, description and keywords.
    - Preview the Zenodo record.

11. **Publish the Zenodo record**
    - Publishing registers/activates the reserved DOI.
    - Verify that the DOI resolves and the archived file corresponds to `v1.0.0-butterfly`.

12. **Final public-version check**
    - Identifying/public metadata: public stable-archive DOI.
    - Final public Data and Code Availability statement cites that archive.
    - The double-anonymous review manuscript remains free of author-identifying repository links.
    - Final proofread and exact-head CI confirmation.

## Do not do

- Do not use a placeholder DOI as though it were real.
- Do not create the final tag before the reserved DOI and author/license metadata are in the release commit if the tagged source is intended to contain those values.
- Do not expose the public author-identifying archive in the blinded manuscript.
- Do not reintroduce the v0.1 five-figure submission route.
- Do not add response-driven ecological analyses to satisfy release administration.

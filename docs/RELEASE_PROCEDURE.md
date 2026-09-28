# Final release procedure — v1.0.0-butterfly

This procedure is administrative only. The scientific claim map, species panels, thresholds, hypotheses, result receipts and manuscript statistics are frozen.

## Why this release uses a manual Zenodo draft

For this release, the public DOI should appear inside the exact tagged source snapshot. Zenodo allows a DOI to be reserved before publication for a manual draft deposit, but Zenodo's GitHub integration does not allow pre-reserving a DOI before the GitHub release.

Therefore, use a **manual Zenodo draft deposit** for v1.0.0-butterfly.

## Inputs already complete

- Package and CITATION version: `1.0.0`
- Intended Git tag: `v1.0.0-butterfly`
- Byte-exact frozen figure inputs: `data/frozen/figure_sources/`
- Blinded review DOCX builder and anonymity checks
- De-identified anonymous reviewer code bundle builder and identity scan
- Release preflight: `python scripts/paper/release_preflight.py`

## Final sequence

1. **Choose the software license**
   - Add `LICENSE`.
   - Use the chosen SPDX identifier consistently in any archive metadata.

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
   - Insert the same DOI in the separate title page.
   - Insert the same DOI in the final public Data and Code Availability statement.
   - Do not put the public DOI in the double-anonymous review manuscript.

5. **Create the anonymous reviewer URL**
   - Upload the already-built de-identified reviewer ZIP to a host whose landing page and file metadata do not identify the authors.
   - Insert that URL in place of `[ANONYMIZED REVIEW LINK]`.
   - Rebuild the final blinded DOCX and rerun anonymity/metadata checks.

6. **Run release preflight**
   - `python scripts/paper/release_preflight.py`
   - It must report `READY`.

7. **Run final CI on the exact release commit**
   - Paper test suite.
   - Artifact-free five-figure regeneration.
   - Blinded DOCX build/anonymity checks.
   - Anonymous reviewer bundle build/identity scan/tests/figure regeneration.

8. **Create the Git tag**
   - Tag the exact green release commit as `v1.0.0-butterfly`.
   - Do not move or reuse the tag later.

9. **Create an archive from that exact tag**
   - Download/export the exact tagged source archive.
   - Confirm that it contains the frozen inputs and the reserved DOI metadata.

10. **Upload the exact tag archive to the existing Zenodo draft**
    - Upload the archive to the draft that already owns the reserved DOI.
    - Check version = `1.0.0`.
    - Check creators/order, license, title, description and keywords.
    - Preview the Zenodo record.

11. **Publish the Zenodo record**
    - Publishing registers/activates the reserved DOI.
    - Verify that the DOI resolves and that the archived file corresponds to `v1.0.0-butterfly`.

12. **Final submission check**
    - Public/title-page metadata: public Zenodo DOI.
    - Blinded manuscript: anonymous reviewer URL, not the public author-identifying archive.
    - Final 20-page DOCX rebuilt after URL insertion.
    - Final proofread.

## Do not do

- Do not use a placeholder DOI as though it were real.
- Do not create the final tag before the reserved DOI and author/license metadata are in the release commit if the tagged source is intended to contain those values.
- Do not expose the public author-identifying archive in the blinded manuscript.
- Do not add response-driven ecological analyses to satisfy release administration.

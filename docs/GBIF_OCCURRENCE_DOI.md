# GBIF occurrence DOI for the occurrence-validation panel

## Why a DOI is required

The butterfly occurrence validation used records retrieved synchronously through the GBIF occurrence search API. Search-API results do not receive an automatic DOI. GBIF citation guidance recommends citing GBIF-mediated data using a DOI.

The exact analysis set has been recovered from the frozen historical workflow artifact:

- source repository: `zuizui0223/TTF`
- workflow run: `36270581743`
- artifact: `10920120240` (`butterfly-climate-release-postgate-summary-v01`)
- source file: `occurrences.csv`
- records: **53,434**
- unique GBIF IDs: **53,434**
- butterfly species with records: **31**
- contributing GBIF datasets: **313**
- `occurrences.csv` SHA-256: `1dc71938fd5fbed48756cc8bfd0e4bb9580f7d20eebfb652547febcc2893dcad`

The datasetKey contribution counts are frozen at
`data/frozen/gbif/gbif_dataset_key_record_counts_v0.1.csv`.

## Preferred DOI route: exact GBIF-ID occurrence download

GBIF currently supports occurrence downloads filtered by `GBIF_ID`, and a single download predicate may contain up to 101,000 items. The 53,434-record analysis set therefore fits in one exact-ID download.

This route is preferable to registering a derived dataset here because it generates a standard GBIF occurrence-download DOI without requiring an external persistent source URL.

### 1. Recover the historical occurrence artifact

Download GitHub Actions artifact `10920120240` from TTF workflow run `36270581743`. The expected ZIP contains `occurrences.csv`.

### 2. Generate the exact download request

```bash
python scripts/prepare_gbif_occurrence_doi_request.py \
  --occurrences butterfly-climate-release-postgate-summary-v01.zip \
  --email YOUR_GBIF_EMAIL \
  --output gbif_occurrence_download_request.json
```

The script refuses to proceed unless the recovered `occurrences.csv` has the expected SHA-256 and exactly 53,434 unique GBIF IDs.

### 3. Submit the request to GBIF

A GBIF.org account is required. Use the GBIF username, not the email address, for HTTP Basic Authentication.

```bash
curl --include \
  --user YOUR_GBIF_USERNAME:YOUR_GBIF_PASSWORD \
  --header "Content-Type: application/json" \
  --data @gbif_occurrence_download_request.json \
  https://api.gbif.org/v1/occurrence/download/request
```

Save the returned download key.

### 4. Wait for GBIF to finish the download and record the DOI

```bash
curl -Ss https://api.gbif.org/v1/occurrence/download/YOUR_DOWNLOAD_KEY
```

When status is `SUCCEEDED`, the response contains the download DOI and download URL.

### 5. Update the manuscript

Replace `GBIF_OCCURRENCE_DOWNLOAD_DOI_PLACEHOLDER` everywhere with the DOI URL, e.g.
`https://doi.org/10.15468/dl.xxxxx`.

The submission preflight must reject the placeholder and accept only a real GBIF occurrence-download DOI.

## Citation wording

Use the citation supplied by the GBIF download landing page as authoritative. The working manuscript template is:

> GBIF.org. 2026. GBIF occurrence download for the butterfly occurrence-validation record set. https://doi.org/GBIF_OCCURRENCE_DOWNLOAD_DOI_PLACEHOLDER

Do not substitute the generic GBIF website citation for the data DOI.

## Fallback

If an exact-ID occurrence download cannot be issued, register a GBIF derived dataset using the 313 datasetKey/count pairs in `data/frozen/gbif/gbif_dataset_key_record_counts_v0.1.csv`. GBIF derived-dataset registration requires authentication and a persistent URL for the extracted data, so the exact-ID occurrence download is simpler for this analysis.

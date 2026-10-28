# CMIP6 monthly raster layout correction (no ecological retuning)

2026-10-28. Initial [GitHub Actions run #40611056393](https://github.com/zuizui0223/chocho/actions/runs/40611056393) **successfully downloaded all 18 required source archives**, and successfully processed the WorldClim 2.1 **baseline 1970–2000 monthly GeoTIFFs**. It then stopped when the first future CMIP6 GeoTIFF arrived with a different archive layout.

The original parser assumed 12 one-band TIFF members per ZIP for both historical and future. WorldClim **CMIP6 future** 10-minute ZIPs instead contain **one 12-band TIFF** encoding twelve calendar months. Correct the *input decoder only*: accept either exactly 12 single-band TIFFs identified by ending month or exactly one 12-band TIFF with band 1–12 corresponding to January–December. Reject any other combination, mismatched raster origins, missing pixels, temperature unit encodings, or source species/coordinates. Fixed source years (2050, 2090), two GCMs, SSPs, 11.5/33°C, 323 degree-days, host-leaf persistence scenarios and original monarch larval photo filters are unchanged.

No biological results were opened or used to select alternatives; the initial run did not compute a single future scenario. Offline regression tests pass after the parser correction. The next same-source run is the first eligible climate analysis.

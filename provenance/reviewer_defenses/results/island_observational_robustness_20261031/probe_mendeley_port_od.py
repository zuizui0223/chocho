#!/usr/bin/env python3
"""Discover actual publishable port-to-port source file identities (metadata only).

Verschuur et al. (2022), DOI 10.17632/vzzy3b9gg4.1. Never convert model
allocated freight, sea-network graph segments or country-pair LSBCI into
observed port voyages. No bulk route bytes are downloaded in this probe.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = "https://api.data.mendeley.com/datasets/vzzy3b9gg4"
URLS = {
    "dataset": ROOT + "?version=1",
    "files": ROOT + "/files?version=1",
    "public_files": "https://api.data.mendeley.com/datasets/publics/vzzy3b9gg4/files?version=1",
    "folders": ROOT + "/folders?version=1",
}
MAX_ITEMS = 500


def probe(url):
    req = Request(url, headers={
        "Accept": "application/json",
        "User-Agent": "chocho-science-source-audit/1.0 (public metadata only)",
    })
    try:
        with urlopen(req, timeout=20) as response:
            raw = response.read(500000)
            data = json.loads(raw)
            return {"http_status": response.status, "data": data}
    except HTTPError as e:
        return {"http_status": e.code, "error": str(e.reason)}
    except (URLError, ValueError, TimeoutError) as e:
        return {"http_status": None, "error": type(e).__name__ + ": " + str(e)[:350]}


def summarize(value):
    if isinstance(value, dict):
        # Metadata should not leak transient S3 presigned download URLs.
        fields = {k: value.get(k) for k in (
            "id", "name", "filename", "folder_id", "size",
            "version", "doi", "description", "total_count", "count"
        ) if k in value}
        for k in ("files", "results", "items", "folders", "data"):
            if isinstance(value.get(k), list):
                fields[k] = [
                    {n: x.get(n) for n in ("id", "name", "filename",
                                           "folder_id", "size") if n in x}
                    for x in value[k][:MAX_ITEMS] if isinstance(x, dict)
                ]
                fields[k + "_reported_count"] = len(value[k])
        return fields
    if isinstance(value, list):
        return {"items": [
            {k: x.get(k) for k in ("id", "name", "filename", "folder_id", "size")
             if k in x}
            for x in value[:MAX_ITEMS] if isinstance(x, dict)],
            "reported_count": len(value)}
    return {"kind": type(value).__name__}


def main():
    results = {}
    for key, url in URLS.items():
        r = probe(url)
        results[key] = {
            "http_status": r["http_status"],
            "error": r.get("error"),
            "summary": summarize(r.get("data")) if "data" in r else None,
        }
    found = []
    for entry in results.values():
        s = entry.get("summary") or {}
        for key in ("items", "files", "results", "folders"):
            found += [x for x in s.get(key, []) if isinstance(x, dict)]
    interesting = sorted({
        (x.get("filename") or x.get("name") or "").strip()
        for x in found
        if any(t in (x.get("filename") or x.get("name") or "").lower()
               for t in ("port", "network", "trade", "flow"))
    })
    result = {
        "schema": "public_mendeley_global_port_od_metadata_probe_v1",
        "source": "Verschuur, Koks and Hall 2022",
        "source_doi": "10.17632/vzzy3b9gg4.1",
        "license_claimed_by_dataset_page": "CC BY 4.0",
        "recorded_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "probe_state": ("PUBLIC_METADATA_AVAILABLE"
                        if any(x["http_status"] == 200 for x in results.values())
                        else "PUBLIC_API_NOT_ACCESSIBLE_FROM_RUNNER"),
        "metadata": results,
        "interesting_names_not_yet_data": interesting,
        "limitations": [
            "This probe reads metadata only; no route edge or butterfly association is tested.",
            "Port-to-port modeled trade/flows are not directly observed AIS vessel trips.",
            "2022 model network cannot identify 2001-2010 transported host plants.",
            "Missing public API access is an access result, not evidence of zero shipping.",
        ],
    }
    path = Path(
        "provenance/reviewer_defenses/results/island_observational_robustness_20261031/"
        "port_to_port_source_discovery_20261031.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf8")
    print(json.dumps({k: v for k, v in result.items()
                      if k != "metadata"}, ensure_ascii=False))


if __name__ == "__main__":
    main()

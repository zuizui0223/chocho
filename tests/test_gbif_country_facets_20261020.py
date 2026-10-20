"""Contract tests for full-index COUNTRY facets, never occurrence samples."""
from __future__ import annotations

import importlib.util
import json
import threading
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import parse_qs, urlparse

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/acquire_gbif_country_facets_20261020.py"


def module():
    assert SCRIPT.is_file(), "The country-facet acquisition implementation is missing"
    spec = importlib.util.spec_from_file_location("country_facets_20261020", SCRIPT)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def original_manifest():
    rows = []
    for index in range(239):
        status = (
            "COMPLETE" if index < 135 else
            "PARTIAL_FETCH_NOT_ELIGIBLE" if index < 231 else
            "TAXON_UNRESOLVED"
        )
        name = "Chlosyne lacinia" if index == 135 else f"Butterfly species{index:03d}"
        key = 1914465 if index == 135 else 10000 + index
        rows.append({
            "species": name, "status": status,
            "taxon": {"key": key} if index < 231 else {"reason": "NO_EXACT_SPECIES_MATCH"},
        })
    rows[229].update(species="Oeneis bore", taxon={"key": 9804158})
    rows[230].update(species="Oeneis polixenes", taxon={"key": 9804158})
    return {
        "schema": "chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01",
        "species_requested": 239,
        "species_statuses": {
            "COMPLETE": 135, "PARTIAL_FETCH_NOT_ELIGIBLE": 96,
            "TAXON_UNRESOLVED": 8,
        },
        "species_receipts": rows,
    }


def payload(count=10, countries=(('US', 7), ('JP', 2))):
    return json.dumps({
        "offset": 0, "limit": 0, "count": count, "results": [],
        "facets": [{"field": "COUNTRY", "counts": [
            {"name": country, "count": number} for country, number in countries
        ]}],
    }).encode()


def test_exact_original_names_and_keys_define_all_924_queries(tmp_path):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    source = api.load_original_manifest(path)
    queries = api.build_queries(source)
    assert len(queries) == 924
    assert len({q["species"] for q in queries}) == 231
    assert len({q["taxon_key"] for q in queries}) == 230
    assert len({q["query_id"] for q in queries}) == 924
    assert {q["taxon_key"] for q in queries if q["species"] == "Chlosyne lacinia"} == {1914465}
    assert not ({r["species"] for r in source["species_receipts"] if r["status"] == "TAXON_UNRESOLVED"}
                & {q["species"] for q in queries})
    for query in queries:
        params = parse_qs(urlparse(query["url"]).query)
        assert params["taxonKey"] == [str(query["taxon_key"])]
        assert params["limit"] == ["0"]
        assert params["facet"] == ["country"]
        assert params["facetLimit"] == ["300"]
        assert params["year"] in (["2010,2017"], ["2018,2025"])
        assert set(params["basisOfRecord"]) == set(api.BASIS_SCOPES[query["basis_scope"]])


def test_missing_source_key_and_status_count_drift_fail_before_network(tmp_path):
    api = module()
    data = original_manifest()
    data["species_receipts"][135]["taxon"] = None
    path = tmp_path / "source.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="taxon key.*Chlosyne lacinia"):
        api.load_original_manifest(path)
    data = original_manifest()
    data["species_receipts"][135]["status"] = "COMPLETE"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="status counts"):
        api.load_original_manifest(path)


def test_country_counts_retain_missing_country_remainder():
    result = module().parse_country_facet(payload())
    assert result == {
        "reported_count": 10, "country_count": 2,
        "country_counts": {"JP": 2, "US": 7},
        "country_records": 9, "missing_country_count": 1,
    }


def test_truncated_json_never_becomes_a_complete_query():
    api = module()
    with pytest.raises(api.FacetValidationError, match="JSON"):
        api.parse_country_facet(payload()[:-4])


def test_facet_limit_is_not_mistaken_for_complete_country_coverage():
    api = module()
    countries = [(chr(65 + n // 26) + chr(65 + n % 26), 1) for n in range(300)]
    with pytest.raises(api.FacetValidationError, match="facetLimit"):
        api.parse_country_facet(payload(300, countries))


@pytest.mark.parametrize("count,countries", [
    (10, (("US", 6), ("US", 1))),
    (10, (("USA", 1),)),
    (10, (("us", 1),)),
    (10, (("US", -1),)),
    (10, (("US", True),)),
    (10, (("US", 11),)),
    (True, (("US", 1),)),
])
def test_invalid_facets_fail_closed(count, countries):
    api = module()
    with pytest.raises(api.FacetValidationError):
        api.parse_country_facet(payload(count, countries))


def test_zero_records_are_a_complete_query_without_invented_country_rows():
    api = module()
    result = api.parse_country_facet(payload(0, ()))
    assert result["reported_count"] == 0
    assert result["country_counts"] == {}
    assert result["missing_country_count"] == 0


def test_nonempty_results_rejects_an_accidental_occurrence_download():
    api = module()
    data = json.loads(payload())
    data["results"] = [{"key": 123}]
    with pytest.raises(api.FacetValidationError, match="occurrence"):
        api.parse_country_facet(json.dumps(data).encode())


def test_species_basis_requires_both_windows_and_failed_queries_remain_explicit(tmp_path):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    source = api.load_original_manifest(path)
    queries = api.build_queries(source)
    chosen = [q for q in queries if q["species"] == "Chlosyne lacinia"]
    receipts = []
    for query in chosen:
        if query["basis_scope"] == "core" and query["window"] == "2018_2025":
            receipts.append({**query, "status": "FAILED", "complete": False, "error": "HTTP 503"})
        else:
            receipts.append({**query, "status": "COMPLETE", "complete": True,
                             **api.parse_country_facet(payload())})
    result, rows = api.assemble_outputs(source, queries, receipts)
    chlosyne = next(r for r in result["species_receipts"] if r["species"] == "Chlosyne lacinia")
    assert chlosyne["basis_scopes"]["all_accepted_bases"]["status"] == "COMPLETE"
    assert chlosyne["basis_scopes"]["core"]["status"] == "PARTIAL_FETCH_NOT_ELIGIBLE"
    assert chlosyne["basis_scopes"]["core"]["windows"]["2018_2025"]["status"] == "FAILED"
    assert chlosyne["status"] == "PARTIAL_FETCH_NOT_ELIGIBLE"
    assert {r["basis_scope"] for r in rows} == {"all_accepted_bases"}
    assert len(rows) == 4
    assert len(result["species_receipts"]) == 239
    assert len(result["query_receipts"]) == 924
    assert result["query_statuses"] == {"COMPLETE": 3, "FAILED": 1, "MISSING": 920}
    unresolved = [r for r in result["species_receipts"] if r["status"] == "TAXON_UNRESOLVED"]
    assert len(unresolved) == 8
    assert result["accepted_key_collisions"] == [{
        "taxon_key": 9804158, "species": ["Oeneis bore", "Oeneis polixenes"],
    }]


def test_checkpoint_reuse_requires_exact_query_and_raw_response_hash(tmp_path):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    query = api.build_queries(api.load_original_manifest(path))[0]
    raw = payload()
    receipt = api.save_complete_query(tmp_path, query, raw, attempts=[{"http_status": 200}])
    assert api.read_complete_checkpoint(tmp_path, query)["body_sha256"] == receipt["body_sha256"]
    (tmp_path / receipt["raw_path"]).write_bytes(raw[:-4])
    assert api.read_complete_checkpoint(tmp_path, query) is None
    api.save_complete_query(tmp_path, query, raw, attempts=[{"http_status": 200}])
    changed = {**query, "taxon_key": query["taxon_key"] + 1}
    assert api.read_complete_checkpoint(tmp_path, changed) is None


@pytest.mark.parametrize("malformed", [[], None, 7])
def test_nonobject_checkpoint_cannot_crash_either_resume_mode(tmp_path, monkeypatch, malformed):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    query = api.build_queries(api.load_original_manifest(path))[0]
    checkpoint = tmp_path / "checkpoints" / f"{query['query_id']}.json"
    checkpoint.parent.mkdir()
    checkpoint.write_text(json.dumps(malformed))
    assert api.read_complete_checkpoint(tmp_path, query) is None
    monkeypatch.setattr(api.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(
        stdout=payload() + b"\n200", stderr=b"", returncode=0,
    ))
    assert api.fetch_query(tmp_path, query, api.GlobalRateLimiter(1), tries=1)["status"] == "COMPLETE"
    checkpoint.write_text(json.dumps(malformed))
    assert api.main(["--original-manifest", str(path), "--output-dir", str(tmp_path), "--assemble-only"]) == 2


def test_shared_key_totals_are_explicitly_deduplicated(tmp_path):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    source = api.load_original_manifest(path)
    queries = api.build_queries(source)
    selected = [q for q in queries if q["species"] in ("Oeneis bore", "Oeneis polixenes")]
    receipts = [{**q, "status": "COMPLETE", "complete": True,
                 **api.parse_country_facet(payload())} for q in selected]
    result, _ = api.assemble_outputs(source, queries, receipts)
    totals = result["record_totals"]["all_accepted_bases"]["2010_2017"]
    assert totals["per_source_name"]["reported_count"] == 20
    assert totals["per_unique_accepted_key"]["reported_count"] == 10
    assert totals["per_unique_accepted_key"]["source_name_representatives"] == ["Oeneis bore"]
    assert result["accepted_key_count_discrepancies"] == []


def test_retry_after_accepts_seconds_or_http_date_from_final_response():
    api = module()
    assert api.parse_retry_after(b"HTTP/1.1 429 Too Many Requests\r\nRetry-After: 9\r\n\r\n") == 9
    now = datetime(2026, 10, 20, 12, 0, 0, tzinfo=timezone.utc)
    assert api.parse_retry_after(
        b"HTTP/1.1 200 Connection established\r\n\r\n"
        b"HTTP/2 429\r\nRetry-After: Tue, 20 Oct 2026 12:00:17 GMT\r\n\r\n", now=now,
    ) == 17
    assert api.parse_retry_after(b"HTTP/2 429\r\nRetry-After: invalid\r\n\r\n") is None


def test_shared_cooldown_delays_subsequent_request_starts(monkeypatch):
    api = module()
    clock = [100.0]
    monkeypatch.setattr(api.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(api.time, "sleep", lambda delay: clock.__setitem__(0, clock[0] + delay))
    limiter = api.GlobalRateLimiter(1)
    limiter.wait()
    limiter.cooldown(9)
    limiter.wait()
    assert clock[0] >= 109


def test_new_cooldown_reaches_workers_already_waiting_for_the_next_slot():
    api = module()
    clock = [100.0]
    sleep_entered, release_sleep, cooldown_returned = (threading.Event() for _ in range(3))
    started = []

    def controlled_sleep(delay):
        sleep_entered.set()
        release_sleep.wait(2)
        clock[0] += delay

    api.time = SimpleNamespace(monotonic=lambda: clock[0], sleep=controlled_sleep)
    limiter = api.GlobalRateLimiter(1)
    limiter.wait()

    def pending_start():
        limiter.wait()
        started.append(clock[0])

    def publish_cooldown():
        limiter.cooldown(9)
        cooldown_returned.set()

    pending = threading.Thread(target=pending_start)
    cooldown = threading.Thread(target=publish_cooldown)
    pending.start()
    try:
        assert sleep_entered.wait(2)
        cooldown.start()
        assert cooldown_returned.wait(1), "A sleeping worker prevents shared cooldown publication"
    finally:
        release_sleep.set()
        pending.join(2)
        if cooldown.ident is not None:
            cooldown.join(2)
    assert started and started[0] >= 109


@pytest.mark.parametrize("http_status,expected_attempts", [(429, 2), (400, 1)])
def test_transport_honors_retry_after_and_retries_only_transient_errors(tmp_path, monkeypatch, http_status, expected_attempts):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    query = api.build_queries(api.load_original_manifest(path))[0]
    clock, starts = [100.0], []
    monkeypatch.setattr(api.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(api.time, "sleep", lambda delay: clock.__setitem__(0, clock[0] + delay))

    def response(command, **kwargs):
        starts.append(clock[0])
        first = len(starts) == 1
        status = http_status if first else 200
        header_path = Path(command[command.index("--dump-header") + 1])
        header_path.write_bytes(f"HTTP/2 {status}\r\nRetry-After: 9\r\n\r\n".encode())
        body = b'"Throttled"' if first else payload()
        return SimpleNamespace(stdout=body + f"\n{status}".encode(), stderr=b"", returncode=0)

    monkeypatch.setattr(api.subprocess, "run", response)
    receipt = api.fetch_query(tmp_path, query, api.GlobalRateLimiter(1), tries=3)
    assert len(receipt["attempts"]) == expected_attempts
    if http_status == 429:
        assert receipt["status"] == "COMPLETE"
        assert starts[1] - starts[0] >= 9
        assert receipt["attempts"][0]["retry_after_seconds"] == 9
    else:
        assert receipt["status"] == "FAILED"


def test_resumed_failed_query_preserves_prior_failure_audit(tmp_path, monkeypatch):
    api = module()
    path = tmp_path / "original.json"
    path.write_text(json.dumps(original_manifest()))
    query = api.build_queries(api.load_original_manifest(path))[0]
    responses = iter([(400, b'"Invalid request"'), (200, payload())])

    def response(command, **kwargs):
        status, body = next(responses)
        return SimpleNamespace(stdout=body + f"\n{status}".encode(), stderr=b"", returncode=0)

    monkeypatch.setattr(api.subprocess, "run", response)
    failed = api.fetch_query(tmp_path, query, api.GlobalRateLimiter(1), tries=1)
    assert failed["status"] == "FAILED"
    resumed = api.fetch_query(tmp_path, query, api.GlobalRateLimiter(1), tries=1)
    assert resumed["status"] == "COMPLETE"
    assert [a["http_status"] for a in resumed["attempts"]] == [400, 200]
    archived = list((tmp_path / "failed_query_receipts").glob("*.json"))
    assert len(archived) == 1
    assert json.loads(archived[0].read_bytes()) == failed

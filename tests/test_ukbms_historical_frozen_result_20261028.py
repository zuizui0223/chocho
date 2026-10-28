"""Frozen source-only GBIF UKBMS fields do not establish real occupancy absences."""
from pathlib import Path
from hashlib import sha256
import json

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'provenance/archive/host_hindcast/results/ukbms_historical_event_source_20261028.json'
EXPECTED_SHA='e6972f8c3bbca1f19cc9f7e9f8ab1fea57354de6909391a4f21d9826ce766daf'


def test_fixed_real_source_and_nondetection_qualification():
    b=SOURCE.read_bytes()
    assert sha256(b).hexdigest()==EXPECTED_SHA
    r=json.loads(b)
    assert r['schema']=='chocho_UKBMS_historical_event_and_absence_GBIF_source_audit_v01'
    assert r['status']=='UKBMS_GBIF_DATASET_HISTORICAL_SOURCE_AUDITED'
    assert r['dataset_registration']=='1e266c3d-92ef-4d5a-8e4a-c04742c772c3'
    assert r['original_occurrence_identifiers_locations_photographs_exported'] is False
    assert r['inference']['explicit_same_site_survey_visits_with_non_detections_verified'] is False
    assert r['inference']['host_climate_hindcast_estimated'] is False


def test_old_and_late_total_source_presence_only():
    r=json.loads(SOURCE.read_text())
    a=r['windows']['1970_1990'];b=r['windows']['2010_2024']
    assert a['presence_records_in_GBIF_index']==595688
    assert b['presence_records_in_GBIF_index']==6293165
    assert a['explicit_GBIF_ABSENT_occurrence_records']==0
    assert b['explicit_GBIF_ABSENT_occurrence_records']==0
    for w in (a,b):
        sample=w['source_field_availability_from_fixed_first_12_rows']
        assert sample['returned_sample_rows']==12
        for f in ('eventID','samplingProtocol','individualCount'):
            assert sample['number_with_field_nonblank'][f]==0
        assert sample['number_with_field_nonblank']['locationID']==12
        assert sample['number_with_field_nonblank']['eventDate']==12


def test_source_does_not_export_original_occurrence_rows():
    raw=SOURCE.read_text()
    for x in ('"decimalLatitude":','"decimalLongitude":','"recordedBy":',
              '"gbifID":','"occurrenceId":','"results":','"site_id":'):
        assert x not in raw

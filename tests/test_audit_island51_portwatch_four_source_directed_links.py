import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_portwatch_four_source_directed_links import (
    PORTS,SOURCE_SAMPLE,ZENODO_VERIFIED,EXPECTED,BUDGET_COUNTS,
    fixed_sources,geo_km,main,MODELED_LAYER,MAX_FEATURE_EXAMPLES
)

def test_original_source_id_map_and_previous_2020_WPI_edge_list():
    if not (PORTS.exists() and SOURCE_SAMPLE.exists() and ZENODO_VERIFIED.exists()):
        pytest.skip("Original WPI files omitted from anonymous source-review snapshot")
    lookup,older=fixed_sources()
    assert len(lookup)==4 and len(older)==94
    assert set(lookup)==set(EXPECTED)
    assert len({x["provider_id"] for x in lookup.values()})==4
    assert len({x["island_id"] for x in lookup.values()})==4
    assert all(geo_km(x["wpi_lat"],x["wpi_lon"],
                      x["portwatch_registry_lat"],x["portwatch_registry_lon"])<=25
               for x in lookup.values())
    assert BUDGET_COUNTS==20 and MAX_FEATURE_EXAMPLES==12
    assert MODELED_LAYER.endswith("Spillover_Simulator_Maritime_Connections/FeatureServer/0")

def test_bounded_query_shape_and_stable_12_source_pairs(monkeypatch,tmp_path):
    if not (PORTS.exists() and SOURCE_SAMPLE.exists() and ZENODO_VERIFIED.exists()):
        pytest.skip("Original WPI and 2020 source files omitted intentionally from anonymous review")
    import audit_island51_portwatch_four_source_directed_links as module
    lookup,_=fixed_sources()
    byprovider={x["provider_id"]:x for x in lookup.values()}
    calls=[]
    def fake(query):
        assert query.get("f")=="json"
        assert query.get("returnGeometry") is None or query.get("returnGeometry")=="false"
        where=query["where"]
        calls.append(query)
        assert all(i in where for i in ["port"]) or False
        if query.get("returnCountOnly")=="true":
            return {"count":1},{"success":True,"attempts":[{"status":200}]}
        assert query["resultRecordCount"]=="1"
        assert len(calls)<=32
        names=where.replace("'","").split(" AND ")
        a=names[0].split("=")[1]
        b=names[1].split("=")[1]
        src=byprovider[a];dst=byprovider[b]
        attr={
          "from_id":a,"to_id":b,"from_name":src["wpi_name"],"to_name":dst["wpi_name"],
          "from_iso3":src["provider_ISO3"],"to_iso3":dst["provider_ISO3"],
          "from_lat":src["portwatch_registry_lat"],"from_lon":src["portwatch_registry_lon"],
          "to_lat":dst["portwatch_registry_lat"],"to_lon":dst["portwatch_registry_lon"],
          "distance":55.0,"transit_days":2.0
        }
        return {"features":[{"attributes":attr}]},{"success":True,"attempts":[{"status":200}]}
    monkeypatch.setattr(module,"fetch_public",fake)
    got=main(tmp_path)
    assert got["status"]=="ALL_12_FIXED_DIRECTED_MODEL_PAIR_QUERIES_AND_8_DEGREES_COMPLETED"
    assert len(got["paired_source_record_comparisons"])==12
    assert len(got["source_four_port_in_out_degree_count_queries"])==8
    assert got["n_pair_count_queries"]==12
    assert got["n_provider_degree_count_queries"]==8
    assert got["n_sample_feature_rows_requested"]==12
    assert len(calls)==32
    assert len({(x["original_WPI_origin"],x["original_WPI_destination"]) for x in
                got["paired_source_record_comparisons"]})==12
    assert all(x["PortWatch_count_query_complete"] for x in got["paired_source_record_comparisons"])
    assert all(x["portwatch_static_model_directed_record_count"]==1 for x in got["paired_source_record_comparisons"])

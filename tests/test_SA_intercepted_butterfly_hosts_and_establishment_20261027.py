"""Synthetic original sample evidence and metadata checks for actual butterfly transports."""
from __future__ import annotations
import sys
from pathlib import Path
from collections import Counter
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_SA_intercepted_butterfly_hosts_and_establishment_20261027 as m


def original_sample(insect,crop,year,commodity,regime="audit"):
    return {
      "Insecta: Lepidoptera":insect,
      "crop.scientificName":crop,
      "year":str(year),
      "crop.commodity":commodity,
      "sample.source":regime
    }


def source_metadata():
    return [
     {"scientificName":"Spodoptera frugiperda",
      "order":"Lepidoptera",
      "SA.occurenceStatus.historical":"absent",
      "SA.occurenceStatus.current":"present (2017)"},
     {"scientificName":"Pieris brassicae",
      "order":"Lepidoptera",
      "SA.occurenceStatus.historical":"absent",
      "SA.occurenceStatus.current":"absent"},
    ]


def original_rows():
    out=[
       original_sample("Pieris brassicae (Pieridae)","Brassica oleracea",2013,"Cuttings: leafy"),
       original_sample("Lampides boeticus (Lycaenidae)","Cajanus cajan",2011,"Grain"),
       original_sample("Spodoptera frugiperda (Noctuidae)","Zea mays",2014,"Seeds"),
       original_sample("Lepidoptera larva","Triticum aestivum",2012,"Grain"),
    ]
    out += [original_sample("Lepidoptera larva","Unknown",2014,"Fresh fruit") for _ in range(174)]
    assert len(out)==178
    return out


def test_original_taxon_parser_rejects_generic_larva_as_butterfly_species():
    assert m.parsed_species("Pieris brassicae (Pieridae)")=="Pieris brassicae"
    assert m.parsed_species("Lampides boeticus (Lycaenidae)")=="Lampides boeticus"
    assert m.parsed_species("Lepidoptera larva") is None
    assert m.parsed_species("Pyralidae larva") is None
    assert m.parsed_species("cf. Papilio demoleus") is None


def test_stage_flag_separate_from_taxon_name_and_shipping():
    assert m.explicit_stage_label("Pyralidae larva")=="EXPLICIT_LARVAL_STAGE"
    assert m.explicit_stage_label("Pieris brassicae (Pieridae)")=="LIFE_STAGE_UNSPECIFIED"
    assert m.explicit_stage_label("Lepidoptera egg")=="EXPLICIT_EGG_OR_OVIPOSITION"
    assert m.explicit_stage_label("Lepidoptera pupa")=="EXPLICIT_PUPAL_STAGE"


def test_crop_source_name_matching_never_equates_with_larval_feeding():
    r=m.interpret(original_rows(),source_metadata(),{
        "Pieris brassicae":{"Brassica oleracea"},
        "Lampides boeticus":{"Cajanus cajan"},
    })
    assert r["intercepted_exact_original_239_butterfly_names"]==[
        "Lampides boeticus","Pieris brassicae"]
    candidates=[
        x for x in r["intercepted_species_and_metadata_recent_establishment"]
        if x["in_original_GEB_239_butterflies"]
    ]
    assert len(candidates)==2
    assert sum(x["positive_original_plant_sample_rows"] for x in candidates)==2
    assert sum(x["original_crop_binomial_exact_to_original_known_HOSTS_plant_pair"]
               for x in candidates)==2
    known_status=next(q for q in candidates
                      if q["identified_lepidoptera_name"]=="Pieris brassicae")
    assert known_status["source_metadata_SA_species_status_available"]
    assert known_status["source_metadata_SA_historical_record_status"]=="absent"
    unknown_status=next(q for q in candidates
                        if q["identified_lepidoptera_name"]=="Lampides boeticus")
    assert not unknown_status["source_metadata_SA_species_status_available"]
    assert unknown_status["source_metadata_SA_historical_record_status"] is None
    assert r["original_SA_taxon_data_without_link_to_specific_established_descendants"]


def test_source_metadata_later_establishment_is_not_traced_lot():
    r=m.interpret(original_rows(),source_metadata(),{
        "Pieris brassicae":set(),"Lampides boeticus":set(),
    })
    assert r["moth_has_unambiguous_original_sample_positive"]
    x=next(z for z in r["intercepted_species_and_metadata_recent_establishment"]
           if z["identified_lepidoptera_name"]=="Spodoptera frugiperda")
    assert x["source_positive_samples_pre_2017"]==1
    assert x["source_earliest_positive_year"]==2014
    assert x["source_2017_year_in_current_status_is_metadata_not_specimen_date"]


def test_invasion_metadata_with_no_exact_sample_stays_unlinked():
    r=m.interpret([
        *original_rows()[:2],
        *[original_sample("Lepidoptera larva","Unknown",2014,"Grain")
          for _ in range(176)]
    ],source_metadata(),{
        "Pieris brassicae":{"Brassica oleracea"},
        "Lampides boeticus":{"Cajanus cajan"}
    })
    assert r["moth_has_unambiguous_original_sample_positive"] is False


def test_unverified_original_establishment_status_fails():
    meta=source_metadata()
    meta[0]["SA.occurenceStatus.current"]="absent"
    with pytest.raises(ValueError,match="status changed"):
        m.interpret(original_rows(),meta,{
            "Pieris brassicae":set(),"Lampides boeticus":set()
        })


def test_larval_stage_by_commodity_never_imputed_from_taxon_label():
    r=m.interpret(original_rows(),source_metadata(),{
        "Pieris brassicae":{"Brassica oleracea"},
        "Lampides boeticus":{"Cajanus cajan"},
    })
    by_type=r["all_positive_Lepidoptera_stages_by_original_commodity"]
    assert by_type["Grain"]["EXPLICIT_LARVAL_STAGE"]==1
    assert by_type["Grain"]["LIFE_STAGE_UNSPECIFIED"]==1
    assert by_type["Cuttings: leafy"]["LIFE_STAGE_UNSPECIFIED"]==1
    assert sum(sum(x.values()) for x in by_type.values())==178
    source=r["all_positive_Lepidoptera_stages_by_material_class"]
    assert source["SEED_GRAIN"]["EXPLICIT_LARVAL_STAGE"]==1
    assert "EXPLICIT_LARVAL_STAGE" not in source["LIVE_PROPAGATION"]

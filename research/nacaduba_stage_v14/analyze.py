"""Small independent, standard-library audit of source-dated history and larval stages."""
from pathlib import Path
import json

def audit(path=None):
    source=Path(path) if path is not None else Path(__file__).with_name("evidence.json")
    d=json.loads(source.read_text(encoding="utf-8"))
    years=d["historical_documentation"]
    a=years["acacia_singapore_first_official_year"]
    b=years["butterfly_singapore_earliest_retrospective_photo"]
    c=years["butterfly_public_species_checklist_year"]
    f=years["larval_host_rearing_documentation_year"]
    assert a<b<=c<f and years["documentation_lag_not_biological_lag"] is True
    stages=d["singapore_rearing_day_ranges"]
    assert len(stages)==7 and len(set(x["stage"] for x in stages))==7
    assert all(0<x["min"]<=x["max"] for x in stages)
    larvae=[x for x in stages if x["requires_floral_resource"]]
    assert {x["stage"] for x in larvae}=={"instar1","instar2","instar3","instar4"}
    egg=next(x for x in stages if x["stage"]=="egg")
    flow=[sum(x[k] for x in larvae) for k in ("min","max")]
    egg_through_feeding=[flow[0]+egg["min"],flow[1]+egg["max"]]
    full=[sum(x[k] for x in stages) for k in ("min","max")]
    hosts=d["host_evidence"]
    wild=sorted(x["plant"] for x in hosts if x["field_host"])
    nonfield=sorted(x["plant"] for x in hosts if not x["field_host"])
    assert len(wild)==3 and all(x["family"]=="Fabaceae" and x["mimosoid"] for x in hosts)
    src=d["v13_source_sensitivity"]
    base=src["original"]
    more=src["add_two_locally_documented_field_host_species_to_global_association"]
    assert base["native"]+base["introduced_only_added"]==base["contemporary"]
    assert more["native"]+more["introduced_only_added"]==more["contemporary"]
    return {
      "documentation_gaps_years":{"plant_to_retrospective_photo":b-a,"plant_to_checklist":c-a,"plant_to_larval_record":f-a},
      "estimated_biological_arrival_lag_years":None,
      "larval_floral_feeding_days":flow,
      "egg_to_larval_feeding_end_days":egg_through_feeding,
      "egg_to_adult_reared_days":full,
      "field_larval_host_plants":wild,
      "not_field_larval_hosts":nonfield,
      "source_hypothetical_resource_delta":{k:more[k]-base[k] for k in base},
      "field_larval_survival_rate":None,
      "flower_phenology_effect":None,
      "colonization_causal_effect":None
    }

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,ensure_ascii=False))

#!/usr/bin/env python3
"""Literature-grounded post-hoc non-equivalence of an introduced host for two Battus butterflies."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from collections import defaultdict

HOST_ID="2651580"
FOCAL="Battus philenor"
OTHER="Battus polydamas"

def rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def distributions(path):
    by=defaultdict(set)
    for r in rows(path):
        host=r["accepted_plant_name_id"].strip()
        area=r["area_code_l3"].strip()
        if host and area:by[host].add(area)
    return by

def main():
    p=argparse.ArgumentParser()
    for key in ("descriptors_csv","insect_host_csv","native_distribution_csv","contemporary_distribution_csv","output_json"):
        p.add_argument("--"+key.replace("_","-"),type=Path,required=True)
    args=p.parse_args()
    desc=rows(args.descriptors_csv)
    panel={r["species"] for r in desc if float(r["host_family_count"])>0 and int(r["host_wgsrpd3_unit_count"])>0}
    if len(panel)!=239:raise RuntimeError("fixed 239 species panel drift")
    by=defaultdict(set)
    for r in rows(args.insect_host_csv):
        by[r["insect_species"].strip()].add(r["accepted_plant_name_id"].strip())
    native=distributions(args.native_distribution_csv)
    contemporary=distributions(args.contemporary_distribution_csv)
    if FOCAL not in panel or OTHER not in panel:raise RuntimeError("focal butterflies missing")
    if HOST_ID not in by[FOCAL] or HOST_ID not in by[OTHER]:raise RuntimeError("both original host links not found")
    alien_areas=contemporary[HOST_ID]-native[HOST_ID]
    if len(native[HOST_ID])!=12 or len(contemporary[HOST_ID])!=60 or len(alien_areas)!=48:
        raise RuntimeError("frozen host plant geographic status drift")
    if "FLA" not in alien_areas:raise RuntimeError("Florida should be an alien range of Aristolochia littoralis")

    def footprint(sp,history,discount_intro=False):
        result=set()
        for host in by[sp]:
            allowed=history.get(host,set())
            if discount_intro and host==HOST_ID:
                allowed=native.get(HOST_ID,set())
            result.update(allowed)
        return result

    species={}
    for sp in (FOCAL,OTHER):
        n=footprint(sp,native)
        c=footprint(sp,contemporary)
        no_alien=footprint(sp,contemporary,discount_intro=True)
        extra=c-n
        dependent=(c-no_alien)-n
        species[sp]={
            "documented_host_plant_species":len(by[sp]),
            "native_known_resource_regions":len(n),
            "contemporary_known_resource_regions":len(c),
            "added_regions":len(extra),
            "added_regions_unique_to_alien_Aristolochia_littoralis":len(dependent),
            "fraction_added_regions_unique_to_alien_Aristolochia_littoralis":len(dependent)/len(extra) if extra else None,
            "dependency_region_codes":sorted(dependent),
            "native_resource_present_in_Florida":int("FLA" in n),
            "introduced_A_littoralis_present_in_Florida":int("FLA" in alien_areas)
        }

    co_users=sorted(s for s in panel if HOST_ID in by[s])
    pair_loss={}
    for other in co_users:
        if other==FOCAL:continue
        shared=by[FOCAL]&by[other]
        if HOST_ID not in shared:continue
        n=set().union(*(native[h] for h in shared))
        c=set().union(*(contemporary[h] for h in shared))
        post=set().union(*((native[h] if h==HOST_ID else contemporary[h]) for h in shared))
        newly_shared=c-n
        lost_new=newly_shared-(post-n)
        pair_loss[other]={
            "shared_known_host_species":len(shared),
            "native_exact_shared_host_regions":len(n),
            "contemporary_exact_shared_host_regions":len(c),
            "added_shared_host_regions":len(newly_shared),
            "added_shared_host_regions_dependent_on_B_philenor_Aristolochia_littoralis_link":len(lost_new),
            "lost_added_shared_host_region_codes":sorted(lost_new)
        }
    if species[FOCAL]["added_regions_unique_to_alien_Aristolochia_littoralis"]!=30:
        raise RuntimeError("expected 30 dependent regions for B.philenor")
    if species[OTHER]["added_regions_unique_to_alien_Aristolochia_littoralis"]!=36:
        raise RuntimeError("expected 36 dependent regions for B.polydamas")
    if sum(x["added_shared_host_regions_dependent_on_B_philenor_Aristolochia_littoralis_link"] for x in pair_loss.values())!=95:
        raise RuntimeError("expected 95 pair-region units dependent on B.philenor A.littoralis link")
    payload={
      "schema":"chocho_battus_exotic_host_fitness_sign_contrast_v0.1",
      "status":"LITERATURE_GROUNDED_POSTHOC_STRUCTURAL_SENSITIVITY",
      "plant":{"taxon":"Aristolochia littoralis","synonym":"Aristolochia elegans","wcvp_plant_id":HOST_ID,"native_WGSRPD3_regions":12,"contemporary_WGSRPD3_regions":60,"introduced_WGSRPD3_regions":48,"introduced_in_Florida":True},
      "frozen_butterfly_panel":239,
      "host_users_in_panel":co_users,
      "species":species,
      "pair_loss":pair_loss,
      "total_new_exact_shared_host_pair_region_units_removed_if_B_philenor_alien_link_unsuitable":sum(x["added_shared_host_regions_dependent_on_B_philenor_Aristolochia_littoralis_link"] for x in pair_loss.values()),
      "evidence":{
        "UF_Battus_polydamas_exotic_plant_feeding":"https://ask.ifas.ufl.edu/publication/IN219",
        "UF_Battus_philenor_exotic_plant_trap":"https://ask.ifas.ufl.edu/publication/IN1170",
        "UF_Aristolochia_introduced_Florida":"https://plant-directory.ifas.ufl.edu/plant-directory/aristolochia-littoralis/",
        "independent_2021_B_polydamas_eating_GloBI_original":"https://www.inaturalist.org/observations/91645851"
      },
      "interpretation":{
        "fitness":"B. polydamas larvae use the exotic host in Florida, whereas B. philenor females readily oviposit but larvae typically fail to survive according to UF/IFAS reviews. No survival is estimated across 48 global introduced regions.",
        "scope":"The 30 B. philenor potential regional opportunities and 95 exact shared-host pair-region units are model-dependent structural possibilities, NOT measured trap populations or demographic losses.",
        "novelty":"Differential feeding/survival of this exact plant for these two butterflies is already in published research; the added analysis diagnoses why a shared host edge is not equivalent to shared positive fitness opportunity."
      }
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"plant":payload["plant"],"species_summary":{k:{x:v for x,v in val.items() if x!="dependency_region_codes"} for k,val in species.items()},"pair_loss_summary":{k:{x:v for x,v in val.items() if x!="lost_added_shared_host_region_codes"} for k,val in pair_loss.items()}},indent=2))

if __name__=="__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,statistics
from pathlib import Path
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--overlap-json",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args()
    o=json.loads(a.overlap_json.read_text())
    rows=o["species"]
    outside=[r for r in rows if r["observed_outside_native"]>0]
    rescued=[r for r in outside if r["outside_native_units_explained_by_introduced_host_ranges"]>0]
    fr=[r["fraction_outside_native_explained_by_introduced"] for r in outside if r["fraction_outside_native_explained_by_introduced"] is not None]
    recovered=sum(r["outside_native_units_explained_by_introduced_host_ranges"] for r in rows)
    outside_units=sum(r["observed_outside_native"] for r in rows)
    p={
      "schema":"chocho_butterfly_occurrence_resource_validation_v0.1",
      "status":"INDEPENDENT_PANEL_OCCURRENCE_VALIDATION",
      "panel_species":len(rows),
      "species_with_occurrence_records":sum(r["occurrence_records"]>0 for r in rows),
      "species_with_observed_units_outside_native_host_envelope":len(outside),
      "species_with_at_least_one_outside_native_unit_recovered_by_introduced_hosts":len(rescued),
      "outside_native_species_x_units":outside_units,
      "outside_native_species_x_units_recovered_by_introduced_hosts":recovered,
      "fraction_outside_native_species_x_units_recovered":recovered/outside_units,
      "median_species_fraction_recovered_among_species_with_outside_native_units":statistics.median(fr),
      "species_with_all_outside_native_units_recovered":sum(v==1 for v in fr),
      "species_with_at_least_half_outside_native_units_recovered":sum(v>=0.5 for v in fr),
      "claim_boundary":"Presence overlap validates envelope relevance but does not prove larval use at each occurrence, causality of range expansion, or true absence in unrecovered units."
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(p,indent=2,sort_keys=True)+"\n")
    print(json.dumps(p,indent=2,sort_keys=True))
if __name__=="__main__":
    main()

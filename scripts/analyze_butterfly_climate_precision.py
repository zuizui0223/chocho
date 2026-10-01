#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--primary-result-json",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args()
    p=json.loads(a.primary_result_json.read_text(encoding="utf-8"))
    n=int(p["primary_test"]["n_species"])
    r=float(p["primary_test"]["observed_partial_spearman"])
    controls=1
    se=1/math.sqrt(n-controls-3)
    z=math.atanh(max(-0.999999,min(0.999999,r)))
    lo=math.tanh(z-1.959963984540054*se)
    hi=math.tanh(z+1.959963984540054*se)
    mde=math.tanh((1.6448536269514722+0.8416212335729143)*se)
    out={
      "schema":"chocho_butterfly_climate_precision_v0.1",
      "status":"POSTHOC_PRECISION_DIAGNOSTIC",
      "n_species":n,
      "controls":controls,
      "observed_partial_spearman":r,
      "approximate_fisher_95_percent_interval":[lo,hi],
      "approximate_one_sided_alpha_0_05_80_percent_power_detectable_absolute_partial_correlation":mde,
      "interpretation":"The independent panel rules out neither moderate negative nor modest positive associations. With n=24 and one control, effects around |partial rho|=0.50 are approximately the scale required for 80% power under Fisher-z correlation theory.",
      "claim_boundary":"This is an approximate Fisher-z precision calculation for a rank-based partial correlation, not a replacement for the preregistered permutation p-value."
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":main()

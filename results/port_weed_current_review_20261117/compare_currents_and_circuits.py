"""Compare dated model approaches with undirected circuit resistance on the same island pairs."""
import argparse
from pathlib import Path
import pandas as pd

def compare(source, output):
    def read(pattern):
        path, = source.rglob(pattern)
        return pd.read_csv(path)
    currents = read("*directed_source_nearshore_matrix.csv")
    resistance = read("*uniform_finvol_electrical_resistance.csv")
    flora = read("*strict_WPI2020_geographic_pairs.csv")
    geo = resistance.merge(flora, on=["source_region", "destination_region"], validate="one_to_one")
    keep = ["source_region", "destination_region", "source_verified_original_GSHHG_uniform_finite_volume_effective_resistance","centroid_geodesic_km"]
    geo = geo[keep].set_axis(["source","destination","resistance","km"], axis=1)
    both = pd.concat([geo,geo.rename(columns={"source":"destination","destination":"source"})])
    col = "NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"
    runs = currents.dropna(subset=[col]).merge(both,on=["source","destination"],validate="one_to_one")
    cases = runs[runs[col]>0].copy()
    reverse = runs.set_index(["destination","source"])[col]
    cases["reverse_approaches"] = [reverse.get((r.source,r.destination)) for r in cases.itertuples()]
    for factor in (1.5,2):
        for i,row in cases.iterrows():
            controls = runs[(runs.source==row.source)&(runs[col]==0)&runs.km.between(row.km/factor,row.km*factor)]
            cases.loc[i,f"controls_{factor}"] = len(controls)
            if len(controls):
                cases.loc[i,f"control_median_resistance_{factor}"] = controls.resistance.median()
    cases[["source","destination","km","resistance","reverse_approaches","controls_1.5","control_median_resistance_1.5","controls_2","control_median_resistance_2"]].to_csv(output,index=False)
    print("directions",len(runs),"positive",len(cases),"reverse",cases.reverse_approaches.value_counts(dropna=False).to_dict())
    for factor in (1.5,2):
        supported=cases[cases[f"controls_{factor}"]>=3]
        print(factor,len(supported),sum(supported.resistance<supported[f"control_median_resistance_{factor}"]),supported.source.nunique())

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("source",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    compare(args.source,args.output)

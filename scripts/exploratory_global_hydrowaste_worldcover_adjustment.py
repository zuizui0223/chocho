#!/usr/bin/env python3
"""2021 ESA WorldCover 10m built-up at all WPI1596 ports, urbanization sensitivity.

Exploratory cross-sectional controls only. Source WPI2019/HydroWASTE ~2022
and GBIF plant records 2001-2024 cannot infer historical nutrients or arrival.
No raw WorldCover COG tiles, HydroWASTE or protected plant links are vendored.
"""
from __future__ import annotations

import argparse,hashlib,json,math,sys
from collections import OrderedDict
from pathlib import Path
import numpy as np,pandas as pd
import rasterio
from rasterio.errors import RasterioIOError
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

SCRIPT_DIR=Path(__file__).resolve().parent
sys.path.insert(0,str(SCRIPT_DIR))
from exploratory_port_worldcover_urban import classify,tile_label,tile_url,CLASSES


class BoundedWorldCoverTiles:
    """LRU cache of remote COG readers; avoids thousands of open raster sockets."""
    def __init__(self,max_open=16):
        self.cache=OrderedDict()
        self.errors={}
        self.max_open=max_open
        self.opened=0
    def get(self,label):
        if label in self.cache:
            self.cache.move_to_end(label)
            return self.cache[label]
        if label in self.errors:return None
        try:
            src=rasterio.open(tile_url(label))
            if src.crs is None or src.crs.to_epsg()!=4326 or src.count!=1:
                src.close()
                raise ValueError(f"ESA tile CRS or bands mismatch {label}")
        except (RasterioIOError,ValueError) as exc:
            self.errors[label]=str(exc)[:300]
            return None
        self.cache[label]=src
        self.opened+=1
        if len(self.cache)>self.max_open:
            _,old=self.cache.popitem(last=False)
            old.close()
        return src
    def close(self):
        for src in self.cache.values():src.close()
        self.cache.clear()


def source_1596(path):
    raw=pd.read_csv(path,keep_default_na=False)
    if (len(raw),raw.port_id.nunique(),raw.country.nunique())!=(1596,1596,142):
        raise ValueError("WPI source-locked 1596/142 original countries no longer match")
    if (raw.in_within_country_comparison.sum(),
        raw[raw.in_within_country_comparison].country.nunique())!=(1368,52):
        raise ValueError("Original country fixed-effect study support changed")
    if raw.longitude.isna().any() or raw.latitude.isna().any():
        raise ValueError("Missing original WPI port positions")
    return raw


def sample_port_urban_features(source, *,radius_m=500,max_open=16):
    # Geographic grouping lowers remote COG open count and avoids a global GDAL
    # file-handle leak across roughly 1600 points.
    source=source.copy()
    source["_esa_tile"]=source.apply(
        lambda r:tile_label(float(r.longitude),float(r.latitude)),axis=1)
    sites=[]
    tiles=BoundedWorldCoverTiles(max_open=max_open)
    try:
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                          CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                          GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
                          GDAL_HTTP_CONNECTTIMEOUT="20",GDAL_HTTP_TIMEOUT="70"):
            for i,r in enumerate(source.sort_values(["_esa_tile","port_id"]).itertuples(),start=1):
                fake={"island_id":"NA","port_id":int(r.port_id),"country_original":r.country,
                      "port_sample_lon":float(r.longitude),"port_sample_lat":float(r.latitude)}
                q=classify(fake,"port",tiles,radius_m)
                sites.append({"port_id":int(r.port_id),
                    "esa_tile":r._esa_tile if hasattr(r,"_esa_tile") else tile_label(float(r.longitude),float(r.latitude)),
                    "quality":q["quality"],
                    "built_fraction_2021":q["built_fraction_among_land_pixels"],
                    "crop_fraction_2021":q["cropland_fraction_among_land_pixels"],
                    "tree_fraction_2021":q["tree_fraction_among_land_pixels"],
                    "grass_fraction_2021":q["grass_fraction_among_land_pixels"],
                    "terrestrial_pixels":q["count_terrestrial_nonwater_pixels"],
                    "circle_pixels":q["count_all_pixels_in_circle"],
                    "fraction_nonwater":q["nonwater_fraction_of_circle_pixels"],
                    "tile_missing":q["source_tiles_missing"]})
                if i%100==0 or i==len(source):
                    complete=sum(x["quality"]=="COMPLETE" for x in sites)
                    print("WORLDCOVER_WPI1596_PROGRESS",i,len(source),"complete",complete,
                          "COG_tiles_opened",tiles.opened,"missing_tiles",len(tiles.errors),flush=True)
    finally:
        tiles.close()
    if len(sites)!=1596 or len({q["port_id"] for q in sites})!=1596:
        raise ValueError("Missing source port land-cover data output rows")
    return pd.DataFrame(sites),dict(tiles_opened=tiles.opened,
              tile_failures=tiles.errors)


def country_model(source,*,urbancovariates):
    z=source[source.in_within_country_comparison].copy()
    z=z[z.quality=="COMPLETE"].copy()
    # Condition on urban data completeness but do not restrict response levels.
    eligible_countries=z.groupby("country").agg(n=("port_id","size"),
                                                  xvariation=("WWTP_count_25km","nunique"))
    eligible=eligible_countries[(eligible_countries.n>=5)&(eligible_countries.xvariation>=2)].index
    z=z[z.country.isin(eligible)].copy()
    if len(z)<500 or z.country.nunique()<18:
        return {"status":"INSUFFICIENT_URBAN_SOURCE_SUPPORT","n":len(z),
                "country_count":z.country.nunique()}
    z["abs_lat"]=z.latitude.abs()
    z["sinlon"]=np.sin(np.deg2rad(z.longitude))
    z["coslon"]=np.cos(np.deg2rad(z.longitude))
    z["ln_total"]=np.log1p(z.all_taxonomically_indexed)
    z["log_ratio"]=np.log((z.alien_candidates+.5)/(z.native_candidates+.5))
    terms=["ln_WWTP_25km","abs_lat","sinlon","coslon","ln_total","C(country)"]
    terms+=urbancovariates
    f=smf.ols("log_ratio ~ "+" + ".join(terms),data=z).fit(
        cov_type="cluster",cov_kwds={"groups":z.country})
    return {"status":"ESTIMATED_STATIC_2021_LANDCOVER_ASSOCIATION",
        "n_ports":len(z),"n_countries":z.country.nunique(),
        "beta":float(f.params["ln_WWTP_25km"]),
        "CI95_country_cluster":list(map(float,f.conf_int().loc["ln_WWTP_25km"])),
        "country_cluster_nominal_p":float(f.pvalues["ln_WWTP_25km"]),
        "included_urban_covariates":urbancovariates}


def run(port_path,*,radius_m=500):
    source=source_1596(port_path)
    pixels,tiles=sample_port_urban_features(source,radius_m=radius_m)
    combined=source.merge(pixels,on="port_id",validate="one_to_one")
    if len(combined)!=1596:raise ValueError("Port × urban source join wrong grain")
    good=combined[combined.quality=="COMPLETE"]
    rho=lambda a,b:float(spearmanr(good[a],good[b]).statistic) if len(good)>=4 else None
    summaries={
       "no_urban_on_identical_urban_complete_ports":country_model(combined,urbancovariates=[]),
       "adjust_2021_built_fraction":country_model(combined,urbancovariates=["built_fraction_2021"]),
       "adjust_2021_built_and_crop_fractions":country_model(combined,urbancovariates=["built_fraction_2021","crop_fraction_2021"]),
    }
    return {
       "schema":"chocho_HydroWASTE_WPI1596_ESA_WorldCover_2021_global_urban_confounding_v1",
       "status":"OBSERVATIONAL_2021_URBAN_PROXY_NOT_HISTORICAL_SEWAGE_NITROGEN",
       "source_1596_port_features_sha256":hashlib.sha256(Path(port_path).read_bytes()).hexdigest(),
       "source_WorldCover":"ESA WorldCover 2021 v200 native 10m GeoTIFF COG",
       "source_url":"https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/",
       "radius_m":radius_m,
       "source_cog_stats":tiles,
       "port_source_population":1596,
       "ports_with_complete_2021_landcover":len(good),
       "countries_with_2021_landcover":good.country.nunique(),
       "incomplete_port_ids":sorted(combined.loc[combined.quality!="COMPLETE","port_id"].astype(int).tolist()),
       "median_built_fraction_across_ports":float(good.built_fraction_2021.median()) if len(good) else None,
       "spearman_WWTP25_vs_urban_built_at_ports":rho("WWTP_count_25km","built_fraction_2021"),
       "spearman_WWTP25_vs_cropland_at_ports":rho("WWTP_count_25km","crop_fraction_2021"),
       "sensitivity_model_same_urban_eligible_population":summaries,
       "limitations":[
           "2021 WorldCover 500m terrestrial pixel fraction is not in-situ urbanization or 2001–2010 land use.",
           "GBIF host species queried within 5km and ESA built fraction within 500m: different footprints.",
           "The WorldCover sampling excludes sea and nodata from fractions, may be absent in tiny islands and polar regions.",
           "WWTP data represent treatment infrastructure and potential discharged water, not total dissolved nitrogen or soil enrichment.",
           "Urban landcover and port infrastructure can be strongly co-located; conditioning on present-day landcover can over-control or miss historical urbanization.",
           "Country FE, survey proxy and 2021 built fraction do not establish nitrogen-caused butterfly-host colonization.",
       ]},combined


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ports",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--table",type=Path,required=True)
    p.add_argument("--radius-m",type=int,default=500,choices=(250,500))
    a=p.parse_args()
    result,port=run(a.ports,radius_m=a.radius_m)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
    a.table.parent.mkdir(parents=True,exist_ok=True)
    port.to_csv(a.table,index=False)
    print("GLOBAL_WORLDCOVER_WWTP_URBAN_CONFOUNDING_GATES",json.dumps({
        "n_raster_valid":result["ports_with_complete_2021_landcover"],
        "n_countries_raster_valid":result["countries_with_2021_landcover"],
        "WWTP_urban_rho":result["spearman_WWTP25_vs_urban_built_at_ports"],
        "models":result["sensitivity_model_same_urban_eligible_population"]},sort_keys=True),flush=True)

if __name__=="__main__":
    main()

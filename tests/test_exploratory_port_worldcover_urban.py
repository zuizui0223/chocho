"""Offline testing of true ESA WorldCover native COG sampling, no web."""
import importlib.util, sys, numpy as np
from pathlib import Path
import rasterio
from rasterio.transform import from_origin

f=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_worldcover_urban.py"
sp=importlib.util.spec_from_file_location("worldcover_port",f)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_tile_names_account_for_negative_lower_left():
    assert m.tile_label(1.1,51.1)=="N51E000"
    assert m.tile_label(-1.1,51.1)=="N51W003"
    assert m.tile_label(-1.1,-52.1)=="S54W003"
    assert m.tile_label(139.8,35.7)=="N33E138"
    assert len(m.window_tiles(0.001,51.0,500))>=1

def test_land_cover_pixel_counts_and_land_only_denominator(tmp_path):
    fn=tmp_path/"mock_cog.tif"
    a=np.full((60,60),50,dtype=np.uint8)
    a[0:10]=80
    a[20:30]=40
    with rasterio.open(fn,"w",driver="GTiff",height=60,width=60,count=1,
        dtype="uint8",crs="EPSG:4326",
        transform=from_origin(-.003,.003,.0001,.0001),nodata=0) as out:
        out.write(a,1)
    with rasterio.open(fn) as ds:
        result=m.disk_stat(ds,0,0,250)
    assert result["circle_pixels"]>0
    assert result["counts"].get(50,0)>0
    assert sum(result["counts"].values())==result["circle_pixels"]

def test_exact_island_27_pair_source_quality_from_original_document():
    observations=[]
    for island in range(1,7):
        for port in range(1, 6 if island < 6 else 3):
            for side in ("port","nonport"):
                for period in ("early_2001_2010","later_2011_2024"):
                    observations.append({"island_id":str(island),"port_id":port,
                        "side":side,"period":period,"plant_records":30,
                        "N_scored_candidate_host_species":2,
                        "N_scored_other_plant_species":5})
    z={"site_period_facets":312,"predeclared_source_pairs":78,"observations":observations}
    assert len(m.focal_27_eligibility(z))==27

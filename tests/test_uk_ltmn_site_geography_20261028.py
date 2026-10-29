"""Pure site and source survey protocol info checks, no original station data."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_uk_ltmn_site_geography_20261028 import evaluate

def test_site_can_be_individually_crosswalked_without_dropping_ungc_source():
    src={"S1":[20,25,{1988,1989,1990},{2008,2009,2010}],
         "S2":[7,31,{1980},{2010,2011,2012}],
         "S3":[20,30,{1985,1986,1987},{2013,2014,2015}]}
    cross={"S1":(300000,230000,("grassland","arable")),
           "S2":(300100,230500,("","")),
           "S3":(100000,350000,("woodland",""))}
    x=evaluate(src,cross)
    assert x["repeated_any_1976_95_vs_2005_22"]==3
    assert x["strict_year_and_visit_matched_transects"]==2
    assert x["source_geocoded_transects_matching_strict_panel"]==2
    assert x["habitats_on_geocoded_strict_panel"]=={"arable":1,"grassland":1,"woodland":1}
    assert x["original_site_numbers_easting_northing_and_habitat_by_site_exported"] is False

def test_geography_missing_cannot_be_inferred():
    src={"S1":[20,25,{1988,1989,1990},{2008,2009,2010}]}
    x=evaluate(src,{})
    assert x["strict_year_and_visit_matched_transects"]==1
    assert x["source_geocoded_transects_matching_strict_panel"]==0

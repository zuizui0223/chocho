"""Independent BNG UK hectad source tests; no original locations or hosts."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_UK_host_hectad_coverage_20261029 import (
  easting_northing_hectad,exact_bsbi_accepted_files)

def test_known_UK_bng_grid_os_quadrant_not_adjacent_country():
    assert easting_northing_hectad(530000,180000)=="TQ38"
    assert easting_northing_hectad(325000,674000)=="NT27"
    assert easting_northing_hectad(408000,287000)=="SP08"
    assert easting_northing_hectad(-1,180000) is None

def test_only_exact_original_species_binomial_file_source():
    names=[
       {"taxonName":"Urtica dioica","ddbidLookup":"w1","rank":"species"},
       {"taxonName":"Urtica dioica","ddbidLookup":"w2","rank":"subspecies"},
       {"taxonName":"Salix caprea","ddbidLookup":"w3","rank":"species"},
    ]
    source={"hectadDistributions/w1.csv","hectadDistributions/w2.csv"}
    m,q=exact_bsbi_accepted_files(names,source)
    assert m=={"Urtica dioica":"hectadDistributions/w1.csv"}
    assert q["lookup_missing_exact_taxon_csv"]==1
    assert "Salix caprea" not in m

"""Non-network tests: name matching never guesses historical presence."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_original_hosts_name_match_20261029 import canonical

def test_accept_exact_taxon_names_and_reject_underspecified_genus():
    assert canonical(" Urtica   dioica ")==canonical("Urtica dioica")
    assert canonical("Urtica sp.")!=canonical("Urtica dioica")

def test_original_whitespace_not_synonym_or_cultivation():
    assert canonical("Rumex hydrolapathum")!=canonical("Rumex japonicus")
    assert canonical("Brassica oleracea")!=canonical("Brassica rapa")

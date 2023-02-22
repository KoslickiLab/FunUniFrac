import os
import sys
from scripts.reproducibility.generate_ko_hierarchy import KEGG_KO_Extraction

def test_get_ko00001_tree():
    ko_extractor = KEGG_KO_Extraction()
    substree_dict = ko_extractor.run('br:ko00001')
    assert len(substree_dict) == 25651
import os, sys
import data
import json
# add parent directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))


class _KeggDatabase:
    """Provides Kegg related data
    Normally refers to local files.
    Let us not instantiate this, but use for a singleton variable.
    """
    def __init__(self) -> None:
        self._brites = None
        self._brite_names = None
        self._koids_human_names = None
    
    @property
    def brites(self):
        if self._brites is not None:
            return self._brites
        else:
            loc = "kegg/brites.txt"
            with open(data.get_data_abspath(loc), 'r') as f:
                self._brites = [
                    str(l).strip() for l in f.readlines()    
                ]
            return self._brites
    
    @property
    def brite_names(self):
        if self._brite_names is not None:
            return self._brite_names
        else:
            loc = "kegg/brite_names.json"
            with open(data.get_data_abspath(loc), 'r') as f:
                dic = json.loads(" ".join(f.readlines()))
                dic = {k.replace("br:",str()): v for k, v in dic.items()}
                self._brite_names = dic
            return self._brite_names
    
    @property
    def koids_human_names(self):
        if self._koids_human_names is not None:
            return self._koids_human_names
        else:
            loc = "kegg/koids_human_names.txt"
            with open(data.get_data_abspath(loc), 'r') as f:
                dic = dict(
                    l.strip().split('\t') for l in f.readlines()
                )
                dic = {k.replace("ko:",str()): v for k, v in dic.items()}
                self._koids_human_names = dic
            return self._koids_human_names

# singleton is enough
instance = _KeggDatabase()

import os, glob
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL_DIR = "data"

def get_data_abspath(path, raise_if_not_found=True):
    abspath = os.path.join(ROOT_DIR, REL_DIR, path)
    if not os.path.exists(abspath) and raise_if_not_found:
        raise Exception(f"data not found: {path} does not exist in {REL_DIR} folder. Full: {abspath}")
    return abspath


def get_data_abspaths(pattern, raise_if_not_found=True):
    abspattern = os.path.join(ROOT_DIR, REL_DIR, pattern)
    abspaths = glob.glob(abspattern)
    if len(abspaths)==0 and raise_if_not_found:
        raise Exception(f"no data found for pattern: {pattern} in {REL_DIR} folder")
    return abspaths



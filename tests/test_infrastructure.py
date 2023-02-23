

def test_root():
    import os, sys
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    directories = os.listdir(ROOT_DIR)
    if set(['data', 'src', 'scripts']).issubset(directories):
        assert True
    else:    
        # directory structure may have been changed.
        assert False
    


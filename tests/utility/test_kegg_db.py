import src.utility.kegg_db as kegg_db


def test__kegg_db_loading():
    data = kegg_db.instance.brites
    data = kegg_db.instance.brite_names
    data = kegg_db.instance.koids_human_names

#!/usr/bin/env python3
"""
Definition
"""


import argparse
import maxi_tree as mt


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Update the MaxiTree object from the database")
    parser.add_argument("-url", metavar="<str>", required=True, help="The endpoint with the name of the database where the MaxiTree object is saved")
    parser.add_argument("-taxid", metavar="<str>", required=True, help="Taxonomy ID to update or to add")

    return parser.parse_args()


if __name__ == '__main__':
    PARAM = agrs_gestion()
    maxitree = mt.MaxiTree.from_database(PARAM.url)
    maxitree.insert(PARAM.taxid)
    try:
        os.mkdir("./treeDB_data/")
    except OSError:
        print("The directory treeDB_data exists, the created pickle file will be added to this folder")

    maxitree.dump("./treeDB_data/maxi_tree.p")

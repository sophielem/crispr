"""
Definition
"""

import argparse
import requests
import pickle
import datetime
import getpass
import pycouch.wrapper as wrapper
from ete3 import NCBITaxa

def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description = "Check if taxonomy information given by user are correct")

    parser.add_argument("-taxid", metavar="<str>",
                        help="The taxonomy ID", required=True)
    parser.add_argument("-gcf", metavar="<str>",
                        help="The GCF associated to the taxonomy ID", required=True)
    parser.add_argument("-out", metavar="<str>",
                        help="The path to the pickle file to fill the taxon database", required=True)
    parser.add_argument("-url", metavar="<str>",
                        help="End point of the taxon database", required=True)
    parser.add_argument("-dbName", metavar="<str>",
                        help="Name of the taxon database", required=True)
    return parser.parse_args()


def valid_taxid(taxid):
    """
    Check if the taxon id given by the user is in the NCBI taxonomy
    database
    """
    ncbi = NCBITaxa()
    try:
        ncbi.get_lineage(taxid)
        return taxid
    except Exception as err:
        print("Program terminated&The taxon id given ({}) \
is not in the NCBI taxonomy database !".format(taxid))
        sys.exit()


def check_taxid_exists(taxid, db_name):
    """
    Definition
    """
    req = wrapper.couchGetRequest(db_name + "/" + taxid)
    if not wrapper.docNotFound(req):
        print("Be careful&Taxon ID ({}) already exists : ".format(taxid), end=" ")
        return req
    return None


def check_gcf(gcf, list_gcf):
    """
    Definition
    """
    if gcf in list_gcf:
        if gcf == list_gcf[0]:
            print("The given GCF is the current GCF")
        else:
            print("The given GCF is older than the current GCF")
    print("Unknow GCF for this taxonomy ID")


if __name__ == '__main__':
    PARAM = args_gestion()
    wrapper.setServerUrl(PARAM.url)
    # Taxonomy ID
    valid_taxid(PARAM.taxid)
    TAXID_EXIST = check_taxid_exists(PARAM.taxid, PARAM.dbName)

    # GCF
    if TAXID_EXIST:
        check_gcf(PARAM.gcf, TAXID_EXIST["GCF"])
    LIST_GCF = [PARAM.gcf] + TAXID_EXIST["GCF"] if TAXID_EXIST else [PARAM.gcf]

    # Fasta file

    # Write info_file
    JSON_DATA = {PARAM.taxid: {"GCF": LIST_GCF, "Date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                               "User": getpass.getuser(), "current": PARAM.gcf}
                }
    pickle.dump(JSON_DATA, open(PARAM.out, "wb"), protocol=3)

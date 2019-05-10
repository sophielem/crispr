"""
Check if a taxonomy ID given exists in the NCBI database. Then, check if this taxonomy ID is
present in a Taxonomy database than you can request with the package pyCouch and check
if the GCF given for this taxonomy ID already exists in the Taxonomy database
"""

import argparse
import requests
import pickle
import datetime
import getpass
import pycouch.wrapper as couchdb
from ete3 import NCBITaxa


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description = "Check if taxonomy information given by user are correct")

    parser.add_argument("-taxid", metavar="<str>",
                        help="The taxonomy ID", required=True)
    parser.add_argument("-gcf", metavar="<str>",
                        help="The GCF associated to the taxonomy ID", required=True)
    parser.add_argument("-url", metavar="<str>",
                        help="End point of the taxon database", required=True)
    parser.add_argument("-dbName", metavar="<str>",
                        help="Name of the taxon database", required=True)
    parser.add_argument("-plasmid", type=str2bool,
                        help="The genome is a plasmid", required=True)
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
    Check if taxon ID exists in the database
    If yes, return the content else return None
    """
    req = couchdb.couchGetRequest(db_name + "/" + taxid)
    if not couchdb.docNotFound(req):
        msg = "Be careful&Taxon ID ({}) already exists : ".format(taxid)
        return req, msg
    return None, None


def check_gcf(gcf, list_gcf, msg):
    """
    Check if the given GCF is in the database and print a message
    Return True if the GCF is different from the current GCF
    """
    if gcf in list_gcf:
        if gcf == list_gcf[0]:
            print("Program terminated&The given GCF is the current GCF for the taxon ID. \
Please contact the support if you really really really want to change this genome")
            return False
        else:
            print("{} The given GCF is older than the current GCF".format(msg))
    print("{} Unknow GCF for this taxonomy ID".format(msg))
    return True


if __name__ == '__main__':
    PARAM = args_gestion()
    couchdb.setServerUrl(PARAM.url)
    if not couchdb.couchPing():
        print("Program terminated&Can't connect to the database with URL {}".format(PARAM.url))
        sys.exit()
    # Taxonomy ID
    if not PARAM.plasmid:
        valid_taxid(PARAM.taxid)
    TAXID_EXIST, MSG = check_taxid_exists(PARAM.taxid, PARAM.dbName)
    # GCF
    if TAXID_EXIST:
        check_gcf(PARAM.gcf, TAXID_EXIST["GCF"], MSG)

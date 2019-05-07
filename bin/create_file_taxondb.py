"""
Create pickle file to inset into Taxonomy database
It can create from a GCF and Taxonomy ID given or from the json file genome_ref_taxid.json file

OUTUPUT:
taxon_dt={"1234": {"GCF": ["GCF_111", "GCF_112", "GCF_113"], "date": "19-04-2019", "User": "Queen"},
          "2345": {"GCF": ["GCF_211", "GCF_212", "GCF_213"], "date": "19-04-2019", "User": "Queen"}}
"""

import os
import sys
import datetime
import json
import pickle
import argparse
import pycouch.wrapper as couchdb


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    command = sys.argv[1] if len(sys.argv) > 1 else None
    parser = argparse.ArgumentParser(description="Convert the genomre_ref_taxid.json to a pickle\
                                                  file ready to be insert into database")
    subparsers = parser.add_subparsers(help="commands")

    # Parser from SCRATCH
    scratch_parser = subparsers.add_parser("scratch", help="Create files to insert from\
                                                            the genomre_ref_taxid.json file")
    scratch_parser.add_argument("-file", metavar="<str>", required=True,
                                help="The json file to convert")
    scratch_parser.add_argument("-user", metavar="<str>", nargs="?", help="The user",
                                default="MMSB")
    scratch_parser.add_argument("-out", metavar="<str>", nargs="?",
                                help="The name of the outputfile", default="taxon_dt.p")

    # Parser for a SINGLE genome
    single_parser = subparsers.add_parser("single", help="Create file to insert")
    single_parser.add_argument("-user", metavar="<str>", help="The user", nargs="?", default="MMSB")
    single_parser.add_argument("-gcf", metavar="<str>", required=True, help="GCF id")
    single_parser.add_argument("-taxid", metavar="<str>", required=True, help="Taxonomy id")
    single_parser.add_argument("-url", metavar="<str>", required=True,
                               help="End pooint for taxon Database with the name of the database")
    single_parser.add_argument("-out", metavar="<str>", nargs="?",
                               help="The name of the outputfile", default="taxon_dt.p")

    return parser.parse_args(), command


def set_env(json_file=None):
    """
    Create a folder taxonDB_data and check if the json_file given exists
    """
    try:
        os.mkdir("./taxonDB_data/")
    except OSError:
        print("Be careful : The directory taxonDB_data exists")

    if json_file and not os.path.isfile(json_file):
        sys.exit("The file does not exist, check your path")


def init_taxondt(gcfs, user):
    """
    Initialize the dictionary of taxon and return it
    """
    tmp_dic = {}
    tmp_dic["GCF"] = gcfs
    tmp_dic["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    tmp_dic["user"] = user
    tmp_dic["current"] = gcfs[0]

    return tmp_dic


if __name__ == '__main__':
    PARAM, COMMAND = args_gestion()
    # Insert from SCRATCH
    if COMMAND == "scratch":
        set_env(PARAM.file)
        DIC_REF = json.load(open(PARAM.file, "r"))

        TAXON_DT = {}
        DUPLICATE = []
        for org in DIC_REF:
            taxid = DIC_REF[org][1]
            gcf = "_".join(DIC_REF[org][0].split("_")[:-1])
            list_gcf = [gcf] if taxid not in TAXON_DT.keys() else TAXON_DT[taxid]["GCF"] + [gcf]
            if len(list_gcf) > 1: DUPLICATE.append(taxid)
            TAXON_DT[taxid] = init_taxondt(list_gcf, PARAM.user)

        if DUPLICATE:
            with open("duplicate_taxon.log", "w") as filout:
                [filout.write("{}\n".format(tax)) for tax in DUPLICATE]
    # Insert for a SINGLE genome
    else:
        set_env()
        TAXON_DT = {}
        couchdb.setServerUrl(PARAM.url)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the taxon database with the URL : {}".format(PARAM.url))
            sys.exit()
        DOC = couchdb.couchGetDoc("", PARAM.taxid)
        LIST_GCF = PARAM.gcf + DOC["GCF"] if DOC else [PARAM.gcf]
        TAXON_DT[PARAM.taxid] = init_taxondt(LIST_GCF, PARAM.user)

    pickle.dump(TAXON_DT, open("./taxonDB_data/" + PARAM.out, "wb"), protocol=3)

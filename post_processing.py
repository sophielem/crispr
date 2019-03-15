#!/usr/bin/env python3
"""
Definition
"""

import sys
import os
import time
import argparse
import wordIntegerIndexing as decoding
import display_result as dspl
import pycouch.wrapper as couchDB


def valid_file(parser, filename):
    """
    Check if the file exists
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        parser.error("The file {} does not exist !".format(filename))
    # If all is good, return the filename
    return filename


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Post-processing results")
    parser.add_argument("-file", metavar="FILE", type=lambda x: valid_file(parser, x),
                        help="The path to the fasta file")
    parser.add_argument("-sl", metavar="<int>",
                        help="The length of the sgrna, excluding pam")
    parser.add_argument("-pam", metavar="<str>",
                        help="The pam motif",
                        required=True)
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-gni", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    return parser.parse_args()


def couchdb_search(sgrna_list, genomes_in):
    """
    Definition
    """
    couchDB.setServerUrl('http://localhost:1234')
    if couchDB.couchPing():
        dic_seq = {}
        start = time.time()
        results = couchDB.bulkRequestByKey(sgrna_list, "crispr_alpha")
        dspl.eprint("Time for the request : " + str(time.time() - start))
        # For each request
        for rslt in results["results"]:
            request_sgrna = rslt["id"]
            doc = rslt["docs"][0]
            # The request if ok
            if list(doc.keys())[0] == 'ok':
                dic_seq[request_sgrna] = {}
                # For each key in request result, so in json file
                for org_name in doc["ok"]:
                    if org_name not in ["_id", "_rev"] and org_name in genomes_in:
                        dic_seq[request_sgrna].update(doc["ok"][org_name])
        return dic_seq
    else:
        sys.exit("Error : can't connect to the database")


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gi.split("&")
    with open(PARAM.file, "r") as index_file:
        LIST_INDEX = index_file.read().splitlines()
    if LIST_INDEX:
        # Decoding of each index into sequence
        LIST_WORDS = [decoding.decode(rank, ["A", "T", "C", "G"], int(PARAM.sl) + int(len(PARAM.pam)))
                      for rank in LIST_INDEX]
        # Search coordinates for each sgrna in each organism
        DIC_SEQ = couchdb_search(LIST_WORDS, GENOMES_IN)
        # Display the result for the navigator
        dspl.display_hits(DIC_SEQ, GENOMES_IN, GENOMES_NOTIN,
                          PARAM.pam, int(PARAM.sl), ".")
    else:
        print("Program terminated & No hits remain")

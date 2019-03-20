#!/usr/bin/env python3
"""
Definition
"""

import sys
import os
import time
import argparse
import requests
import re
import wordIntegerIndexing as decoding
import display_result as dspl


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
    # parser.add_argument("-file", metavar="FILE", type=lambda x: valid_file(parser, x),
    #                     help="The path to the fasta file")
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
    parser.add_argument("-r", metavar="<str>",
                        help="The end point",
                        required=True)
    parser.add_argument("-c", metavar="<int>",
                        help="The length of the slice for the request",
                        required=True)
    parser.add_argument("-f", metavar="<str>",
                        help="The file index",
                        required=True)               
    parser.add_argument('--no-proxy', action='store_true')

    return parser.parse_args()


def couchdb_search(sgrna_list, genomes_in, end_point, len_slice, noPoxyBool):
    """
    Definition
    """
    reqFunc = requests
    if noPoxyBool:
        reqFunc = requests.Session()
        reqFunc.trust_env = False

    joker = 0
    results = {"request": {}}
    try :
        res = reqFunc.get(end_point + "/handshake")
        data = res.json()
        dspl.eprint("HANDSHAKE PACKET:\n" + str(data) + "\n")
    except Exception as e:
        dspl.eprint("Could not perform handshake, exiting")
        print("Program terminated&No handshake")
        sys.exit(1)

    while True:
        try:
            for i in range(0, len(sgrna_list), len_slice):
                dspl.eprint(i)
                request_sliced = {"keys" :sgrna_list[i : i + len_slice]}
                results["request"].update(reqFunc.post(end_point + "/bulk_request",
                                                        json=request_sliced).json()["request"])
        except Exception as e:
            dspl.eprint("Something wrong append, retrying time", str(joker))
            dspl.eprint("Error LOG is ", str(e))
            joker += 1
            if joker > 50:
                print("Program terminated&after 50 tries to access to the database")
                sys.exit(1)
            time.sleep(5)
            continue
        break

    dspl.eprint("motif-broker ans:\n", str(results))

    dic_seq = {}
    for sgrna in results["request"]:
        dic_seq[sgrna] ={}
        for org_name in results["request"][sgrna]:
            noBACKSLASH_org_name = org_name.replace('/','_')
            if noBACKSLASH_org_name in genomes_in:
                dic_seq[sgrna][noBACKSLASH_org_name] = results["request"][sgrna][org_name]

    #dspl.eprint("filtered word dict:\n" + str(dic_seq))
    return dic_seq


def parse_setcompare_out(fileIn):
    """
    Definition
    """
    with open(fileIn, "r") as filin:
        text = filin.readlines()
     #   dspl.eprint("--->" + str(text))

    nb_hits = re.search("[0-9]+", text[-2]).group()
    #dspl.eprint("nb hits: " + str(nb_hits))
    if int(nb_hits) == 0:
        return [], 0
    index_line = text[-1].split(",")
    index_line[-1] = index_line[-1].strip()
    index_line = [int(index) for index in index_line]
    #dspl.eprint("\nindex line: " + str(index_line))
    return index_line, nb_hits


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    # with open(PARAM.file, "r") as index_file:
    #     LIST_INDEX = index_file.read().splitlines()
    LIST_INDEX, NB_HITS = parse_setcompare_out(PARAM.f)
    if LIST_INDEX:

        dspl.eprint("LIST INDEX & NB HITS\n" + str(LIST_INDEX) +  "\n" + str(NB_HITS))
        # Decoding of each index into sequence
        LIST_WORDS = [decoding.decode(rank, ["A", "T", "C", "G"], int(PARAM.sl) + int(len(PARAM.pam)))
                      for rank in LIST_INDEX]
        dspl.eprint("LIST_WORDS\n" + str(LIST_WORDS) )
        # Search coordinates for each sgrna in each organism
        DIC_SEQ = couchdb_search(LIST_WORDS, GENOMES_IN, PARAM.r, int(PARAM.c), PARAM.no_proxy)
        dspl.eprint("DIC:: filtered word dict:\n" + str(DIC_SEQ))
        # Display the result for the navigator
        dspl.display_hits(DIC_SEQ, GENOMES_IN, GENOMES_NOTIN,
                          PARAM.pam, int(PARAM.sl), ".")

        print(','.join(GENOMES_NOTIN))
        print("TASK_KEY")
        print(len(DIC_SEQ))
        dspl.eprint("??? => " + str(NB_HITS) )

    else:
        print("Program terminated&No hits remain")

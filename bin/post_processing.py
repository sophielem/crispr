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
from collections import OrderedDict

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
    parser.add_argument("-nb_top", metavar="<int>",
                        help="The top hits to download",
                        default=1000)
    return parser.parse_args()


def couchdb_search(dic_hits, genomes_in, end_point, len_slice, noPoxyBool):
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
            sgrna_list = list(dic_hits.keys())
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


    for sgrna in results["request"]:
        dic_seq = {}
        for org_name in results["request"][sgrna]:
            noBACKSLASH_org_name = org_name.replace('/','_')
            if noBACKSLASH_org_name in genomes_in:
                dic_seq[noBACKSLASH_org_name] = results["request"][sgrna][org_name]
        dic_hits[sgrna].set_genomes_dict(dic_seq)

    return dic_hits


def check_find_database(dic_hits):
    """
    Check if all sgrna were found in database
    """
    for sgrna in dic_hits:
        if not dic_hits[sgrna].genomes_dict:
            dspl.eprint("Not find in database : {}".format(sgrna))


def parse_setcompare_out(fileIn, nb_top):
    """
    Definition
    """
    with open(fileIn, "r") as filin:
        text = filin.readlines()

    nb_hits = re.search("[0-9]+", text[-2]).group()
    if int(nb_hits) == 0:
        return {}, 0

    index_dic = OrderedDict()
    i = 0
    for rank_occ in text[-1].strip().split(","):
        if i == nb_top: break
        index_dic[int(rank_occ.split(":")[0])] = rank_occ.split(":")[1]
        i += 1
    # index_dic = {int(rank_occ.split(":")[0]) : rank_occ.split(":")[1] for rank_occ in text[-1].strip().split(",")}
    return index_dic, nb_hits


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    DIC_INDEX, NB_HITS = parse_setcompare_out(PARAM.f, int(PARAM.nb_top))
    if DIC_INDEX:

        dspl.eprint("DIC_INDEX & NB HITS\n" + str(DIC_INDEX.keys()) +  "\n" + str(NB_HITS))
        # Decoding of each index into sequence
        len_seq = int(PARAM.sl) + int(len(PARAM.pam))
        DIC_HITS = OrderedDict()
        for rank in DIC_INDEX:
            sequence = decoding.decode(rank, ["A", "T", "C", "G"], len_seq)
            DIC_HITS[sequence] = dspl.Hit(sequence, DIC_INDEX[rank])
        # DIC_HITS = {decoding.decode(rank, ["A", "T", "C", "G"], len_seq) : dspl.Hit(decoding.decode(rank, ["A", "T", "C", "G"], len_seq), DIC_INDEX[rank]) for rank in DIC_INDEX}

        # Search coordinates for each sgrna in each organism
        DIC_HITS = couchdb_search(DIC_HITS, GENOMES_IN, PARAM.r, int(PARAM.c), PARAM.no_proxy)

        # Display the result for the navigator
        dspl.display_hits(DIC_HITS, GENOMES_IN, GENOMES_NOTIN,
                          PARAM.pam, int(PARAM.sl), ".", int(PARAM.nb_top), True, list(DIC_HITS.keys()))

        print(','.join(GENOMES_NOTIN))
        print("TASK_KEY")
        print(len(DIC_HITS))
        dspl.eprint("??? => " + str(NB_HITS) )

    else:
        print("Program terminated&No hits remain")

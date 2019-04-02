#!/usr/bin/env python3
"""
Post-processing after the SetCompare script. Retrieve coordinates of sgRNA with requests into
database and display it for the webservice.
"""

import sys
import os
import time
import argparse
import re
from collections import OrderedDict
import requests
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
    parser.add_argument("-nb_top", metavar="<int>",
                        help="The top hits to download",
                        default=1000)
    return parser.parse_args()


def couchdb_search(sgrna_list, end_point, len_slice, no_poxy_bool):
    """
    Check if it can connect to the dabase
    Requests a list of sgrna and return the result
    """
    req_func = requests
    if no_poxy_bool:
        req_func = requests.Session()
        req_func.trust_env = False

    joker = 0
    results = {"request": {}}
    try:
        res = req_func.get(end_point + "/handshake")
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
                results["request"].update(req_func.post(end_point + "/bulk_request",
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

    # dspl.eprint("motif-broker ans:\n", str(results))
    return results


def treat_db_search_20(dic_hits, genomes_in, end_point, len_slice, no_poxy_bool):
    """
    Treat directly sequences with 20 length and return an update of a dictionary containing
    Hit objects
    """
    results = couchdb_search(list(dic_hits.keys()), end_point, len_slice, no_poxy_bool)
    for sgrna in results["request"]:
        dic_seq = {}
        for org_name in results["request"][sgrna]:
            nobackslash_org_name = org_name.replace('/', '_')
            if nobackslash_org_name in genomes_in:
                dic_seq[nobackslash_org_name] = results["request"][sgrna][org_name]
        dic_hits[sgrna].set_genomes_dict(dic_seq)

    return dic_hits


def treat_db_search_other(dic_hits, dic_index, genomes_in, end_point, len_slice, no_poxy_bool):
    """
    Find the result for word_20 associated to a word_15. Then merge results and update the
    Hit object
    """
    sgrna_list = [word for w15 in dic_index for word in dic_index[w15][1]]
    results = couchdb_search(sgrna_list, end_point, len_slice, no_poxy_bool)
    # Loop on word_15 to find their word_20 associated
    for sgrna in dic_hits:
        dic_seq = {}
        for word_20 in dic_index[sgrna][1]:
            for org_name in results["request"][word_20]:
                # To remove when the database is clean
                nobackslash_org_name = org_name.replace('/', '_')
                if nobackslash_org_name in genomes_in:
                    if nobackslash_org_name in dic_seq:
                        dic_seq[nobackslash_org_name].update(results["request"][word_20][org_name])
                    else:
                        dic_seq[nobackslash_org_name] = results["request"][word_20][org_name]
        dic_hits[sgrna].set_genomes_dict(dic_seq)
    return dic_hits


def check_find_database(dic_hits):
    """
    Check if all sgrna were found in database
    """
    for sgrna in dic_hits:
        if not dic_hits[sgrna].genomes_dict:
            dspl.eprint("Not find in database : {}".format(sgrna))


def parse_setcompare_other(rank_w15):
    """
    Return a tuple of occurence and list of word_20
    """
    rankw20_occ = rank_w15[1].split(",")
    return (rankw20_occ[0], rankw20_occ[1].strip().split(" "))


def parse_setcompare_out(output_c, nb_top, word_len):
    """
    Parse the output of setCompare C script and return a dictionary with the coding word
    and its occurence if the length is 20 or a tuple containing its occurence and
    a list of word_20 associated
    """
    with open(output_c, "r") as filin:
        for line in filin:
            regex_nb_hits = re.search("^# ([0-9]+)", line)
            if regex_nb_hits:
                nb_hits = int(regex_nb_hits.group(1))
                if nb_hits == 0:
                    print("Program terminated&No hits remain")
                    sys.exit()
                break

        index_dic = OrderedDict()
        i = 0
        for rank_occ in filin:
            if i == nb_top: break
            rank_splitted = rank_occ.split(":")
            index_dic[int(rank_splitted[0])] = rank_splitted[1] if word_len == 20 else parse_setcompare_other(rank_splitted)
            i += 1
    return index_dic, nb_hits


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    DIC_INDEX, NB_HITS = parse_setcompare_out(PARAM.f, int(PARAM.nb_top), int(PARAM.sl))

    LEN_SEQ = int(PARAM.sl) + int(len(PARAM.pam))
    DIC_HITS = OrderedDict()
    # Decoding of each index into sequence
    for rank in DIC_INDEX:
        sequence = decoding.decode(rank, ["A", "T", "C", "G"], LEN_SEQ)
        if int(PARAM.sl) == 20:
            DIC_HITS[sequence] = dspl.Hit(sequence, DIC_INDEX[rank])
        else:
            DIC_HITS[sequence] = dspl.Hit(sequence, DIC_INDEX[rank][0])

    # Search coordinates for each sgrna in each organism
    if int(PARAM.sl) == 20:
        DIC_HITS = treat_db_search_20(DIC_HITS, GENOMES_IN, PARAM.r, int(PARAM.c),
                                      PARAM.no_proxy)
    else:
        DIC_HITS = treat_db_search_other(DIC_HITS, DIC_INDEX, GENOMES_IN, PARAM.r,
                                         int(PARAM.c), PARAM.no_proxy)

    # Display the result for the navigator
    dspl.display_hits(DIC_HITS, GENOMES_IN, GENOMES_NOTIN,
                      PARAM.pam, int(PARAM.sl), ".", int(PARAM.nb_top),
                      True, list(DIC_HITS.keys()))

    print(','.join(GENOMES_NOTIN))
    print("TASK_KEY")
    print(len(DIC_HITS))

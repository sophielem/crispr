#!/usr/bin/env python3
"""
    Create files for node
"""
DEBUG = True
import os
import argparse
import pickle
import common_functions as cf
import allgenomes as ag


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Save nodes union and intersection")
    parser.add_argument("-leaves", metavar="<str>",
                        help="List of leaves for the node", required=True)
    parser.add_argument("-nodes", metavar="<str>",
                        help="List of nodes for the node", required=True)
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the reference genome folder")
    args = parser.parse_args()
    return args


def setup_application(param):
    """
    Split list of leaves and nodes
    """
    return param.leaves.split("&"), param.nodes.split("&")


if __name__ == '__main__':
    PARAM = args_gestion()
    LIST_LEAVES, LIST_NODES = setup_application(PARAM)
    if DEBUG: print(LIST_LEAVES); print(LIST_NODES)

    workdir = os.getcwd()
    fasta_path = workdir + "/reference_genomes/fasta"
    dict_organism_code = cf.read_json_dic(PARAM.rfg +
                                          '/genome_ref_taxid.json')
    if LIST_LEAVES != [""]:
            dic_leaves = ag.construction(fasta_path, "NGG", 20, LIST_LEAVES, [],
                                         dict_organism_code, workdir, PARAM.rfg)

    if LIST_NODES:
        dic_nodes = {}

    # Fusion des dictionnaires
    if dic_nodes and dic_leaves:
        keys_leaves = set(dic_leaves.keys())
        keys_nodes = set(dic_nodes.keys())
        common_keys = list(keys_leaves & keys_nodes)
        dic_fusion = {}
        for sgrna in common_keys:
            dic_fusion[sgrna] = {}
            dic_fusion[sgrna] = {dic_leaves[sgrna], dic_nodes[sgrna]}
    elif dic_nodes:
        dic_fusion = dic_nodes
    else:
        dic_fusion = dic_leaves

    if DEBUG: cf.eprint(str(len(dic_fusion)) + " hits remain after ")

    # Grand dictionnaire contenant le dic_fusion et une clef pour les metadata

    #
    # if not os.path.isdir(PARAM.rfg + "/node"):
    #     os.mkdir(PARAM.rfg + "/node")
    # pickle.dump(dic_inter, open(PARAM.rfg + "/node/trois.p", "wb"), protocol=3)

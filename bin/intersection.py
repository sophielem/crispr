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
    leaves = param.leaves.split("&")
    nodes = param.nodes.split("&")
    if leaves == [""]: leaves = []
    if nodes == [""]: nodes = []
    return leaves, nodes


def fusion_dic(dic_leaves, dic_nodes):
    """
    Fusion of sequences sgrna found in nodes and in leaves. Then return the
    results
    """
    if dic_nodes and dic_leaves:
        # Only keep common keys between these 2 dictionnaries
        keys_leaves = set(dic_leaves.keys())
        keys_nodes = set(dic_nodes.keys())
        common_keys = list(keys_leaves & keys_nodes)
        dic_fusion = {}
        # Fusion of dict for common sequences
        for sgrna in common_keys:
            dic_fusion[sgrna] = {}
            dic_fusion[sgrna] = {dic_leaves[sgrna], dic_nodes[sgrna]}
        return dic_fusion
    # No leaves to add to the intersection
    elif dic_nodes:
        return dic_nodes
    # No nodes to add to the intersection
    else:
        return dic_leaves


if __name__ == '__main__':
    PARAM = args_gestion()
    LIST_LEAVES, LIST_NODES = setup_application(PARAM)
    if DEBUG: print(LIST_LEAVES); print(LIST_NODES)

    WORKDIR = os.getcwd()
    FASTA_PATH = WORKDIR + "/reference_genomes/fasta"
    DICT_ORGANISM_CODE = cf.read_json_dic(PARAM.rfg +
                                          '/genome_ref_taxid.json')
    dic_leaves = {}
    dic_nodes = {}

    if LIST_LEAVES:
            dic_leaves = ag.construction(FASTA_PATH, "NGG", 20, LIST_LEAVES, [],
                                         DICT_ORGANISM_CODE, WORKDIR, PARAM.rfg)
    if LIST_NODES:
        pass

    dic_fusion = fusion_dic(dic_leaves, dic_nodes)

    if DEBUG: cf.eprint(str(len(dic_fusion)) + " hits remain after ")

    # Grand dictionnaire contenant le dic_fusion et une clef pour les metadata
    dic_inter = {}
    dic_inter["metadata"] = LIST_LEAVES + LIST_NODES
    dic_inter["data"] = dic_fusion
    name_node = "trois"
    if not os.path.isdir(PARAM.rfg + "/node"):
        os.mkdir(PARAM.rfg + "/node")
    pickle.dump(dic_inter, open(PARAM.rfg + "/node/" + name_node + ".p", "wb"), protocol=3)

#!/usr/bin/env python3
"""
    Create files for node
"""
DEBUG = True
import os
import argparse
import pickle
import time
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
                        help="The path to the reference genome folder",
                        required=True)
    args = parser.parse_args()
    return args


def setup_application(param):
    """
    Split list of leaves and nodes
    """
    leaves = param.leaves.split("&")
    nodes = param.nodes.split("&")
    if leaves == [""]: leaves = []
    if nodes == [""]:
        nodes = []
        nodes_path = []
    return leaves, nodes


def fusion_nodes(dic_ref, list_nodes):
    """
    Retrieve common keys of each dictionnary and keep them.
    A fusion dictionnary is created, using common keys.
    """
    # The list of path is empty
    if not list_nodes: return dic_ref

    node_path = list_nodes[0]
    dic_node = pickle.load(open(node_path, "rb"))["data"]
    # Retrieve keys for the new node and old nodes
    keys_node = set(dic_node.keys())
    keys_ref = set(dic_ref.keys())
    # Common keys between new and old nodes
    common_keys = list(keys_ref & keys_node)
    dic_inter_nodes = {}
    # Fusion of dict for common sequences
    for sgrna in common_keys:
        dic_inter_nodes[sgrna] = dict(dic_ref[sgrna], **dic_node[sgrna])
    return fusion_nodes(dic_inter_nodes, list_nodes[1: ])


def fusion_node_leaf(list_leaves, dic_node, fasta_path, dict_org_code, workdir,
                     uncompressed_file_path):
    """
    Uncompressed fasta files in a temporary directory, then sort genomes of
    leaves by size. Then, with all these included genomes, the principal
    search is launched
    """
    nb_file = 4
    cf.unzip_files(uncompressed_file_path, list_leaves, dict_org_code, workdir)
    sorted_genomes = cf.sort_genomes(list_leaves, fasta_path, dict_org_code,
                                     False)
    # All genomes are include for the intersection
    sorted_genomes = [(genome, "in") for genome in sorted_genomes]
    # Write sgRNA sequences common to all nodes in several files for bowtie
    list_fasta = cf.write_to_fasta_parallel(dic_node, nb_file, workdir)
    dic_seq = ag.principal_search(sorted_genomes, list_fasta, dict_org_code,
                                  dic_node, "NGG", 20, workdir, time.time())
    cf.delete_used_files(workdir)
    return dic_seq


def write_node_file(list_leaves, list_nodes, dic_fusion, path_to_write):
    """
    Intersection dictionnary, containing the list of nodes and leaves necessary
    to construct this node and the dictionnary containing sgRNA sequences and
    coordinates of each organism being part of the node
    """
    dic_inter = {}
    dic_inter["metadata"] = list_leaves + list_nodes
    dic_inter["data"] = dic_fusion
    name_node = "one"
    if not os.path.isdir(path_to_write):
        os.mkdir(path_to_write)
    pickle.dump(dic_inter, open(path_to_write + name_node + ".p", "wb"),
                                protocol=3)


def union_dic(list_elmt, dic_union, dict_org_code, is_leaf, param):
    """
    Definition
    """
    for elmt in list_elmt:
        # Define where search the pickle file
        if is_leaf:
            ref = dict_org_code[elmt][0]
            elmt = elmt.replace("/", "_")
            path_to_file = param.rfg + "/pickle/" + elmt + "." + ref + ".p"
            dic_elmt = pickle.load(open(path_to_file, "rb"))
        else:
            path_to_file = elmt
            dic_elmt = pickle.load(open(path_to_file, "rb"))["data"]

        # Retrieve all keys, so all sgRNA sequences
        keys_elmt = set(dic_elmt.keys())
        keys_ref = set(dic_union.keys())
        # Keep all distinct sgRNA sequences
        unique_keys = keys_ref.union(keys_elmt)
        for sgrna in unique_keys:
            # If sgRNA sequences is in both dict, merge
            if sgrna in keys_ref and sgrna in keys_elmt:
                dic_union[sgrna] = dict(dic_union[sgrna], **dic_elmt[sgrna])
            # If it is a new sequence, add it
            elif sgrna in keys_elmt:
                dic_union[sgrna] = dic_elmt[sgrna]

    return dic_union


def intersection(uncompressed_file_path, workdir, list_leaves, list_nodes, dict_org_code):
    """
    Definition
    """
    fasta_path = workdir + "/reference_genomes/fasta"
    nodes_path = [uncompressed_file_path + "/node/inter/" + node + ".p" for node in list_nodes]

    if list_leaves and not list_nodes:
        cf.eprint("Intersection of leaves")
         return ag.construction(fasta_path, "NGG", 20, list_leaves, [],
                                     dict_org_code, workdir, uncompressed_file_path)
    else:
        cf.eprint("Intersection of nodes")
        dic_node = fusion_nodes(pickle.load(open(nodes_path[0], "rb"))["data"],
                                nodes_path[1: ])
        if list_leaves:
            cf.eprint("Intersection of nodes and leaves")
            return fusion_node_leaf(list_leaves, DIC_NODE, fasta_path,
                                          dict_org_code, workdir, uncompressed_file_path)
        else:
             return dic_node


def union(arg):
    """
    Definition
    """
    pass


if __name__ == '__main__':
    PARAM = args_gestion()
    LIST_LEAVES, LIST_NODES = setup_application(PARAM)
    WORKDIR = os.getcwd()
    DICT_ORGANISM_CODE = cf.read_json_dic(PARAM.rfg +
                                          '/genome_ref_taxid.json')

    # DIC_FUSION = (PARAM.rfg, WORKDIR, LIST_LEAVES, LIST_NODES, DICT_ORGANISM_CODE)
    FASTA_PATH = WORKDIR + "/reference_genomes/fasta"
    NODES_PATH = [PARAM.rfg + "/node/inter/" + node + ".p" for node in LIST_NODES]

    if LIST_LEAVES and not LIST_NODES:
        cf.eprint("Intersection of leaves")
        dic_fusion = ag.construction(FASTA_PATH, "NGG", 20, LIST_LEAVES, [],
                                     DICT_ORGANISM_CODE, WORKDIR, PARAM.rfg)
    else:
        cf.eprint("Intersection of nodes")
        DIC_NODE = fusion_nodes(pickle.load(open(NODES_PATH[0], "rb"))["data"],
                                NODES_PATH[1: ])
        if LIST_LEAVES:
            cf.eprint("Intersection of nodes and leaves")
            dic_fusion = fusion_node_leaf(LIST_LEAVES, DIC_NODE, FASTA_PATH,
                                          DICT_ORGANISM_CODE, WORKDIR, PARAM.rfg)
        else:
            dic_fusion = DIC_NODE

    if DEBUG: cf.eprint(str(len(dic_fusion)) + " hits remain after ")

    write_node_file(LIST_LEAVES, LIST_NODES, dic_fusion, PARAM.rfg + "/node/inter/")

    NODES_PATH = [PARAM.rfg + "/node/union/" + node + ".p" for node in LIST_NODES]
    dic_union = {}

    if LIST_LEAVES:
        dic_union = union_dic(LIST_LEAVES, dic_union, DICT_ORGANISM_CODE, True, PARAM)
    if LIST_NODES:
        dic_union = union_dic(NODES_PATH, dic_union, DICT_ORGANISM_CODE, False, PARAM)


    if DEBUG: cf.eprint(str(len(dic_union)) + " hits remain after ")
    write_node_file(LIST_LEAVES, LIST_NODES, dic_union, PARAM.rfg + "/node/union/")

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
    else:
        nodes_path = [param.rfg + "/node/" + node + ".p" for node in nodes]
    return leaves, nodes, nodes_path


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
        dic_inter_nodes[sgrna] = {}
        dic_inter_nodes[sgrna] = {dic_ref[sgrna], dic_node[sgrna]}
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
    name_node = "three"
    if not os.path.isdir(path_to_write):
        os.mkdir(path_to_write)
    pickle.dump(dic_inter, open(path_to_write + name_node + ".p", "wb"),
                                protocol=3)


if __name__ == '__main__':
    PARAM = args_gestion()
    LIST_LEAVES, LIST_NODES, NODES_PATH = setup_application(PARAM)
    if DEBUG: print(LIST_LEAVES); print(LIST_NODES)

    WORKDIR = os.getcwd()
    FASTA_PATH = WORKDIR + "/reference_genomes/fasta"
    DICT_ORGANISM_CODE = cf.read_json_dic(PARAM.rfg +
                                          '/genome_ref_taxid.json')

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

    write_node_file(LIST_LEAVES, LIST_NODES, dic_fusion, PARAM.rfg + "/node/")

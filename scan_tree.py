#!/usr/bin/env python3
"""
Definition
"""


import sys
import json
import argparse
import pickle
from ete3 import NCBITaxa
# import intersection as inter

DEBUG = True

def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Allgenomes program.")
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-gni", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the reference genome folder")
    return parser.parse_args()


def setup_application(parameters, dict_organism_code):
    """
    Processes the arguments to be usable for research
    """
    org_selected = parameters.gi.split('&')
    org_excluded = parameters.gni.split('&')
    org_selected = [i for i in org_selected if i in dict_organism_code]
    org_excluded = [i for i in org_excluded if i in dict_organism_code]
    # List of taxon ID of included genome
    list_taxid_in = [dict_organism_code[sp][1] for sp in org_selected]
    # List of taxon ID of excluded genome
    list_taxid_notin = [dict_organism_code[sp][1] for sp in org_excluded]
    return list_taxid_in, list_taxid_notin


def create_topo_tree(dict_org, param):
    """
    Create a topology tree with the organism dictionary given in parameter
    """
    # List of taxon ID of all genomes
    all_spc = [DICT_ORG[sp][1] for sp in param.gi.split('&')] if DEBUG else [DICT_ORG[sp][1] for sp in DICT_ORG]
    # Create the entire topology tree
    ncbi = NCBITaxa()
    return ncbi.get_topology(all_spc)


def traversing_tree(subtree, list_taxid):
    """
    Check if children of a given node are all entered by user. If the child is
    a subnode, then the same thing is done for its children.
    If all leaves have been entered, return true else return false
    """
    # Retrieve name of all sisters leaves
    names_sister = [pop.name for pop in subtree.children if pop.is_leaf()]
    # Check all names leaves have been entered by user else return false, this node cannot be used
    if not set(names_sister).issubset(list_taxid): return False

    # Retrieve all sub-nodes
    names_subnodes = [pop for pop in subtree.children if not pop.is_leaf()]
    if names_subnodes:
        for subnode in names_subnodes:
            # Check all leaves of sub_nodes have been entered by user
            if not traversing_tree(subnode, list_taxid):
                return False
    # The node can be used
    return True


def issubnode(node_ref, nodes_name):
    """
    Definition
    """
    # Pas le nom, il faut l'objet noeud !!!!!!!
    children_name = [child.name for child in node_ref.children if not child.is_leaf()]
    # Enlever le noeud enfant de la liste et ajouter ce noeud
    to_remove = [child for child in children_name if child in nodes_name]
    nodes_name = [name for name in nodes_name if name not in to_remove]

    if not hasattr(node_ref.up, "name"):
        nodes_name.append(node_ref.name)
    elif node_ref.up.name not in nodes_name and node_ref.name not in nodes_name: # Ne pas ajouter le noeud dans la liste
        nodes_name.append(node_ref.name)
    return nodes_name


if __name__ == '__main__':
    PARAM = args_gestion()
    DICT_ORG = json.load(open(PARAM.rfg + '/genome_ref_taxid.json','r'))
    list_taxid_in, list_taxid_notin = setup_application(PARAM, DICT_ORG)
    tree_topo = create_topo_tree(DICT_ORG, PARAM)

    nodes_in = []
    nodes_notin = []
    for leaf_id in list_taxid_in:
        if traversing_tree(tree_topo.get_leaves_by_name(leaf_id)[0].up, list_taxid_in):
            if DEBUG: print(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)
            nodes_in = issubnode(tree_topo.get_leaves_by_name(leaf_id)[0].up, nodes_in)
            
            # dic_node_in = inter.inter_ndoes(pickle.load(open(nodes_in_name[0], "rb"))["data"], nodes_in_name[1:])
    for leaf_id in list_taxid_notin:
        if traversing_tree(tree_topo.get_leaves_by_name(leaf_id)[0].up, list_taxid_notin):
            if DEBUG: print(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)
            nodes_notin.append(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)
            dic_node_notin = inter.inter_ndoes(pickle.load(open(nodes_notin[0], "rb"))["data"], nodes_notin[1:])
# Exemple commande
# ./scan_tree.py -rfg ../reference_genomes_pickle -gi "Streptosporangium roseum DSM 43021 GCF_000024865.1&Catenulispora acidiphila DSM 44928 GCF_000024025.1&Streptomyces violaceusniger Tu 4113 GCF_000147815.2&Streptomyces bingchenggensis BCW-1 GCF_000092385.1&Archangium gephyra GCF_001027285.1" -gni ""

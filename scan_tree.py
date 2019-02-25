#!/usr/bin/env python3
"""
Definition
"""


import sys
import json
import argparse
from ete3 import NCBITaxa

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
    list_taxid_in = [DICT_ORG[sp][1] for sp in org_selected]
    # List of taxon ID of excluded genome
    list_taxid_notin = [DICT_ORG[sp][1] for sp in org_excluded]
    return list_taxid_in, list_taxid_notin


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


def create_topo_tree(dict_org, param):
    """
    Create a topology tree with the organism dictionary given in parameter 
    """
    # List of taxon ID of all genomes
    all_spc = [DICT_ORG[sp][1] for sp in param.gi.split('&')] if DEBUG else [DICT_ORG[sp][1] for sp in DICT_ORG]
    # Create the entire topology tree
    ncbi = NCBITaxa()
    return ncbi.get_topology(all_spc)


if __name__ == '__main__':
    PARAM = args_gestion()
    DICT_ORG = json.load(open(PARAM.rfg + '/genome_ref_taxid.json','r'))
    list_taxid_in, list_taxid_notin = setup_application(PARAM, DICT_ORG)
    tree_topo = create_topo_tree(DICT_ORG, PARAM)

    for leaf_id in list_taxid_in:
        if traversing_tree(tree_topo.get_leaves_by_name(leaf_id)[0].up, list_taxid_in):
            print(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)

    for leaf_id in list_taxid_notin:
        if traversing_tree(tree_topo.get_leaves_by_name(leaf_id)[0].up, list_taxid_notin):
            print(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)

# Exemple commande
# ./scan_tree.py -rfg ../reference_genomes_pickle -gi "Streptosporangium roseum DSM 43021 GCF_000024865.1&Catenulispora acidiphila DSM 44928 GCF_000024025.1&Streptomyces violaceusniger Tu 4113 GCF_000147815.2&Streptomyces bingchenggensis BCW-1 GCF_000092385.1&Archangium gephyra GCF_001027285.1" -gni ""

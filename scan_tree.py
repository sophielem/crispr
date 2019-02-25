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
    organisms_selected = parameters.gi.split('&')
    organisms_excluded = parameters.gni.split('&')
    organisms_selected = [i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded = [i for i in organisms_excluded if i in dict_organism_code]
    return organisms_selected, organisms_excluded


def traversing_tree(subtree, list_taxid):
    """
    Definition
    """
    # Récupérer le nom de toutes les feuilles soeurs
    names_sister = [pop.name for pop in subtree.children if pop.is_leaf()]
    # Plus élégant
    # Vérifier que tous les noms des feuilles ont été entrés par l'utilisateur
    # Sinon retourné faux, ce noeud ne peut pas être utilisé
    for sister in names_sister:
        if sister not in list_taxid:
            return False
    # Récupérer les sous-noeuds
    names_subnodes = [pop for pop in subtree.children if not pop.is_leaf()]

    if names_subnodes:
        for subnode in names_subnodes:
            # Vérifier que les feuilles des sous-noeuds ont été entrés par l'utilisateur
            if not traversing_tree(subnode, list_taxid):
                return False
    # Le noeud peut être utilisé
    return True


if __name__ == '__main__':
    PARAM = args_gestion()
    DICT_ORG = json.load(open(PARAM.rfg + '/genome_ref_taxid.json','r'))
    org_selected, org_excluded = setup_application(PARAM, DICT_ORG)
    # Liste des taxon ID des génomes inclus
    list_taxid_in = [DICT_ORG[sp][1] for sp in org_selected]
    # Liste des taxon ID des génomes exclus
    list_taxid_notin = [DICT_ORG[sp][1] for sp in org_excluded]
    # Liste des taxon ID de tous les génomes
    all_spc = [DICT_ORG[sp][1] for sp in org_selected] if DEBUG else [DICT_ORG[sp][1] for sp in DICT_ORG]
    # Construction de l'arbre entier
    ncbi = NCBITaxa()
    tree_topo = ncbi.get_topology(all_spc)

    for leaf_id in list_taxid_in:
        if traversing_tree(tree_topo.get_leaves_by_name(leaf_id)[0].up, list_taxid_in):
            print(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)

    for leaf_id in list_taxid_notin:
        if traversing_tree(tree_topo.get_leaves_by_name(leaf_id)[0].up, list_taxid_notin):
            print(tree_topo.get_leaves_by_name(leaf_id)[0].up.name)

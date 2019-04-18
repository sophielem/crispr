#!/usr/bin/env python3
"""
Add new genomes from a given node in an existing tree in json format
"""

import os
import sys
import json
import argparse
import shutil

def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Retrieve all children from a given node and create configure file to add them to the database")
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the reference genome folder",
                        required=True)
    parser.add_argument("-tree", metavar="<str>",
                        help="The path to the json tree",
                        required=True)
    parser.add_argument("-node", metavar="<str>",
                        help="The node to search after",
                        required=True)
    parser.add_argument("-dir", metavar="<str>",
                        help="The workdir where configure files will be created",
                        nargs="?", default="./")
    args = parser.parse_args()
    return args


def create_folder(path_folder):
    """
    Create folder to create configure files
    """
    if os.path.exists(path_folder):
        shutil.rmtree(path_folder)
    os.makedirs(path_folder)


def search_subtree(tree, value):
    """
    Return a subtree from a given node
    """
    # Node name is the searched node
    if tree["text"] == value: return tree
    # No children, so return None
    if "children" not in list(tree.keys()): return None
    # Traverse the tree via children
    for subref in tree["children"]:
        res = search_subtree(subref, value)
        if res != None:
            return res


def search_leaves(tree, list_child):
    """
    Return a list of all children from a given tree
    """
    # The node has not children, so it is a leave else traverse children
    if "children" not in list(tree.keys()):
        list_child.append(tree["text"])
        return list_child
    # Traverse the tree via children
    for subref in tree["children"]:
        list_child = search_leaves(subref, list_child)
    return list_child


def create_config_file(list_leaves, rfg, path_folder):
    """
    Copy fasta file and create config file for each genome to add
    """
    path_fasta = rfg + "/genome_fasta/"

    try:
        with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
            dic_ref = json.load(json_data)
    except:
        sys.exit("No genome_ref_taxid file in the path {}, can't copy fasta file".format(rfg))


    for leaf in list_leaves:
        try:
            ref = dic_ref[leaf][0]
            shutil.copy(path_fasta + ref + "_genomic.fna", path_folder)
            taxid = dic_ref[leaf][1]
            asm = ref.split("_")[-1]
            gcf = "_".join(ref.split("_")[ :-1])
            with open(path_folder + ref, "w") as config_file:
                config_file.write("GCF\tASM\tTaxon ID\n")
                config_file.write("{}\t{}\t{}".format(gcf, asm, taxid))
        except KeyError:
            print("The organism {} is not in the genome_ref_taxid file, can't find its reference".format(leaf))
        except FileNotFoundError:
            print("No fasta file found for {} at {}".format(leaf, path_fasta + ref + "_genomic.fna")
        except:
            print("Problem with the organism {}".format(leaf))


if __name__ == '__main__':
    PARAM = args_gestion()
    PATH_FOLDER = PARAM.dir + "/genomes_add/"
    create_folder(PATH_FOLDER)
    TREE = json.load(open(PARAM.tree, "r"))
    SUBTREE = search_subtree(TREE, PARAM.node)
    if not SUBTREE:
        sys.exit("No subtree found, check the name of the node")
    LIST_LEAVES = search_leaves(SUBTREE, [])
    with open("genomes_from_node.log", "w") as filout:
        for leaf in LIST_LEAVES:
            filout.write(leaf + "\n")
    create_config_file(LIST_LEAVES, PARAM.rfg, PATH_FOLDER)

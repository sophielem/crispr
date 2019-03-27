#!/usr/bin/env python3
"""
Add new genomes from a given node in an existing tree in json format
"""

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
    parser.add_argument("-node", metava="<str>",
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

    with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)

    for leaf in list_leaves:
        ref = dic_ref[leaf][0]
        shutil.copy(path_fasta + ref + "_genomic.fna", path_folder)

        gcf = re.search("GCF.*\_", ref).group()
        asm = ref.replace(gcf, "")
        gcf = gcf[:-1]
        ids = ".".join([gcf, asm, dic_ref[leaf][1]])
        with open(path_folder + "/" +  ref + "_genomic.conf", "r") as conf_file:
            conf_file.write(ids)


if __name__ == '__main__':
    PARAM = args_gestion()
    PATH_FOLDER = PARAM.dir + "/genomes_add"
    create_folder(PATH_FOLDER)

    SUBTREE = search_subtree(PARAM.tree, PARAM.node)
    LIST_LEAVES = search_leaves(SUBTREE, [])

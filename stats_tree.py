#!/usr/bin/env python3
"""
Count each prefix of a plain text of sequences
"""

import os
import copy
import pickle
import math
import argparse
import ete3
from tqdm import tqdm


def valid_file(parser, filename):
    """
    Check if the file exists else display an error
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        parser.error("The file {} does not exist !".format(filename))
    return filename


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Weighted tree")
    parser.add_argument("-file", metavar="FILE", type=lambda x: valid_file(parser, x),
                        help="The path to the file to load", required=True)
    parser.add_argument("-height", metavar="<int>", help="The depth for the tree",
                        required=True)
    parser.add_argument("-load", action="store_true",
                        help="Load the tree file if this flag is present")
    args = parser.parse_args()
    return args


def init_weighted_tree(height):
    """
    Initialize a weighted tree. If we describe depth tuple, then the first
    corresponds to the letter A, the second the letter C, the third the letter T
    and the last one the letter G. The first item of the list is the letter A
    and so on....
    """
    weighted_tree = ((([0], [0], [0], [0]),
                      ([0], [0], [0], [0]),
                      ([0], [0], [0], [0]),
                      ([0], [0], [0], [0])
                     ), (([0], [0], [0], [0]),
                         ([0], [0], [0], [0]),
                         ([0], [0], [0], [0]),
                         ([0], [0], [0], [0])
                        ), (([0], [0], [0], [0]),
                            ([0], [0], [0], [0]),
                            ([0], [0], [0], [0]),
                            ([0], [0], [0], [0])
                           ), (([0], [0], [0], [0]),
                               ([0], [0], [0], [0]),
                               ([0], [0], [0], [0]),
                               ([0], [0], [0], [0])
                              ))
    # Construct the tree. If we increase the depth by one, we quadruple
    # the number of leaves previous
    for i in range(3, height):
        weighted_tree = (copy.deepcopy(weighted_tree), copy.deepcopy(weighted_tree),
                         copy.deepcopy(weighted_tree), copy.deepcopy(weighted_tree))
    return weighted_tree


def construct_tuple_tree(filename, weighted_tree, height):
    """
    Read the file, and for each sequence, count increase the number of
    prefix. Then write this tuple in a pickle file. A pickle file must be
    created for each height
    """
    # Know the number of total word
    nb_word = 0
    aa_idx = {"A": 0, "C": 1, "T": 2, "G": 3}
    with open(filename, "r") as sequence_file:
        for word in tqdm(sequence_file):
            nb_word += 1
            # For each word, we retrieve the list of indices to load to reach
            # the leaf, then we add it 1
            index = [aa_idx[word[i]] for i in range(height)]
            tmp_tree = weighted_tree
            for i in index:
                tmp_tree = tmp_tree[i]
            tmp_tree[0] += 1
    # Write the tuple tree in a pickle file
    print("--- Save the tree in a pickle file ----")
    write_dic = {}
    write_dic["tree"] = weighted_tree
    write_dic["nb_word"] = nb_word
    pickle.dump(write_dic, open("weighted_tree_" + str(height) + ".p", "wb"))

    return nb_word, weighted_tree


def construct_nw_tree(cpt, aa_list, string_nw, seq, dic_weight, wgh_tree, nb_word, height):
    """
    Construct the newick string to be readable by ete3 package
    """
    # Leaf reached
    if cpt == height:
        for i, weight in enumerate(wgh_tree):
            letter = aa_list[i]
            dic_weight[seq + letter] = weight[0] / nb_word
            string_nw += seq + letter + ","
        cpt = 1
        return string_nw

    cpt += 1
    for i, subtree in enumerate(wgh_tree):
        string_nw += "("
        seq += aa_list[i]
        string_nw = construct_nw_tree(cpt, aa_list, string_nw, seq,
                                      dic_weight, subtree, nb_word, height)
        string_nw = string_nw[:-1]
        string_nw += "),"
        seq = seq[:-1]
    return string_nw


def add_weight_to_tree(tree_nw, dic_weight, height):
    """
    Add feature like weight to the Tree object to display it
    """
    factor = pow(10, height - 2)
    colors = ["red", "orange", "gold", "yellow", "greenyellow", "aquamarine",
              "mediumaquamarine", "seagreen", "green", "limegreen"]
    for sgrna in dic_weight:
        leaf = tree_nw.get_leaves_by_name(sgrna)[0]
        leaf.add_feature("weight", dic_weight[sgrna])
        leaf.add_face(ete3.TextFace(dic_weight[sgrna]), column=1,
                      position="branch-right")
        idx = math.trunc(dic_weight[sgrna] * factor)
        leaf.add_face(ete3.TextFace(leaf.name, fgcolor=colors[idx]), column=0)
    return tree_nw


if __name__ == '__main__':
    PARAM = args_gestion()
    HEIGHT = int(PARAM.height)

    print("--- Initialization of the tree ----")
    WEIGTHED_TREE = init_weighted_tree(HEIGHT)

    print("--- Load the file ----")
    if PARAM.load:
        # Load the tuple tree representing the weighted tree
        DIC_TMP = pickle.load(open(PARAM.file, "rb"))
        WEIGTHED_TREE = DIC_TMP["tree"]
        NB_WORD = DIC_TMP["nb_word"]
    else:
        NB_WORD, WEIGTHED_TREE = construct_tuple_tree(PARAM.file, WEIGTHED_TREE, HEIGHT)

    DIC_WEIGHT = {}
    AA_LIST = ["A", "C", "T", "G"]
    print("--- Write the newick string ----")
    NW_STRING = construct_nw_tree(1, AA_LIST, "(", "", DIC_WEIGHT,
                                  WEIGTHED_TREE, NB_WORD, HEIGHT)
    NW_STRING = NW_STRING[:-1]
    NW_STRING += ");"

    TREE_NW = ete3.Tree(NW_STRING)
    TREE_NW = add_weight_to_tree(TREE_NW, DIC_WEIGHT, HEIGHT)

    TS = ete3.TreeStyle()
    TS.show_leaf_name = False
    TS.min_leaf_separation = 20
    # TS.rotation = 90
    TREE_NW.render("weighted_tree_" + str(HEIGHT) + ".png", tree_style=TS)

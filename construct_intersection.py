#!/usr/bin/env python3
"""
Construct nodes intersection and union from scratch
"""


import json
from ete3 import NCBITaxa
import scan_tree as st


if __name__ == '__main__':
    DICT_ORG = json.load(open(sys.argv[1] + "/genome_ref_taxid.json", "r"))
    TREE_TOPO = st.create_topo_tree(DICT_ORG)

    list_terminal_leaves = []
    for node in TREE_TOPO:
        if node.is_leaf():
            flag = True
            for sister in node.get_sisters():
                if not sister.is_leaf():
                    flag = False
            if flag:
                ancestor = TREE_TOPO.get_common_ancestor(list_terminal_leaves + node)
                if ancestor.is_root():
                    list_terminal_leaves.append(node)

    for node in list_terminal_leaves:
        name_node = node.up.name
        list_leaves = [sister.name for sister in node.get_sisters()]
        list_leaves = "&".join(list_leaves)
        'python scan_tree -rfg ' + sys.argv[1] + ' -nodes "' + list_leaves + '" -leaves " " -name "' + name_node +'"'
        fname(node.up)
        node = node.up

def fname(elmt):
    name_node = node.name
    list_leaves = [child.name for child in elmt.get_children() if child.is_leaf()]
    list_leaves = "&".join(list_leaves) if list_leaves else " "

    list_nodes = [child.name for child in elmt.get_chidlren() if not child.is_leaf()]
    for node_sister in list_nodes:
        # Si le subnode n'a pas encore été créé, il va falloir le générer
        if not os.path.isfile(sys.argv[1] + "/node/inter/" + node_sister + ".p"):
            fname(node_sister)

    list_nodes = "&".join(list_nodes) if list_nodes else " "

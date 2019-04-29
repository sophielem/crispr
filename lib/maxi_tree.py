"""
Definition
"""

import os
import sys
import re
import json
import pickle
import argparse
import datetime
import pickle
import pycouch.wrapper as couchdb
from ete3 import Tree,NCBITaxa


class MaxiTree(object):
    """docstring for MaxiTree."""
    def __init__(self, dt_tree, all_spc):
        self.tree = dt_tree
        self.all_spc = all_spc

    @classmethod
    def from_database(cls, end_point):
        couchdb.setServerUrl(end_point)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon Tree database")
            sys.exit()
        maxi_tree = couchdb.couchGetDoc("", "maxi_tree")
        try:
            maxi_tree = json.loads(maxi_tree["tree"].replace("'", '"'))
        except json.JSONDecodeError:
            print("Program terminated&Problem with the maxi_tree in the database")
            sys.exit()

        tree_topo = Tree()
        cls.convert_json_to_tree(cls, tree_topo, maxi_tree)
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @staticmethod
    def convert_json_to_tree(cls, node, dic):
        if node.is_root() and node.name == '':
            node.name = dic["text"]
            new_node = node
        else:
            new_node = node.add_child(name=dic["text"])
        try:
            taxid = re.search(" : ([0-9]+)", dic["text"]).group(1)
            new_node.add_feature("taxon", taxid)
        except:
            pass
        if not "children" in dic.keys(): return None
        for i in dic["children"]:
            cls.convert_json_to_tree(cls, new_node, i)

    @classmethod
    def from_maxitree(cls, maxitree_file):
        if not os.path.isfile(maxitree_file): sys.exit("*** File does not exist ***")
        tree_topo = pickle.load(open(maxitree_file, "rb"))
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @classmethod
    def from_gen_ref(cls, gen_ref_file):
        if not os.path.isfile(gen_ref_file): sys.exit("*** File does not exist ***")
        dic_ref = json.load(open(gen_ref_file))
        all_taxid = [int(dic_ref[org][1]) for org in dic_ref]
        tree_topo = cls.construct_tree(cls, all_taxid)
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @staticmethod
    def construct_tree(cls, list_taxid):
        ncbi = NCBITaxa()
        tree_topo = ncbi.get_topology(list_taxid)
        # add the feature taxID
        for node in tree_topo:
            node.add_feature("taxon", node.name)
        couchdb.setServerUrl("http://127.0.0.1:5984/taxon_db")
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()
        # Change name by organism name
        for node in tree_topo.iter_descendants():
            node.name = rename_node(node.name, ncbi)
            if hasattr(node, "taxon"):
                gcf = couchdb.couchGetRequest(node.taxon)["current"]
                node.name = "{} {} : {}".format(node.name, gcf, node.taxon)
        tree_topo.name = ncbi.get_taxid_translator([int(tree_topo.name)])[int(tree_topo.name)]
        return tree_topo

    ### Write AND Representation of tree ###
    def __repr__(self):
        return "Number of leaves : {} from the root \"{}\"".format(len(self.all_spc), self.tree.name)

    def write_json_full(self):
        json_tree = self.get_json(True)
        with open("jsontree_full.json", "w") as filout:
            filout.write(json_tree)

    def write_json(self):
        json_tree = self.get_json(False)
        with open("jsontree.json", "w") as filout:
            filout.write(json_tree)

    def get_json(self, full):
        return str(self.__get_json__(self.tree, full))

    def __get_json__(self, node, full):
        if full:
            json = {"text": node.name}
        else:
            json = {"text": ":".join(node.name.split(":")[:-1])} if hasattr(node, "taxon") else {"text": node.name}
        if node.children:
            json["children"] = []
            for ch in node.children:
                json["children"].append(self.__get_json__(ch, full))
        return json

    def dump(self, filout):
        json_tree = self.get_json(True)
        dic_tree = {}
        dic_tree["maxi_tree"] = {}
        dic_tree["maxi_tree"]["tree"] = json_tree
        dic_tree["maxi_tree"]["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        pickle.dump(dic_tree, open(filout, "wb"), protocol=3)

    ### Modify Tree ###
    def is_member(self, taxid):
        return taxid in self.all_spc

    def insert(self, taxid):
        # Check if can connect to the database
        couchdb.setServerUrl("http://127.0.0.1:5984/taxon_db")
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()

        # Check if it is just a change of GCF
        if self.is_member(taxid):
            node = self.tree.search_nodes(taxon=str(taxid))[0]
            gcf_old = re.search("GCF_[0-9]+\.?[0-9]*", node.name)[0]
            gcf_new = couchdb.couchGetRequest(str(taxid))["current"]
            if gcf_old == gcf_new:
                print("Program terminated&Error, the new GCF is the same than the old GCF")
                sys.exit()
            node.name = node.name.replace(gcf_old, gcf_new)
            return True

        # Construct the topology tree with taxonID
        ncbi = NCBITaxa()
        tree_topo = ncbi.get_topology(self.all_spc + [taxid])
        for node in tree_topo:
            node.add_feature("taxon", node.name)
        for node in tree_topo.iter_descendants():
            node.name = rename_node(node.name, ncbi)
            if hasattr(node, "taxon"):
                if node.taxon == str(taxid):
                    try:
                        gcf = couchdb.couchGetRequest(node.taxon)["current"]
                        node.name = "{} {} : {}".format(node.name, gcf, node.taxon)
                    except:
                        print("Program terminated&Problem with the taxon {} in the Taxon database \
(WiP : insert a new member in the Tree)".format(taxid))
                        sys.exit()
                else:
                    node.name = self.tree.search_nodes(taxon=node.taxon)[0].name
        tree_topo.name = ncbi.get_taxid_translator([int(tree_topo.name)])[int(tree_topo.name)]

        # Rename the new node
        # try:
        #     ncbi = NCBITaxa()
        #     tree_topo = ncbi.get_topology(self.all_spc + [taxid])
        #     new_node = tree_topo.search_nodes(name=taxid)[0]
        #     # Rename the new node
        #     new_node.add_feature("taxon", new_node.name)
        #     new_node.name = rename_node(taxid, ncbi)
        #     # gcf = couchdb.couchGetRequest(taxid)["current"]
        #     gcf = "GCF_1234"
        #     new_node.name = "{} {} : {}".format(new_node.name, gcf, new_node.taxon)
        # except ValueError:
        #     print("Program terminated&This taxon does not exist in the NCBI database")
        #     sys.exit()
        # except Exception as e:
        #     print("Program terminated&Problem with the taxon {} in the Taxon database".format(taxid))
        #     print(e)
        #     sys.exit()
        #
        # self.__set_tree__(new_node, ncbi)
        self.tree = tree_topo
        self.all_spc = self.all_spc + [taxid]
        return True

    # def __set_tree__(self, new_node, ncbi):
    #     parent = new_node.up
    #     if not parent:
    #
    #     # The parent is a parent of another leaf, so exists in the original tree
    #     res = self.tree.search_nodes(name=ncbi.get_taxid_translator([int(parent.name)])[int(parent.name)])
    #
    #     if not res:
    #         # The parent is a leaf from the original tree
    #         res = self.tree.search_nodes(taxon=parent.name)
    #         # Remove the taxid and the GCF from the node because now it becomes a parent and not a leaf
    #         if res:
    #             remove = re.search("GCF_[0-9]+\.?[0-9]* : [0-9]+", res[0].name)[0]
    #             res[0].name = res[0].name.replace(remove, "").strip()
    #             res[0].del_feature("taxon")
    #
    #     if res:
    #          res[0].add_child(new_node)
    #     else:
    #         # The parent does not exist in the original tree
    #         parent.name = rename_node(parent.name, ncbi)
    #         self.__set_tree__(parent, ncbi)


def rename_node(taxid, ncbi):
    """
    Definition
    """
    node_name = ncbi.get_taxid_translator([int(taxid)])[int(taxid)]
    node_name = node_name.replace("'", '')
    node_name = node_name.replace("/", '_')
    return node_name


def args_gestion():
    """
    Definition
    """
    parser = argparse.ArgumentParser(description="Construct the MaxiTree and save it in the database")
    parser.add_argument("-file", metavar="<str>", required=True, help="The genome_ref_taxid.json file")
    return parser.parse_args()

if __name__ == '__main__':
    PARAM = args_gestion()
    TREE = MaxiTree.from_gen_ref(PARAM.file)
    JSON_TREE = TREE.get_json(True)
    DIC_TREE = {}
    DIC_TREE["maxi_tree"] = {}
    DIC_TREE["maxi_tree"]["tree"] = JSON_TREE
    DIC_TREE["maxi_tree"]["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    pickle.dump(DIC_TREE, open("jsontree_full.p", "wb"), protocol=3)

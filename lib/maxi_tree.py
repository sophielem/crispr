import json
import os
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
    def from_maxitree(cls, maxitree_file):
        if not os.path.isfile(maxitree_file): sys.exit("*** File does not exist ***")
        tree_topo = pickle.load(open(maxitree_file, "rb"))
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @classmethod
    def from_gen_ref(cls, gen_ref_file):
        if not os.path.isfile(gen_ref_file): sys.exit("*** File does not exist ***")
        dic_ref = json.load(open(gen_ref_file))
        list_taxid = [int(dic_ref[org][1]) for org in dic_ref]
        tree_topo = cls.construct_tree(cls, list_taxid)
        return cls(tree_topo, list_taxid)

    @staticmethod
    def construct_tree(cls, list_taxid):
        ncbi = NCBITaxa()
        tree_topo = ncbi.get_topology(list_taxid)
        # add the feature taxID
        for node in tree_topo:
            node.add_feature("taxon", node.name)
        # Change name by organism name
        couchdb.setServerUrl("http://127.0.0.1:5984/taxon_db")
        for node in tree_topo.iter_descendants():
            node.name = ncbi.get_taxid_translator([int(node.name)])[int(node.name)]
            node.name = node.name.replace("'", '')
            node.name = node.name.replace("/", '_')
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
        """
        Convert the Tree object to Json object
        """
        if full:
            json = {"text": node.name}
        else:
            json = {"text": ":".join(node.name.split(":")[:-1])} if hasattr(node, "taxon") else {"text": node.name}
        if node.children:
            json["children"] = []
            for ch in node.children:
                json["children"].append(self.__get_json__(ch, full))
        return json

    def dump(self, filin):
        pickle.dump(self.tree, open(filin, "wb"), protocol=3)

    ### Modify Tree ###
    def is_member(self, taxid):
        return taxid in self.all_spc

    def insert(self, taxid):
        if self.is_member(taxid):
            print("Taxon already into database")
            return False
        try:
            ncbi = NCBITaxa()
            ncbi.get_lineage(taxid)
        except Exception as e:
            print("This taxon does not exist in the NCBI database")
            return False

        self.all_spc = self.all_spc + [taxid]
        self.__set_tree__()
        return True

    def __set_tree__(self):
        self.tree = self.construct_tree(self, self.all_spc)

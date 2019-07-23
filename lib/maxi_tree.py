"""
Class MaxiTree object. This tree contains: a tree object from the package ete3 with the name
of the organism, the current GCF and its taxonomy Id if it's present into the CRISPR database. This taxon ID is
referenced in the name separate by a ':' and in a node's feature 'taxon'.
The name of a node is composed like this:
    name (without ' and /) GCF_ID : taxon_ID

Several constructor:
    1. from the database
    2. from the genomre_ref_taxid.json file
    3. from a MaxiTree pickle file
    4. from the taxonomy database

Several methods are implemented to insert a new node, to write a json file of the tree
full (with taxonID) or not, to get the json string full or not...

The main function allow to create this Tree from the genomre_ref_taxid.json file and to create a
pickle MaxiTree file to insert into the taxon_tree_db
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
        # Check if it can connect to the database
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon Tree database")
            sys.exit()
        # Try to get the MaxiTree object by the key maxi_tree
        maxi_tree = couchdb.couchGetDoc("", "maxi_tree")
        try:
            maxi_tree = json.loads(maxi_tree["tree"].replace("'", '"'))
        except json.JSONDecodeError:
            print("Program terminated&Problem with the maxi_tree in the database")
            sys.exit()
        # Convert the json string into a Tree object
        tree_topo = Tree()
        cls.convert_json_to_tree(cls, tree_topo, maxi_tree)
        # List all taxon ID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @staticmethod
    def convert_json_to_tree(cls, node, dic):
        # Give the name to the root
        if node.is_root() and node.name == '':
            node.name = dic["text"]
            new_node = node
        else:
        # Give the name to the node
            new_node = node.add_child(name=dic["text"])
        try:
        # Try to retrieve the taxon ID if it exists
            taxid = re.search(" : ([0-9]+)", dic["text"]).group(1)
            new_node.add_feature("taxon", taxid)
        except:
            pass
        # The current node is a leaf, so it is the end of the process for this subtree
        if not "children" in dic.keys(): return None
        # For each child, traverse it
        for i in dic["children"]:
            cls.convert_json_to_tree(cls, new_node, i)

    @classmethod
    def from_maxitree(cls, maxitree_file):
        # Check if the pickle file exists
        if not os.path.isfile(maxitree_file): sys.exit("*** File does not exist ***")
        try:
            tree_topo = pickle.load(open(maxitree_file, "rb"))
        except:
            sys.exit("Problem to open the pickle file {}".format(maxitree_file))
        # Generate the list of taxonID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @classmethod
    def from_gen_ref(cls, gen_ref_file, end_point):
        # Check if the genome_ref_taxid.json file exists
        if not os.path.isfile(gen_ref_file): sys.exit("*** File does not exist ***")
        try:
            dic_ref = json.load(open(gen_ref_file))
            all_taxid = [int(dic_ref[org][1]) for org in dic_ref]
        except:
            sys.exit("Problem with the format of the file {}".format(gen_ref_file))
        # Construct the Tree object
        # Add the taxonID of plasmid
        tree_topo = cls.construct_tree(cls, all_taxid + [36549], end_point)
        # Generate the list of TaxonID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @staticmethod
    def construct_tree(cls, list_taxid, end_point):
        # Check if it can connect to the database
        couchdb.setServerUrl(end_point)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()

        ncbi = NCBITaxa()
        # Create the topology Tree with taxon ID
        tree_topo = ncbi.get_topology(list_taxid)
        # Add the feature taxon for leaves
        for node in tree_topo.iter_descendants():
            if int(node.name) in list_taxid:
                if node.children:
                    node.add_child(name=node.name)
                else:
                    node.add_feature("taxon", node.name)

            node.name = rename_node(node.name, ncbi)
            if hasattr(node, "taxon"):
                # Add the GCF and taxon ID to the name
                gcf = couchdb.couchGetRequest(node.taxon)["current"]
                node.name = "{} {} : {}".format(node.name, gcf, node.taxon)
        # Name for the root
        tree_topo.name = ncbi.get_taxid_translator([int(tree_topo.name)])[int(tree_topo.name)]
        return tree_topo

    @classmethod
    def from_taxon_database(cls, end_point):
        couchdb.setServerUrl(end_point)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()
        # Retrieve all taxon ID and plasmid name
        all_entries = couchdb.couchGetRequest("_all_docs")
        list_taxon = [row['key'] for row in all_entries['rows']]
        # Only keep taxon id and cast into int
        list_taxond_id = [int(i)  for i in list_taxon if re.match("^[0-9]*$", i)]
        # Only keep plasmid name
        list_plasmids = [i  for i in list_taxon if not re.match("^[0-9]*$", i)]
        # Construct the Tree object and Add the taxonID of plasmid
        tree_topo = cls.construct_tree(cls, list_taxond_id + [36549], end_point)
        # Insert plasmid under node plasmid
        for plasmid in list_plasmids:
            cls.insert_plasmid(cls, plasmid, end_point, tree_topo)
        # Generate the list of TaxonID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

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
            # Remove the taxon ID if it is present
            json = {"text": ":".join(node.name.split(":")[:-1])} if hasattr(node, "taxon") else {"text": node.name}
        # If node has children, create a list of it and traverse it, else do not create attribute
        # children and return the json string
        if node.children:
            json["children"] = []
            for ch in node.children:
                json["children"].append(self.__get_json__(ch, full))
        return json

    def dump(self, filout):
        # Get the Tree in json string with taxonID
        json_tree = self.get_json(True)
        # Create a dictionary with the MaxiTree json and the date
        dic_tree = {}
        dic_tree["maxi_tree"] = {}
        dic_tree["maxi_tree"]["tree"] = json_tree
        dic_tree["maxi_tree"]["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        pickle.dump(dic_tree, open(filout, "wb"), protocol=3)

    ### Modify Tree ###
    def is_member(self, taxid):
        return taxid in self.all_spc

    def insert(self, taxid, end_point):
        # Check if can connect to the database
        couchdb.setServerUrl(end_point)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()

        # If the taxonID is already present, it is just a change of the GCF
        if self.is_member(taxid):
            # search the node in the original tree
            node = self.tree.search_nodes(taxon=str(taxid))[0]
            # Retrieve the GCF and compare it to the GCF given by the user
            gcf_old = re.search("GCF_[0-9]+\.?[0-9]*", node.name)[0]
            gcf_new = couchdb.couchGetRequest(str(taxid))["current"]
            if gcf_old == gcf_new:
                print("Program terminated&Error, the new GCF is the same than the old GCF")
                sys.exit()
            # Replace GCF
            node.name = node.name.replace(gcf_old, gcf_new)
            return True

        # Construct the topology tree with taxonID
        ncbi = NCBITaxa()
        tree_topo = ncbi.get_topology(self.all_spc + [taxid])
        # Rename node
        for node in tree_topo.iter_descendants():
        # Retrieve name from the original tree, faster than create it, just replace the GCF id
        # by the one given by the user which is already is the taxon_db
            if int(node.name) in self.all_spc:
                if node.children:
                    node.add_child(name=node.name)
                    node.name = rename_node(node.name, ncbi)
                else:
                    node.add_feature("taxon", node.name)
                    node.name = self.tree.search_nodes(taxon=node.name)[0].name
            elif node.name == str(taxid):
                if node.children:
                    node.add_child(name=node.name)
                    node.name = rename_node(node.name, ncbi)
                else:
                     try:
                         node.add_feature("taxon", node.name)
                         node.name = rename_node(node.name, ncbi)
                         gcf = couchdb.couchGetRequest(node.taxon)["current"]
                         node.name = "{} {} : {}".format(node.name, gcf, node.taxon)
                     except:
                         print("Program terminated&Problem with the taxon {} in the Taxon database \
    (WiP : insert a new member in the Tree)".format(taxid))
                         sys.exit()
            else:
                node.name = rename_node(node.name, ncbi)

        # Search for the node plasmid and replace it in the new tree by detaching the old plasmid
        # node without children
        plasmid_origin = self.tree.search_nodes(taxon="36549")[0]
        plasmid_new = tree_topo.search_nodes(taxon="36549")[0]
        plasmid_new.detach()
        try:
            tree_topo.add_child(plasmid_origin)
        except:
            pass

        # Name for the root
        tree_topo.name = ncbi.get_taxid_translator([int(tree_topo.name)])[int(tree_topo.name)]
        self.tree = tree_topo
        self.all_spc = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return True

    def insert_plasmid(self, name, end_point, tree=None):
        # Check if can connect to the database
        couchdb.setServerUrl(end_point)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()
        gcf = couchdb.couchGetRequest(name)["current"]
        name = "{} {}".format(name, gcf)
        plasmid_node = tree.search_nodes(taxon="36549")[0] if tree else self.tree.search_nodes(taxon="36549")[0]
        plasmid_node.add_child(name=name)
        new_plasmid = plasmid_node.search_nodes(name=name)[0]
        new_plasmid.add_feature("plasmid", name)
        return True

        ### WiP : keep the original Tree and just insert necessary nodes ###
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
        # self.tree = tree_topo
        # self.all_spc = self.all_spc + [taxid]
        # return True

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
    Rename node by searching the organism name from its taxonomy ID and
    replacing single quote by space and slash by underscore
    """
    node_name = ncbi.get_taxid_translator([int(taxid)])[int(taxid)]
    node_name = node_name.replace("'", '')
    node_name = node_name.replace("/", '_')
    return node_name


def args_gestion():
    """
    Takes and treat arguments given
    """
    parser = argparse.ArgumentParser(description="Construct the MaxiTree and save it in the database")
    parser.add_argument("-file", metavar="<str>", required=False, help="The genome_ref_taxid.json file")
    parser.add_argument("-url", metavar="<str>", required=False, help="Whole end_point for taxon_database")
    return parser.parse_args()

if __name__ == '__main__':
    PARAM = args_gestion()
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        sys.exit("Give a file or the end point for the taxon database")

    TREE = MaxiTree.from_gen_ref(PARAM.file) if PARAM.file else MaxiTree.from_taxon_database(PARAM.url)
    JSON_TREE = TREE.get_json(True)
    DIC_TREE = {}
    DIC_TREE["maxi_tree"] = {}
    DIC_TREE["maxi_tree"]["tree"] = JSON_TREE
    DIC_TREE["maxi_tree"]["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    pickle.dump(DIC_TREE, open("jsontree_full.p", "wb"), protocol=3)

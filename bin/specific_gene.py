#!/usr/bin/env python3
"""
Definition
"""

import sys
import argparse
import pickle
from collections import OrderedDict
import wordIntegerIndexing as decoding
import display_result as dspl
import post_processing as pp
from parse_blast import BlastReport, BlastHit


class ResumeSeq(object):
    """docstring for ResumeSeq."""
    def __init__(self, dic_org=None):
        self.proportion = 0
        self.dic_org = dic_org

    def __repr__(self):
        return "Proportion : {}\tTotal occurences : {}".format(self.proportion, self.total_occ())

    def cal_proportion(self, nb_total):
        """
        Calcul proportion of organism present in the dictionary
        """
        self.proportion = len(self.dic_org) / int(nb_total)
        return self.proportion

    def total_occ(self, org_name=None):
        """
        Calcul the total number of occurences for the sequence or for a given organism
        for a sequence
        """
        occ = 0
        if org_name:
            for coord_seq in self.dic_org[org_name]:
                occ += coord_seq.total_occ()
        else:
            for org_name in self.dic_org:
                occ += self.total_occ(org_name)
        return occ

    def update(self, dict_org):
        """
        Update the dictionary
        """
        if dict_org:
            self.dic_org.update(dict_org)

    def org_names(self):
        """
        Return a list of organism name which have this sequence
        """
        return list(self.dic_org.keys())

    def write(self, genomes_in):
        """
        Return a string of coordinates for each reference for each organism
        given an order by the list
        """
        to_write = ""
        for gen_in in genomes_in:
            if gen_in in self.org_names():
                for coord_seq in self[gen_in]:
                    to_write += coord_seq.ref + ':' + ','.join(coord_seq.list_coord) + ';'
            else:
                to_write += "None ;"
        to_write += self.total_occ()
        return to_write.strip(";")

    def list_ref(self, org_name):
        """
        Return a list of dictionaries containing coordinates of the sequence for
        a reference organism
        """
        list_ref = [{"ref": coord_seq.ref, "coords": coord_seq.list_coord} for coord_seq in self.dic_org[org_name]]
        return list_ref

    def list_occ(self):
        """
        Return a list of dicionaries containing a list_ref for all references for
        each organism
        """
        list_occurences = [{'org': genome, 'all_ref': self.list_ref(genome)} for genome in self.org_names()]
        return list_occurences



class CoordSeq(object):
    """docstring for CoordSeq."""
    def __init__(self, coords=None, ref=None):
        self.list_coord = coords
        self.ref = ref

    def __repr__(self):
        return "Total occurences : {}\tFor {} coordinates".format(self.total_occ(), len(self))

    def __len__(self):
        return len(self.list_coord)

    def update(self, coord):
        """
        Update the coordinate list
        """
        if coord not in self.list_coord:
            self.list_coord.append(coord)
            return True
        return False

    def total_occ(self):
        """
        The number of coordinates for the reference
        """
        return len(self.list_coord)


### METHODS FOR THE ARGUMENTS ###
def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(__file__, description="Specific gene program.")
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-gni", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    parser.add_argument("-sl", metavar="<int>",
                        help="The length of the sgrna, excluding pam",
                        nargs="?", const=20)
    parser.add_argument("-pam", metavar="<str>",
                        help="The pam motif",
                        nargs="?", const="NGG")
    parser.add_argument("-f", metavar="<str>",
                        help="The file index",
                        required=True)
    parser.add_argument("-blast", metavar="<str>",
                        help="The blast output file",
                        required=True)
    parser.add_argument("-nb_top", metavar="<int>",
                        help="The top hits to download",
                        default=1000)
    parser.add_argument('--no-proxy', action='store_true')
    parser.add_argument("-r", metavar="<str>",
                        help="The end point",
                        required=True)
    parser.add_argument("-c", metavar="<int>",
                        help="The length of the slice for the request",
    args = parser.parse_args()
    return args


def on_gene(sgrna, gene):
    """
    Check if a given sgrna is on a given gene
    """
    return gene.start <= sgrna[0] and gene.end >= sgrna[1]


def coord_on_gene(list_coord, list_gene, ref):
    """
    For a reference organism check if all sgrna are on genes
    Return a CoordSeq object if at least one sgrna is on a gene
    """
    coord_seq = CoordSeq(ref=ref)
    for coord in list_coord:
        nb_on_gene = list(filter(lambda x: on_gene(coord, x), list_gene))
        if nb_on_gene:
            coord_seq.update(coord)
    if coord_seq.list_coord:
        return coord_seq
    return None


def check_on_gene(blast_file, dic_index, nb_gi):
    """
    Return a dictionary of ResumeSeq object for each sequence which is at least
    on one gene in one organism.
    """
    blastoutput = pickle.load(open(blast_file, "rb"))
    resume_seq = {}
    # For each sgrna
    for seq in dic_index:
        # Create a ResumeSeq to save info
        resume_seq[seq] = ResumeSeq()
        # For each organism find where an homologous gene has been found
        for org in blastoutput.org_names():
            tmp = {}
            for ref in blastoutput.ref_names(org):
                # Find the list of coordinates on gene
                coord_seq = coord_on_gene(dic_index[seq].genomes_dict[org][ref],
                                          blastoutput.gene_coords(org, ref), ref)
                if coord_seq:
                    tmp[org].append(coord_seq)
            resume_seq[seq].update(tmp)
        # If no coordinate in the ResumeSeq object, then delete the sgrna from the dictionary
        if not resume_seq[seq].dic_org:
            del resume_seq[seq]
        else:
            resume_seq[seq].cal_proportion(nb_gi)
    # If resume_seq is not empty, then at least one sgrna is on one gene
    if resume_seq:
        return resume_seq
    # None sgrna on homologous gene
    print("Progam terminated&No hits on genes")
    sys.exit()


if __name__ == "__main__":
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    DIC_INDEX, NB_HITS = pp.parse_setcompare_out(PARAM.f, None)
    if DIC_INDEX:
        LEN_SEQ = int(PARAM.sl) + int(len(PARAM.pam))
        DIC_HITS = OrderedDict()
        for rank in DIC_INDEX:
            sequence = decoding.decode(rank, ["A", "T", "C", "G"], LEN_SEQ)
            DIC_HITS[sequence] = dspl.Hit(sequence, DIC_INDEX[rank])
        # Search coordinates for each sgrna in each organism
        DIC_HITS = pp.couchdb_search(DIC_HITS, GENOMES_IN, PARAM.r, int(PARAM.c), PARAM.no_proxy)
        # Keep sgrna which are on gene
        RESUME_SEQ = check_on_gene(PARAM.blast, DIC_HITS, len(GENOMES_IN))
        # Sort sgrna by the proportion of organism containing this sgrna on a gene
        RESUME_SEQ = sorted(RESUME_SEQ, key=lambda hit: hit.proportion, reverse=True)
        dspl.display_hits(RESUME_SEQ, GENOMES_IN, GENOMES_NOTIN,
                          PARAM.pam, int(PARAM.sl), ".", int(PARAM.nb_top))
        print(','.join(GENOMES_NOTIN))
        print("TASK_KEY")
        print(len(RESUME_SEQ))
    else:
        print("Progam terminated&No hits remain")

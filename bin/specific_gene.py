#!/usr/bin/env python3

import argparse
import pickle
import wordIntegerIndexing as decoding
import display_result as dspl
import post_processing as pp
from parse_blast import BlastReport, BlastHit


class ResumeSeq(object):
    """docstring for ResumeSeq."""
    def __init__(self, dic_org):
        self.proportion = 0
        self.dic_org = dic_org

    def __repr__(self):
        return "Proportion : {}\tTotal occurences : {}".format(self.proportion, self.total_occ())

    def cal_proportion(self, nb_total):
        self.proportion = len(self.dic_org) / int(nb_total)
        return self.proportion

    def total_occ(self, org_name):
        return self.dic_org[org].total_occ()

    def total_occ(self):
        occ = 0
        for org_name in self.dic_org:
            occ += self.total_occ(org_name)
        return occ

    def update(self, dict_org):
        if dict_org:
            self.dic_org.update(dict_org)

    def org_names(self):
        return list(self.dic_org.keys())

    def write(self, genomes_in):
        to_write = ""
        for gi in genomes_in:
            if gi in self.org_names():
                for ref in self.dic_org[gi]:
                    to_write += ref + ':' + ','.join(self.dic_org[gi][ref]) + ';'
            else:
                to_write += "None ;"
        to_write += self.total_occ()
        return to_write.strip(";")

    def list_ref(self, org_name):
        list_ref = [{"ref": ref, "coords": self.dic_org[org_name][ref].list_coord} for ref in self.dic_org[org_name]]
        return list_ref

    def list_occ(self):
        list_occurences = [{'org': genome, 'all_ref': self.list_ref(genome)} for genome in self.org_names()]
        return list_occurences



class CoordSeq(object):
    """docstring for CoordSeq."""
    def __init__(self):
        self.list_coord = []
        self.list_occ = []

    def __init__(self, coords, occs):
        self.list_coord = coords
        self.list_occ = occs

    def __repr__(self):
        return "Total occurences : {}\tFor {} coordinates".format(self.total_occ(), len(self))

    def update(self, coord, occ):
        if coord not in self.list_coord:
            self.list_coord.append(coord)
            self.list_occ.append(occ)
            return True
        return False

    def total_occ():
        return sum(self.list_occ)

    def __len__(self):
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
    parser.add_argument("-nb_top", metavar="<int>",
                        help="The top hits to download",
                        default=1000)
    args = parser.parse_args()

    return args


def on_gene(sgrna, gene):
    return gene.start <= sgrna[0] and gene.end >= sgrna[1]


def coord_on_gene(list_coord, list_gene):
    """
    Definition
    """
    coord_seq = CoordSeq()
    for coord in list_coord:
        nb_on_gene = len(list(filter(lambda x: on_gene(coord, x), list_gene)))
        if nb_on_gene != 0:
            coord_seq.update(coord, nb_on_gene)
    if len(coord_seq) != 0:
        return coord_seq
    return None


def check_on_gene(blast_file, dic_index, nb_gi):
    """
    La sortie de blast est plus courte : peut être qu'un seul organisme aura un
    gène homologue.
    1. Boucler sur les séquences
    2. Ne garder que les organismes présents dans blast
    3. Vérifier que sgrna est sur un gène tout en checkant la référence
    """
    blastoutput = pikle.load(open(blast_file, "rb"))
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
                coord_seq = coord_on_gene(dic_index[seq].genomes_dict[org][ref], blastoutput.gene_coords(org, ref))
                if coord_seq:
                    if not tmp:
                        tmp[org] = {}
                    tmp[org][ref] = coord_seq
            resume_seq[seq].update(tmp)
        if len(resume_seq[seq].dic_org) == 0:
            del resume_seq[seq]
        else:
            resume_seq[seq].cal_proportion(nb_gi)
    if resume_seq:
        return resume_seq
    else:
        print("Progam terminated&No hits on genes")
        sys.exit()


if __name__ == "__main__":
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    DIC_INDEX, NB_HITS = pp.parse_setcompare_out(PARAM.f, None)
    if DIC_INDEX:
        LEN_SEQ = int(PARAM.sl) + int(len(PARAM.pam))
        for rank in DIC_INDEX:
            sequence = decoding.decode(rank, ["A", "T", "C", "G"], len_seq)
            DIC_HITS[sequence] = dspl.Hit(sequence, DIC_INDEX[rank])
            # Search coordinates for each sgrna in each organism
        DIC_HITS = pp.couchdb_search(DIC_HITS, GENOMES_IN, PARAM.r, int(PARAM.c), PARAM.no_proxy)

        RESUME_SEQ = check_on_gene(PARAM.blast, DIC_HITS, len(GENOMES_IN))
        RESUME_SEQ = sorted(RESUME_SEQ, key=lambda hit:hit.proportion, reverse=True)
        dspl.display_hits(RESUME_SEQ, GENOMES_IN, GENOMES_NOTIN,
                          PARAM.pam, int(PARAM.sl), ".", int(PARAM.nb_top))
        print(','.join(GENOMES_NOTIN))
        print("TASK_KEY")
        print(len(RESUME_SEQ))
    else:
        print("Progam terminated&No hits remain")

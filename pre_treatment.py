"""
Pre-processing to add a new genome to the database
"""

import os
import argparse
import json
import sys
import tarfile
import shutil
import pickle
from Bio import SeqIO
from ete3 import NCBITaxa
import common_functions as cf


class Lineage:
    """
    Object which resume the full lineage
    """
    def __init__(self):
        self.species = "No specie"
        self.genus = "No genus"
        self.family = "No family"
        self.order = "No order"
        self.classe = "No class"
        self.phylum = "No phylum"


def valid_file(parser, filename):
    """
    Check if the file exists
    """
    if not os.path.isfile(filename):
        parser.error("The file {} does not exist !".format(filename))
    else:
        return filename


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Pre-treatment to import new genomes")
    parser.add_argument("-async", action='store_true')
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the reference genome folder")
    parser.add_argument("-file", metavar="FILE", type=lambda x: valid_file(parser, x),
                        help="The path to the fasta file", required=True)
    parser.add_argument("-gcf", metavar="<str>",
                        help="The GCF assembly ID ",
                        required=True)
    parser.add_argument("-asm", metavar="<str>",
                        help="The ASM assembly ID",
                        required=True)
    parser.add_argument("-taxid", metavar="<str>",
                        help="The taxon ID",
                        required=True)
    args = parser.parse_args()
    return args


def setup(ref):
    """
    Create directories for the fasta and the indexation of the new genome
    """
    os.mkdir("reference_genomes")
    os.mkdir("reference_genomes/fasta/")
    os.mkdir("reference_genomes/fasta/" + ref)
    os.mkdir("reference_genomes/index2/")
    os.mkdir("reference_genomes/index2/" + ref)


def set_dic_taxid(param):
    """
    Create dictionnary for json file and gzip the fasta file
    """
    # Retrieve the name of the genome
    seq_record = SeqIO.parse(param.file, "fasta")
    name = next(seq_record).description
    name = name.split(",")[0]
    ref = param.gcf + "_" + param.asm

    with open(param.rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)

    # Check if the reference is not in the dic_ref, so in the database
    for organism in dic_ref:
        if dic_ref[organism][0] == ref:
            sys.exit("ERROR : This genome is already in the database")
        elif dic_ref[organism][1] == param.taxid:
            sys.exit("ERROR : This taxon ID is already in the database")

    dic_taxid = {}
    dic_taxid[ref] = param.taxid
    # Retrieve all taxon id present in the database
    for name_gcf in dic_ref:
        ref_tmp = dic_ref[name_gcf][0]
        dic_taxid[ref_tmp] = dic_ref[name_gcf][1]

    dic_ref[name + ' ' + param.gcf] = [ref, param.taxid]

    setup(ref)
    shutil.copyfile(param.file, "reference_genomes/fasta/" + ref +
                    "/" + ref + "_genomic.fna")
    # Write the new json file with the new genome
    json.dump(dic_ref, open(param.rfg + "/genome_ref_taxid.json", 'w'), indent=4)
    return dic_taxid, ref, name


def index_bowtie_blast(ref_new):
    """
    Build bowtie indexation
    """
    print("INDEX")
    os.system("bowtie2-build reference_genomes/fasta/" + ref_new + "/" +
              ref_new + "_genomic.fna reference_genomes/index2/" + ref_new +
              "/" + ref_new)
    os.system("makeblastdb -in reference_genomes/fasta/" + ref_new + "/" +
              ref_new + "_genomic.fna -dbtype nucl")


def create_lineage_objects(dic_tax):
    """
    Retrieve species, genus, family, order, classe and phylum for each genome :
    new or old.
    Then, return a dictionnary of lineage object containing these informations
    """
    ncbi = NCBITaxa()
    dic_lineage = {}
    count = 0
    for ref in dic_tax:
        lineage_object = Lineage()
        tax_ref = dic_tax[ref]
        try:
            lineage = ncbi.get_lineage(tax_ref)
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            for i in ranks:
                if ranks[i] == 'species':
                    lineage_object.species = names[i]
                elif ranks[i] == 'genus':
                    lineage_object.genus = names[i]
                elif ranks[i] == 'family':
                    lineage_object.family = names[i]
                elif ranks[i] == 'order':
                    lineage_object.order = names[i]
                elif ranks[i] == 'class':
                    lineage_object.classe = names[i]
                elif ranks[i] == 'phylum':
                    lineage_object.phylum = names[i]
            dic_lineage[ref] = (lineage_object, count)
            count += 1
        except Exception as e:
            with open("problem_taxon.log", "a") as filout:
                filout.write(tax_ref + "\n")

    return dic_lineage


def distance_dic(dic_lineage, ref_new, param):
    """
    Calcul the distance between genomes according to their species, genus,
    family, order, classe or phylum.
    """
    with open(param.rfg + "/distance_dic.json", "r") as json_data:
        dic = json.load(json_data)
    dic[ref_new] = {}
    for ref1 in dic_lineage:
        if ref_new == ref1:
            dic[ref_new][ref1] = 0
        elif dic_lineage[ref_new][0].species == dic_lineage[ref1][0].species:
            dic[ref_new][ref1] = 1
        elif dic_lineage[ref_new][0].genus == dic_lineage[ref1][0].genus:
            dic[ref_new][ref1] = 2
        elif dic_lineage[ref_new][0].family == dic_lineage[ref1][0].family:
            dic[ref_new][ref1] = 3
        elif dic_lineage[ref_new][0].order == dic_lineage[ref1][0].order:
            dic[ref_new][ref1] = 4
        elif dic_lineage[ref_new][0].classe == dic_lineage[ref1][0].classe:
            dic[ref_new][ref1] = 5
        elif dic_lineage[ref_new][0].phylum == dic_lineage[ref1][0].phylum:
            dic[ref_new][ref1] = 6
        else:
            dic[ref_new][ref1] = 7
    return dic


def distance_matrix(dic_taxid, ref_new, param):
    """
    Create the distance matrix to know if genomes are close or not
    """
    print('DISTANCE MATRIX')
    dic_lineage = create_lineage_objects(dic_taxid)
    dist_dic = distance_dic(dic_lineage, ref_new, param)
    json.dump(dist_dic, open(param.rfg + "/distance_dic.json", "w"), indent=4)


def json_tree(bdd_path):
    """
    Definition
    """
    print('JSON TREE')
    os.system('python3 bin/tax2json.py ' + bdd_path)


def construct_in(fasta_path, organism, organism_code, param):
    """
    Construct the sequences for first organism,
    with python regular expression research
    """
    pam = "NGG"
    non_pam_motif_length = 20
    fasta_file = (fasta_path + '/' + organism_code + '/' +
                  organism_code + '_genomic.fna')
    sgrna = "N" * non_pam_motif_length + pam

    seq_dict = {}

    for genome_seqrecord in SeqIO.parse(fasta_file,'fasta'):
        genome_seq = genome_seqrecord.seq
        ref = genome_seqrecord.id
        seq_list_forward = cf.find_sgrna_seq(str(genome_seq),
                                             cf.reverse_complement(sgrna))
        seq_list_reverse = cf.find_sgrna_seq(str(genome_seq), sgrna)
        for indice in seq_list_forward:
            end = indice + len(pam) + non_pam_motif_length
            seq = genome_seq[indice:end].reverse_complement()
            seq = str(seq)
            if seq not in seq_dict:
                seq_dict[seq] = {organism: {}}
            if ref not in seq_dict[seq][organism]:
                seq_dict[seq][organism][ref] = []
            seq_dict[seq][organism][ref].append('+(' + str(indice+1) + ',' +
                                                str(end) + ')')

        for indice in seq_list_reverse:
            end = indice + len(pam) + non_pam_motif_length
            seq = genome_seq[indice:end]
            seq = str(seq)
            if seq not in seq_dict:
                seq_dict[seq] = {organism: {}}
            if ref not in seq_dict[seq][organism]:
                seq_dict[seq][organism][ref] = []
            seq_dict[seq][organism][ref].append('-(' + str(indice+1) + ',' +
                                                str(end) + ')')

    pickle_file = (param.rfg + "/pickle/" + organism.replace("/", "_") + "." + organism_code + ".p")
    print(pickle_file)
    pickle.dump(seq_dict, open(pickle_file, "wb"), protocol=3)



def compress(param, ref):
    # Compress the fasta file in the database
    with tarfile.open(param.rfg + "/fasta/" + ref + ".tar.gz", "w:gz") as tar:
        tar.add("reference_genomes/fasta/" + ref + "/" + ref + "_genomic.fna")
        tar.add("reference_genomes/fasta/" + ref + "/" + ref + "_genomic.fna.nhr")
        tar.add("reference_genomes/fasta/" + ref + "/" + ref + "_genomic.fna.nin")
        tar.add("reference_genomes/fasta/" + ref + "/" + ref + "_genomic.fna.nsq")

    # Compress the fasta file in the database
    with tarfile.open(param.rfg + "/index2/" + ref + ".tar.gz", "w:gz") as tar:
        tar.add("reference_genomes/index2/" + ref + "/" + ref + ".1.bt2")
        tar.add("reference_genomes/index2/" + ref + "/" + ref + ".2.bt2")
        tar.add("reference_genomes/index2/" + ref + "/" + ref + ".3.bt2")
        tar.add("reference_genomes/index2/" + ref + "/" + ref + ".4.bt2")
        tar.add("reference_genomes/index2/" + ref + "/" + ref + ".rev.1.bt2")
        tar.add("reference_genomes/index2/" + ref + "/" + ref + ".rev.2.bt2")


if __name__ == '__main__':
    PARAM = args_gestion()
    # Create dictionnary with all taxon ID
    DIC_TAXID, REF_NEW, NAME = set_dic_taxid(PARAM)
    # Index the new genome
    index_bowtie_blast(REF_NEW)
    # Compress the indexation and the fasta file
    compress(PARAM, REF_NEW)
    # Calcul distances between new and old genomes
    distance_matrix(DIC_TAXID, REF_NEW, PARAM)
    construct_in("reference_genomes/fasta", NAME + " " + PARAM.gcf, REF_NEW, PARAM)

    json_tree(PARAM.rfg)
    # Delete temporary directory
    shutil.rmtree("reference_genomes")

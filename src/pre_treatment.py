"""
Pre-processing to add a new genome to the database
"""

import os
import argparse
import json
import sys
import shutil
import pickle
import re
from Bio import SeqIO
from ete3 import NCBITaxa
import wordIntegerIndexing as decoding


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


# ARGUMENTS GESTION
def valid_fasta_file(parser, filename):
    """
    Check if the file exists and if it is a standard fasta file
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        parser.error("The file {} does not exist !".format(filename))
    # Try to open it and check if it is a fasta file
    try:
        fasta = next(SeqIO.parse(filename, "fasta"))
        # Check if the sequence only contains ACTG letters
        if any(letter not in "ACTG" for letter in fasta.seq):
            parser.error("Be careful, the sequence is not a nucleotide sequence !")
    except Exception as err:
        parser.error("The file {} is not a fasta file !".format(filename))
    # If all is good, return the filename
    return filename


def valid_taxid(parser, taxid):
    """
    Check if the taxon id given by the user is in the NCBI taxonomy
    database
    """
    ncbi = NCBITaxa()
    try:
        ncbi.get_lineage(taxid)
        return taxid
    except Exception as err:
        parser.error("The taxon id given ({}) is not in the NCBI taxonomy database !".format(taxid))


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    command = sys.argv[1] if len(sys.argv) > 1 else None
    # Argparsing
    parser = argparse.ArgumentParser(description="Pre-treatment to import new genomes")
    subparsers = parser.add_subparsers(help='commands')

    # Create metafile : pickle and index
    metafile_parser = subparsers.add_parser("metafile", help="Create the pickle file of the genome")
    metafile_parser.add_argument("-rfg", metavar="<str>",
                                 help="The path to the reference genome folder",
                                 required=True)
    metafile_parser.add_argument("-file", metavar="FILE",
                                 type=lambda x: valid_fasta_file(metafile_parser, x),
                                 help="The path to the fasta file",
                                 required=True)
    metafile_parser.add_argument("-gcf", metavar="<str>",
                                 help="The GCF assembly ID ", required=True)
    metafile_parser.add_argument("-asm", metavar="<str>",
                                 help="The ASM assembly ID", required=True)
    metafile_parser.add_argument("-taxid", type=lambda x: valid_taxid(metafile_parser, x),
                                 help="The taxon ID", required=True)

    # Add to the database
    add_bdd_parser = subparsers.add_parser("add", help="Add the pickle file to the database")
    add_bdd_parser.add_argument("-file", metavar="FILE",
                                type=lambda x: valid_taxid(add_bdd_parser, x),
                                help="The path to the pickle file to add to the database",
                                required=True)

    # Create metafile : pickle and index and add the pickle file to the database
    all_parser = subparsers.add_parser("all", help="Create the pickle file of the genome")
    all_parser.add_argument("-rfg", metavar="<str>",
                            help="The path to the reference genome folder",
                            required=True)
    all_parser.add_argument("-file", metavar="FILE", type=lambda x: valid_fasta_file(all_parser, x),
                            help="The path to the fasta file",
                            required=True)
    all_parser.add_argument("-gcf", metavar="<str>",
                            help="The GCF assembly ID ", required=True)
    all_parser.add_argument("-asm", metavar="<str>",
                            help="The ASM assembly ID", required=True)
    all_parser.add_argument("-taxid", type=lambda x: valid_taxid(all_parser, x),
                            help="The taxon ID", required=True)

    args = parser.parse_args()
    if hasattr(args, "taxid"):
        valid_taxid(parser, args.taxid)

    return args, command


# COMPARE NEW GENOME TO OLD TO FIND ITS BRANCH IN THE TOPOLOGY TREE
def set_dic_taxid(filename, gcf, asm, taxid, rfg):
    """
    Create dictionnary for json file and gzip the fasta file
    """
    # Retrieve the name of the genome
    seq_record = SeqIO.parse(filename, "fasta")
    name = next(seq_record).description
    name = name.split(",")[0]
    ref = gcf + "_" + asm

    with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)

    # Check if the reference is not in the dic_ref, so in the database
    references = [ref_ref[0] for ref_ref in dic_ref.values()]
    tax_ids = [id[1] for id in dic_ref.values()]
    if ref in references:
        sys.exit("ERROR : This genome is already in the database")
    if taxid in tax_ids:
        sys.exit("ERROR : This taxon ID is already in the database")

    dic_taxid = {}
    dic_ref[name + ' ' + gcf] = [ref, taxid]
    # Retrieve all taxon id present in the database
    dic_taxid = {dic_ref[name_gcf][0]: dic_ref[name_gcf][1] for name_gcf in dic_ref}

    shutil.copyfile(filename, rfg + "/fasta/" +
                    ref + "_genomic.fna")
    # Write the new json file with the new genome
    json.dump(dic_ref, open(rfg + "/genome_ref_taxid.json", 'w'), indent=4)
    return dic_taxid, ref, name


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
        except Exception as err:
            with open(".problem_taxon.log", "a") as filout:
                filout.write(tax_ref + "\n")

    return dic_lineage


def distance_dic(dic_lineage, ref_new, rfg):
    """
    Calcul the distance between genomes according to their species, genus,
    family, order, classe or phylum.
    """
    with open(rfg + "/distance_dic.json", "r") as json_data:
        dic = json.load(json_data)
    dic[ref_new] = {}
    for ref1 in dic_lineage:
        if ref_new == ref1:
            dic[ref_new][ref1] = 0
        elif dic_lineage[ref_new][0].species == dic_lineage[ref1][0].species:
            dic[ref_new][ref1] = 1
            dic[ref1][ref_new] = 1
        elif dic_lineage[ref_new][0].genus == dic_lineage[ref1][0].genus:
            dic[ref_new][ref1] = 2
            dic[ref1][ref_new] = 2
        elif dic_lineage[ref_new][0].family == dic_lineage[ref1][0].family:
            dic[ref_new][ref1] = 3
            dic[ref1][ref_new] = 3
        elif dic_lineage[ref_new][0].order == dic_lineage[ref1][0].order:
            dic[ref_new][ref1] = 4
            dic[ref1][ref_new] = 4
        elif dic_lineage[ref_new][0].classe == dic_lineage[ref1][0].classe:
            dic[ref_new][ref1] = 5
            dic[ref1][ref_new] = 5
        elif dic_lineage[ref_new][0].phylum == dic_lineage[ref1][0].phylum:
            dic[ref_new][ref1] = 6
            dic[ref1][ref_new] = 6
        else:
            dic[ref_new][ref1] = 7
            dic[ref1][ref_new] = 7
    return dic


def distance_matrix(dic_taxid, ref_new, rfg):
    """
    Create the distance matrix to know if genomes are close or not
    """
    print('DISTANCE MATRIX')
    dic_lineage = create_lineage_objects(dic_taxid)
    dist_dic = distance_dic(dic_lineage, ref_new, rfg)
    json.dump(dist_dic, open(rfg + "/distance_dic.json", "w"), indent=4)


# FIND ALL SGRNA
def build_expression(seq):
    """
    Build the regular expression by replacing the letter by iupac code or
    the letter itself
    """
    result = ''
    iupac_code = {'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
                  'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
                  'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
    for letter in seq:
        if letter in iupac_code:
            result = result + iupac_code[letter]
        else:
            result = result + letter
    return result


def find_sgrna_seq(seq, pam):
    """
    Uses Regular expression matching of the pam motif to the reference genome
    to get the start positions (0-based) of each match
    """
    reg_exp = build_expression(pam)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices


def reverse_complement2(sequence):
    """
    Function for turning a 5'-3' nucleotidic sequence into
    its 5'-3' reverse complement.
    """
    rev_comp = []
    for idx in range(len(sequence) - 1, -1, -1):
        if sequence[idx] == 'A':
            rev_comp = rev_comp + ['T']
        elif sequence[idx] == 'C':
            rev_comp = rev_comp + ['G']
        elif sequence[idx] == 'G':
            rev_comp = rev_comp + ['C']
        elif sequence[idx] == 'T':
            rev_comp = rev_comp + ['A']
        else:
            rev_comp = rev_comp + ['N']
    return "".join(rev_comp)


def construct_in(organism, organism_code, rfg, pam="NGG", non_pam_motif_length=20):
    """
    Construct the sequences for first organism,
    with python regular expression research
    """
    print("SEARCH SGRNA")
    fasta_file = (rfg + '/fasta/' + organism_code + '_genomic.fna')
    sgrna = "N" * non_pam_motif_length + pam

    seq_dict = {}

    for genome_seqrecord in SeqIO.parse(fasta_file, "fasta"):
        genome_seq = genome_seqrecord.seq
        ref = genome_seqrecord.id
        seq_list_forward = find_sgrna_seq(str(genome_seq),
                                          reverse_complement2(sgrna))
        seq_list_reverse = find_sgrna_seq(str(genome_seq), sgrna)
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

    pickle_file = (rfg + "/pickle/" + organism.replace("/", "_") + ".p")
    pickle.dump(seq_dict, open(pickle_file, "wb"), protocol=3)
    return pickle_file


def indexation(name_file, rfg, pickle_file):
    """
    Code each sgrna sequence to a rank. With this, it is easier to compare
    sequences of several genomes
    """
    name_file = name_file.replace("/", "_")
    target_file = rfg + "/index/" + name_file + '.index'

    total = decoding.indexPickle(pickle_file, target_file)
    print("Successfully indexed", total, "words\nfrom:", name_file, "\ninto:", target_file)


# CONSTRUCT THE NEW JSON TOPOLOGY TREE
def json_tree(bdd_path):
    """
    Take the entire list of genomes and create the topology tree in json format
    """
    print('JSON TREE')
    os.system('python3 src/tax2json.py ' + bdd_path)


if __name__ == '__main__':
    PARAM, COMMAND = args_gestion()
    if COMMAND == "metafile" or COMMAND == "all":
        # Create dictionnary with all taxon ID
        DIC_TAXID, REF_NEW, NAME = set_dic_taxid(PARAM.file, PARAM.gcf, PARAM.asm,
                                                 PARAM.taxid, PARAM.rfg)

        # Calcul distances between new and old genomes
        distance_matrix(DIC_TAXID, REF_NEW, PARAM.rfg)
        # the fasta file was copied in the tmp directory ./reference_genomes
        PICKLE_FILE = construct_in(NAME + " " + PARAM.gcf, REF_NEW, PARAM.rfg)
        indexation(NAME + " " + PARAM.gcf, PARAM.rfg, PICKLE_FILE)
        json_tree(PARAM.rfg)

    if COMMAND == "add" or COMMAND == "all":
        pass

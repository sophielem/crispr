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


# ARGUMENTS GESTION
def valid_fasta_file(parser, filename):
    """
    Check if the file exists and if it is a standard fasta file
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        parser.error("Program terminated&The file {} does not exist !".format(filename))
    # Try to open it and check if it is a fasta file
    try:
        fasta = next(SeqIO.parse(filename, "fasta"))
    except Exception as err:
        parser.error("Program terminated&The file {} is not a fasta file !".format(filename))
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
        parser.error("Program terminated&The taxon id given ({})\
                      is not in the NCBI taxonomy database !".format(taxid))


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


# CHECK IF THE TAXON ID IS ALREADY PRESENT AND COPY THE FASTA FILE
def set_dic_taxid(filename, gcf, asm, taxid, rfg):
    """
    Create dictionnary for json file and gzip the fasta file
    """
    # Retrieve the name of the genome
    seq_record = SeqIO.parse(filename, "fasta")
    name = next(seq_record).description
    name = name.split(",")[0]
    name = name.replace("/", "_")
    ref = gcf + "_" + asm

    with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)

    # Check if the reference is not in the dic_ref, so in the database
    references = [ref_ref[0] for ref_ref in dic_ref.values()]
    tax_ids = [id[1] for id in dic_ref.values()]
    if ref in references:
        sys.exit("Program terminated&ERROR : This genome is already in the database")
    if taxid in tax_ids:
        sys.exit("Program terminated&ERROR : This taxon ID is already in the database")

    dic_ref[name + ' ' + gcf] = [ref, taxid]
    # Write the new json file with the new genome
    json.dump(dic_ref, open(rfg + "/genome_ref_taxid.json", 'w'), indent=4)
    shutil.copyfile(filename, rfg + "/fasta/" + ref + "_genomic.fna")

    return ref, name


# FIND ALL SGRNA
def complement_seq(sequence):
    """
    Function for turning a 5'-3' nucleotidic sequence into
    its 5'-3' reverse complement.
    """
    rev_comp = []
    complement = {"A" : "T", "C" : "G", "T" : "A", "G" : "C"}
    rev_comp = [complement[sequence[idx]] if sequence[idx] in complement.keys() else "N" for idx in range(len(sequence) - 1, -1, -1)]
    return "".join(rev_comp)


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


def find_indices_sgrna(seq, pam):
    """
    Uses Regular expression matching of the pam motif to the reference genome
    to get the start positions (0-based) of each match
    """
    reg_exp = build_expression(pam)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices


def find_sgrna_seq(seq_list, len_seq, reverse, str_reverse, seq_dict, genome_seq, organism, ref):
    """
    Uses start index and the length of the sequence to get the sequence and
    fill the dictionary containing all sgRNA sequences with coordinates
    """
    for indice in seq_list:
        end = indice + len_seq
        seq = genome_seq[indice:end] if reverse else genome_seq[indice:end].reverse_complement()
        seq = str(seq)
        if seq not in seq_dict:
            seq_dict[seq] = {organism: {}}
        if ref not in seq_dict[seq][organism]:
            seq_dict[seq][organism][ref] = []
        seq_dict[seq][organism][ref].append(str_reverse + str(indice+1) + ',' +
                                            str(end) + ')')
    return seq_dict


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
        seq_list_forward = find_indices_sgrna(str(genome_seq),
                                              complement_seq(sgrna))
        seq_list_reverse = find_indices_sgrna(str(genome_seq), sgrna)

        seq_dict = find_sgrna_seq(seq_list_forward, len(pam) + non_pam_motif_length,
                                  False, "+(", seq_dict, genome_seq, organism, ref)
        seq_dict = find_sgrna_seq(seq_list_reverse, len(pam) + non_pam_motif_length,
                                  True, "-(", seq_dict, genome_seq, organism, ref)

    pickle_file = (rfg + "/pickle/" + organism + ".p")
    pickle.dump(seq_dict, open(pickle_file, "wb"), protocol=3)
    return pickle_file


def indexation(name_file, rfg, pickle_file):
    """
    Code each sgrna sequence to a rank. With this, it is easier to compare
    sequences of several genomes
    """
    target_file = rfg + "/index/" + name_file + '.index'

    total = decoding.indexPickle(pickle_file, target_file)
    print("Successfully indexed", total, "words\nfrom:", name_file, "\ninto:", target_file)


# CONSTRUCT THE NEW JSON TOPOLOGY TREE
def json_tree(bdd_path):
    """
    Take the entire list of genomes and create the topology tree in json format
    """
    print('JSON TREE')
    os.system('python3 lib/tax2json.py ' + bdd_path)


if __name__ == '__main__':
    PARAM, COMMAND = args_gestion()
    if COMMAND == "metafile" or COMMAND == "all":
        # Create dictionnary with all taxon ID
        REF_NEW, NAME = set_dic_taxid(PARAM.file, PARAM.gcf, PARAM.asm,
                                      PARAM.taxid, PARAM.rfg)
        PICKLE_FILE = construct_in(NAME + " " + PARAM.gcf, REF_NEW, PARAM.rfg)
        indexation(NAME + " " + PARAM.gcf, PARAM.rfg, PICKLE_FILE)
        json_tree(PARAM.rfg)
        print("Created metafile for {}".format(NAME))

    if COMMAND == "add" or COMMAND == "all":
        pass
        print("Add {} to the database. Ready for request".format(NAME))

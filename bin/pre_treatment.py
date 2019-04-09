"""
Pre-processing to add a new genome to the database
"""

import os
import argparse
import json
import sys
import shutil
import subprocess
from Bio import SeqIO
import word_detect
import wordIntegerIndexing as decoding
import display_result as dspl
import get_assembly_infos as gai

DEBUG = True


# ARGUMENTS GESTION
def exists_file(parser, filename):
    """
    Check if the file exists and if it is a standard fasta file
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        parser.error("Program terminated&The file {} does not exist !".format(filename))


def valid_fasta_file(parser, filename):
    """
    Check if the file exists and if it is a standard fasta file
    """
    # Check if the file exists
    exists_file(parser, filename)
    # Try to open it and check if it is a fasta file
    try:
        next(SeqIO.parse(filename, "fasta"))
    except Exception as err:
        parser.error("Program terminated&The file {} is not a fasta file !".format(filename))
    # If all is good, return the filename
    return filename


def check_metafile_exist(rfg, basename_file):
    """
    Check if the pickle and index file exist
    """
    return (os.path.exists(rfg + "/genome_index/" + basename_file + ".index") and
            os.path.exists(rfg + "/genome_pickle/" + basename_file + ".p"))


def parse_arg(subparser, command):
    """
    Treat arguments for subparsers
    """
    subparser.add_argument("-file", metavar="FILE",
                           type=lambda x: valid_fasta_file(subparser, x),
                           help="The path to the pickle file to add to the database")
    subparser.add_argument("-rfg", metavar="<str>",
                           help="The path to the reference genome folder",
                           required=True)
    subparser.add_argument("-url", metavar="<str>",
                           help="The end point", required=True)
    subparser.add_argument("-size", metavar="int",
                           help="Maximal number of entry to add at a time")
    subparser.add_argument("-min", metavar="int",
                           help="Index of the first file to add to database")
    subparser.add_argument("-max", metavar="int",
                           help="Index of the last file to add to database")
    subparser.add_argument("-tree", metavar="<str>",
                           help="Path to the json tree", required=True)
    subparser.add_argument("-m", metavar="<str>",
                           help="DB volumes end-point mapper", required=True)
    if command == "add":
        subparser.add_argument("-dir", metavar="<str>",
                               help="The path to the pickle file to add to the database")
    return subparser


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

    # Add to the database
    add_bdd_parser = subparsers.add_parser("add", help="Add the pickle file to the database")
    add_bdd_parser = parse_arg(add_bdd_parser, command)

    # Create metafile : pickle and index and add the pickle file to the database
    all_parser = subparsers.add_parser("all", help="Create the pickle file of the genome")
    all_parser = parse_arg(all_parser, command)

    args = parser.parse_args()
    if command == "add":
        if args.file and args.dir:
            parser.error("Program terminated&Missing file or directory")
        elif not args.file and not args.dir:
            parser.error("Program terminated&Choose file or directory, not both")
    return args, command


def check_genome_exists(filename, gcf, asm, taxid, rfg):
    """
    Create dictionnary for json file and gzip the fasta file
    """
    # Retrieve the name of the genome
    ref = gcf + "_" + asm

    with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)

    # Check if the reference is not in the dic_ref
    references = [ref_ref[0] for ref_ref in dic_ref.values()]
    tax_ids = [id[1] for id in dic_ref.values()]
    if ref in references:
        print("Program terminated&ERROR : This genome is already in the database")
        sys.exit()
    if taxid in tax_ids:
        print("Program terminated&ERROR : This taxon ID is already in the database")
        sys.exit()

    shutil.copyfile(filename, rfg + "/genome_fasta/" + ref + "_genomic.fna")

    return ref


def create_metafile(fasta_file, pickle_file, index_file, name):
    """
    Create pickle and index files
    """
    dspl.eprint("--- Search sgRNA ---")
    word_detect.construct_in(fasta_file, pickle_file)
    dspl.eprint("--- Indexation ---")
    total = decoding.indexAndOccurencePickle(pickle_file, index_file)
    dspl.eprint("Successfully indexed", total, "words\nfrom:", name, "\ninto:", index_file)
    dspl.eprint("SUCCESS&Created metafile for {}".format(name))


def set_dic_taxid(dic_index_files, error_list, rfg):
    """
    Add new genomes in the genome_reference with the taxon ID. Before, check
    if the genome is not in the error_list
    """
    with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)
    tmp_dic = {os.path.basename(filename).replace(".p", ""): dic_index_files[filename]
               for filename in dic_index_files
               if os.path.basename(filename).replace(".p", "") not in error_list}
    dic_ref.update(tmp_dic)
    # Write the new json file with the new genome
    json.dump(dic_ref, open(rfg + "/genome_ref_taxid.json", 'w'), indent=4)


def set_pickle_file_add(param):
    """
    Set the dictionnary of pickle files to add into the database
    """
    dic_pickle_files = {}
    # Create list with all fasta files
    list_files = os.listdir(param.dir) if param.dir else [param.file]
    # For each fasta file, retrieve the name of the pickle and index files
    # Then, check if they exist and create a list of path to pickle files
    for file_to_add in list_files:
        gcf, asm, taxid, name = gai.get_gcf_taxid(file_to_add)
        if not check_metafile_exist(param.rfg, name):
            print("Program terminated&Metafiles do not exist for {}".format(name))
            sys.exit()
        else:
            dic_pickle_files[param.rfg + "/genome_index/" + name + ".p"] = [gcf + "_" + asm, taxid]
    return dic_pickle_files


def set_pickle_file_all(rfg, name, ref, taxid):
    """
    Set the dictionnary of pickle files to add into the database
    """
    dic_pickle_files = {}
    dic_pickle_files[rfg + "/genome_index/" + name + ".p"] = [ref, taxid]
    return dic_pickle_files


if __name__ == '__main__':
    PARAM, COMMAND = args_gestion()

    if COMMAND == "metafile" or COMMAND == "all":
        GCF, ASM, TAXID, NAME = gai.get_gcf_taxid(PARAM.file)
        # Create dictionnary with all taxon ID
        REF_NEW = check_genome_exists(PARAM.file, GCF, ASM,
                                      TAXID, PARAM.rfg)
        create_metafile(PARAM.rfg + '/genome_fasta/' + REF_NEW + '_genomic.fna',
                        PARAM.rfg + "/genome_pickle/" + NAME + ".p",
                        PARAM.rfg + "/genome_index/" + NAME + '.index',
                        NAME)

    if COMMAND == "add" or COMMAND == "all":
        DIC_PICKLE_FILES = set_pickle_file_all(PARAM.rfg, NAME, REF_NEW, TAXID) if COMMAND == "all" else set_pickle_file_add(PARAM)
        if DEBUG: dspl.eprint(DIC_PICKLE_FILES)
        # Insertion into the database
        CMD_LINE = "python lib/couchBuild.py --url {} \
                    --map {} --data {}".format(PARAM.url, PARAM.m, list(DIC_PICKLE_FILES.keys()))
        PROCESS = subprocess.call(CMD_LINE.split())
        # If error for some organism
        if os.path.isfile("error_add_db.err"):
            with open("error_add_db.err", "r") as filin:
                ERROR_LIST = filin.read().splitlines()
        else:
            ERROR_LIST = []
        # If error to insert all genomes, stop the program
        if len(ERROR_LIST) == len(list(DIC_PICKLE_FILES.keys())):
            print("Program terminated&Problem to insert genome into database")
        else:
            set_dic_taxid(DIC_PICKLE_FILES, ERROR_LIST, PARAM.rfg)
            dspl.eprint("--- JSON TREE ---")
            CMD_LINE = 'python lib/tax2json.py ' + PARAM.rfg + " " + PARAM.tree
            PROCESS = subprocess.call(CMD_LINE.split())
            print("SUCCESS&Genomes are in the database, except : {}".format(ERROR_LIST))

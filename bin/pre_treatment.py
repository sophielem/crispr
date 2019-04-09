"""
Pre-processing to add a new genome to the database
"""

import os
import argparse
import json
import sys
import shutil
import re
from Bio import SeqIO
from tqdm import tqdm
import couchBuild
import word_detect
import wordIntegerIndexing as decoding
import pycouch.wrapper as couchDB
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
    for filename in dic_index_files:
        # The genome has been inserted into the database
        if filename not in error_list:
            attr = dic_index_files[filename]
            # Retrieve the original name of the organism
            filename = os.path.basename(filename)
            name = filename.replace(".index", "")
            # Add it to the reference dictionary
            dic_ref[name] = attr
    # Write the new json file with the new genome
    json.dump(dic_ref, open(rfg + "/genome_ref_taxid.json", 'w'), indent=4)


# def add_to_database(files_to_add, end_point, i_min, i_max, batch_size, rules_file):
#     """
#     Add genomes in the database
#     """
#     if i_max > len(files_to_add):
#         i_max = len(files_to_add)
#
#     couchDB.setServerUrl(end_point)
#     if not couchDB.couchPing():
#         print("Program terminated&Impossible to connect to the database")
#         sys.exit()
#
#     with open(rules_file, "r") as rules_filin:
#         couchDB.setKeyMappingRules(json.load(rules_filin))
#
#     error_files = []
#     for filename in files_to_add[i_min:i_max]:
#         gen = couchBuild.GenomeData(filename)
#         # print("globing", filename, "#items", len(c))
#
#         for i in tqdm(range(0, len(gen), batch_size)):
#             j = i + batch_size if i + batch_size < len(gen) else len(gen)
#             gen_splitted = gen[i:j]
#             res = couchDB.volDocAdd(gen_splitted)
#             for gen_splitted in res:
#                 if not 'ok' in gen_splitted:
#                     dspl.eprint("Error here ==>", str(gen_splitted))
#                     error_files.append(filename)
#     return error_files


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

    if COMMAND == "add":
        DIC_INDEX_FILES = {}
        # Create list with all fasta files
        LIST_FILES = os.listdir(PARAM.dir) if PARAM.dir else [PARAM.file]
        # For each fasta file, retrieve the name of the pickle and index files
        # Then, check if they exist and create a list of path to index files
        for file_to_add in LIST_FILES:
            GCF, ASM, TAXID, NAME = gai.get_gcf_taxid(file_to_add)
            if not check_metafile_exist(PARAM.rfg, NAME):
                print("Program terminated&Metafiles do not exist for {}".format(NAME))
                sys.exit()
            else:
                DIC_INDEX_FILES[PARAM.rfg + "/genome_index/" + NAME + ".index"] = [GCF + "_" + ASM, TAXID]

    elif COMMAND == "all":
        DIC_INDEX_FILES = {}
        DIC_INDEX_FILES[PARAM.rfg + "/genome_index/" + NAME + ".index"] = [REF_NEW, TAXID]

        if DEBUG: print("DIC_INDEX_FILES == > {}".format(DIC_INDEX_FILES))

    if COMMAND == "add" or COMMAND == "all":
        # ERROR_LIST = add_to_database(list(DIC_INDEX_FILES.keys()), PARAM.url,
        #                              PARAM.min, PARAM.max, PARAM.size, PARAM.m)
        os.system("python couchBuild.py --url {} --map {} --data {}".format(PARAM.url, PARAM.m, list(DIC_INDEX_FILES.keys())))
        if os.path.isfile("error_add_db.err"):
            with open("error_add_db.err", "r") as filin:
                ERROR_LIST = filin.read().splitlines()
        else:
            ERROR_LIST = []
        set_dic_taxid(DIC_INDEX_FILES, ERROR_LIST, PARAM.rfg)
        dspl.eprint('JSON TREE')
        os.system('python lib/tax2json.py ' + PARAM.rfg + " " + PARAM.tree)
        dspl.eprint("SUCCESS&Genomes are in the database, except : {}".format(ERROR_LIST))

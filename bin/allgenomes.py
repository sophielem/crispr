#!/usr/bin/env python3
"""
    Find sgRNA sequence for the CRISP/CAS9 method. Find sgRNA common to some
    bacteria and exclude some bacteria.
"""

from __future__ import print_function
import uuid
import time
import argparse
import sys
import os
import shutil
import pickle
from Bio import SeqIO
import common_functions as cf
import pre_treatment as pt

# Global variable
TASK_KEY = str(uuid.uuid1())
UNCOMPRESSED_GEN_DIR = None
ASYNC = False


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Allgenomes program.")
    parser.add_argument("-async", action='store_true')
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-gni", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    parser.add_argument("-sl", metavar="<int>",
                        help="The length of the sgrna, excluding pam")
    parser.add_argument("-pam", metavar="<str>",
                        help="The pam motif",
                        required=True)
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the reference genome folder")
    parser.add_argument("-cah", metavar="<str>",
                        help="The path to cache folder for the webservice")
    parser.add_argument("-add", action="store_true",
                        help="Add a temporary genome")
    parser.add_argument("-file", metavar="FILE", type=lambda x: pt.valid_file(parser, x),
                        help="The path to the fasta file")
    parser.add_argument("-gcf", metavar="<str>",
                        help="The GCF assembly ID ")
    parser.add_argument("-asm", metavar="<str>",
                        help="The ASM assembly ID")
    parser.add_argument("-taxid", metavar="<str>",
                        help="The taxon ID")
    args = parser.parse_args()
    if args.add and (not args.file or not args.gcf or not args.asm or not args.taxid):
        parser.error("-add requires -file, -gcf, -asm and -taxid arguments")
    return args


def set_global_env(param):
    """
    Set global variable
    """
    global TASK_KEY
    global UNCOMPRESSED_GEN_DIR
    global ASYNC  # A switch to properly configure sbatch delagation
    TASK_KEY = str(uuid.uuid1())
    cf.TASK_KEY = TASK_KEY
    UNCOMPRESSED_GEN_DIR = param.rfg
    ASYNC = param.async
    if ASYNC:
        workdir = os.getcwd()
    else:
        workdir = cf.setup_work_space(param)
    return workdir


def setup_application(parameters, dict_organism_code):
    """
    Processes the arguments to be usable for research
    """
    organisms_selected = parameters.gi.split('&')
    organisms_excluded = parameters.gni.split('&')
    organisms_selected = [i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded = [i for i in organisms_excluded if i in dict_organism_code]
    non_pam_motif_length = int(parameters.sl)
    return organisms_selected, organisms_excluded, non_pam_motif_length


def display_no_hits(genome, start_time, workdir, in_exc):
    """
    Print the name of the genome after which there are no more hits and
    exit the program
    """
    print("Program terminated&No hits remain after {} genome {}".format(in_exc, genome))
    end_time = time.time()
    total_time = end_time - start_time
    cf.eprint('TIME', total_time)
    cf.delete_used_files(workdir)
    sys.exit(1)


def sort_hits(hitlist):
    """
    Scoring of the hits found, where positive scores mean stronger
    sgrna constructs.
    Complete Hit objects with score and sort the hits with this scores.
    """
    cf.eprint('-- Sort hits --')
    for hit in hitlist:
        score = 0
        for genome in hit.genomes_dict:
            for ref in hit.genomes_dict[genome]:
                score += len(hit.genomes_dict[genome][ref])
        hit.score = score
    sorted_hitlist = sorted(hitlist, key=lambda hit: hit.score, reverse=True)
    return sorted_hitlist


def add_tmp_genome(param):
    """
    Enter a temporary genome in the database
    """
    shutil.copyfile(UNCOMPRESSED_GEN_DIR + "/genome_ref_taxid.json", UNCOMPRESSED_GEN_DIR + "/genome_ref_taxid_tmp.json")
    shutil.copyfile(UNCOMPRESSED_GEN_DIR + "/distance_dic.json", UNCOMPRESSED_GEN_DIR + "/distance_dic_tmp.json")
    dic_taxid, ref_new, name = pt.set_dic_taxid(param.file, param.gcf, param.asm, param.taxid, UNCOMPRESSED_GEN_DIR)
    pt.index_bowtie_blast(ref_new)
    # the fasta file was copied in the tmp directory ./reference_genomes
    pt.construct_in("reference_genomes/fasta", name + " " + param.gcf, ref_new, UNCOMPRESSED_GEN_DIR)
    pt.compress(UNCOMPRESSED_GEN_DIR, ref_new)
    # Delete temporary directory
    shutil.rmtree("reference_genomes")
    return name, ref_new


def remove_tmp_genome(param, name, ref_new):
    """
    Remove the temporary genome
    """
    # Enlever le dict taxid contenant le genome tmp et le remplacer par le fichier d'origine
    os.remove(UNCOMPRESSED_GEN_DIR + "/genome_ref_taxid.json")
    os.system("mv {}/genome_ref_taxid_tmp.json {}/genome_ref_taxid.json".format(UNCOMPRESSED_GEN_DIR, UNCOMPRESSED_GEN_DIR))
    os.remove(UNCOMPRESSED_GEN_DIR + "/distance_dic.json")
    os.system("mv {}/distance_dic_tmp.json {}/distance_dic.json".format(UNCOMPRESSED_GEN_DIR, UNCOMPRESSED_GEN_DIR))
    # Enlever les fichiers compressés de la BDD
    os.remove(UNCOMPRESSED_GEN_DIR + "/fasta/" + ref_new + ".tar.gz")
    os.remove(UNCOMPRESSED_GEN_DIR + "/index2/" + ref_new + ".tar.gz")
    organism = name + " " + param.gcf
    os.remove(UNCOMPRESSED_GEN_DIR + "/pickle/" + organism.replace("/", "_") + "." + ref_new + ".p")

    
def construction(fasta_path, pam, non_pam_motif_length, genomes_in, genomes_not_in, dict_org_code, workdir):
    """
    Principal algorithm to launch to find sgRNA common
    """
    start_time = time.time()
    num_thread = 4
    num_file = 4
    cf.eprint('## Search for', len(genomes_in), "included genomes and",
              len(genomes_not_in), 'excluded genomes with', num_thread,
              'thread(s) ##')

    cf.unzip_files(UNCOMPRESSED_GEN_DIR, genomes_in + genomes_not_in,
                   dict_org_code, workdir)

    if len(genomes_in) != 1:
        sorted_genomes = cf.sort_genomes(genomes_in, fasta_path, dict_org_code,
                                         False)
    else:
        sorted_genomes = genomes_in

    if len(genomes_not_in) >= 1:
        sorted_genomes_notin = cf.sort_genomes(genomes_not_in, fasta_path,
                                               dict_org_code, True)
    else:
        sorted_genomes_notin = []

    cf.eprint('-- RESEARCH --')
    # Research in first genome
    name_file = sorted_genomes[0].replace("/", "_")
    pickle_file = (UNCOMPRESSED_GEN_DIR + "/pickle/" + name_file +
                   "." + dict_org_code[sorted_genomes[0]][0] + ".p")
    if not os.path.isdir(UNCOMPRESSED_GEN_DIR + "/pickle"): os.mkdir(UNCOMPRESSED_GEN_DIR + "/pickle")
    if os.path.isfile(pickle_file):
        dic_seq = pickle.load(open(pickle_file, "rb"))
    else:
        dic_seq = pt.construct_in(fasta_path, sorted_genomes[0],
                                  dict_org_code[sorted_genomes[0]][0],
                                  UNCOMPRESSED_GEN_DIR)

    cf.eprint(str(len(dic_seq)) + ' hits in first included genome ' +
              sorted_genomes[0])
    list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file, workdir)

    # Order for research
    dist_file = UNCOMPRESSED_GEN_DIR + "/distance_dic.json"
    dist_dic = cf.read_json_dic(dist_file)
    cf.eprint('-- Determinate order for research -- ')
    list_order = cf.order_for_research(sorted_genomes[1:],
                                       sorted_genomes_notin, sorted_genomes[0],
                                       dict_org_code, dist_dic, [])
    # Execute the rest of the search
    for i in list_order:
        genome = i[0]
        if i[1] == 'notin':
            dic_seq = cf.bowtie_multithr(num_thread, list_fasta,
                                         dict_org_code[genome][0], dic_seq,
                                         genome, len(pam) + non_pam_motif_length, workdir, False)
        elif i[1] == 'in':
            dic_seq = cf.bowtie_multithr(num_thread, list_fasta,
                                         dict_org_code[genome][0], dic_seq,
                                         genome, len(pam) + non_pam_motif_length, workdir, True)

        in_exc = "include" if i[1] == "in" else "exclude"
        # No hits remain after exclude or include genome
        if not dic_seq:
            display_no_hits(genome, start_time, workdir, in_exc)

        cf.eprint(str(len(dic_seq)) + " hits remain after " + in_exc +
                  " genome " + genome)
        list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file, workdir)

    cf.delete_used_files(workdir)

    hit_list = cf.construct_hitlist(dic_seq)
    hit_list = sort_hits(hit_list)

    # Put results in local file for access via the interface.
    cf.write_to_file(genomes_in, genomes_not_in, hit_list[:10000], pam,
                     non_pam_motif_length, workdir)

    # Output formatting for printing to interface
    cf.output_interface(hit_list[:100], workdir)


def main():
    """
    Main function. Create the organism - code dictionnary and launch
    comparisons between selected genomes
    """
    start_time = time.time()

    parameters = args_gestion()

    WORKDIR = set_global_env(parameters)
    fasta_path = WORKDIR + '/reference_genomes/fasta'
    # Add a temporary genome in the database
    if parameters.add:
        name, ref_new = add_tmp_genome(parameters)

    # Keys: organism, values: genomic reference (ncbi)
    dict_organism_code = cf.read_json_dic(UNCOMPRESSED_GEN_DIR +
                                          '/genome_ref_taxid.json')
    organisms_selected, organisms_excluded, non_pam_motif_length = setup_application(parameters, dict_organism_code)
    print(','.join(organisms_excluded))
    cf.eprint('SELECTED', organisms_selected)
    cf.eprint('EXCLUDED', organisms_excluded)
    print(TASK_KEY)

    cf.eprint('---- CSTB complete genomes ----')
    cf.eprint('Parallelisation with distance matrix')
    construction(fasta_path, parameters.pam, non_pam_motif_length, organisms_selected,
                 organisms_excluded, dict_organism_code, WORKDIR)

    # Remove the temporary genome from the database
    if parameters.add:
        remove_tmp_genome(parameters, name, ref_new)

    end_time = time.time()
    total_time = end_time - start_time
    cf.eprint('TIME', total_time)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Definition
"""
# -*-coding:Utf-8 -*-

from __future__ import print_function
import uuid
import time
import argparse
import sys
import os
from queue import Queue
from threading import Thread
from Bio import SeqIO
import common_functions as cf

# Global variable
TASK_KEY = str(uuid.uuid1())
WORKDIR = None
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
    args = parser.parse_args()
    return args


def setup_application(parameters, dict_organism_code):
    """
    Processes the arguments to be usable for research
    """
    organisms_selected = parameters.gi.split('&')
    organisms_excluded = parameters.gni.split('&')
    organisms_selected = [i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded = [i for i in organisms_excluded if i in dict_organism_code]
    non_pam_motif_length = int(parameters.sl)

    return organisms_selected, organisms_excluded, parameters.pam, non_pam_motif_length


def construct_in(fasta_path, organism, organism_code, pam, non_pam_motif_length):
    """
    Construct the sequences for first organism,
    with python regular expression research
    """
    fasta_file = (fasta_path + '/' + organism_code + '/' +
                  organism_code + '_genomic.fna')
    sgrna = "N" * non_pam_motif_length + pam

    seq_dict = {}

    for genome_seqrecord in SeqIO.parse(fasta_file, 'fasta'):
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
    return seq_dict


def add_in_parallel(num_thread, list_fasta, organism_code, dic_seq, genome, len_sgrna):
    """
    Launch bowtie alignments for included genomes and treat the results,
    with parallelization (ie if 4 threads are selected, then 4 bowtie will be
    launch at the same time, with 4 subfiles of the beginning file.
    For included genomes, only the sequences matching exactly
    (no mismatch) with genome will be conserved.
    """
    def worker():
        """
        Definition
        """
        while True:
            e = q.get()
            fasta_file = e['input_fasta']
            num_str = str(e['num'])
            result_file = cf.run_bowtie(organism_code, fasta_file, num_str)
            res = open(result_file, 'r')
            dic_result = {}
            for line in res:
                # not a comment
                if line[0] != '@':
                    # can align these 2 reads
                    if line.split('\t')[2] != '*':
                        l_split = line.split('\t')
                        # CIGAR string representation of alignment
                        cigar = l_split[5]
                        if cigar == '23M':
                            # reads quality
                            read_quality = l_split[-2]
                            # Name of reference sequence where alignment occurs
                            ref = l_split[2]
                            if read_quality.split(':')[-1] == '23':
                                seq = l_split[0]
                                if seq not in dic_result:
                                    dic_result[seq] = dic_seq[seq]
                                if genome not in dic_result[seq]:
                                    dic_result[seq][genome] = {}
                                seq_align = l_split[9]
                                if seq != seq_align:
                                    strand = '+'
                                else:
                                    strand = '-'
                                start = l_split[3]
                                end = int(start) + len_sgrna - 1
                                if ref not in dic_result[seq][genome]:
                                    dic_result[seq][genome][ref] = []
                                coord = (strand + '(' + start + ',' +
                                         str(end) + ')')
                                dic_result[seq][genome][ref].append(coord)
            e['results'] = dic_result
            res.close()
            q.task_done()

    q = Queue()
    for i in range(num_thread):
        t = Thread(target=worker)
        t.daemon = True
        t.start()

    for e in list_fasta:
        q.put(e)

    q.join()
    total_results = {}
    for e in list_fasta:
        total_results.update(e['results'])

    return total_results


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


def construction(fasta_path, pam, non_pam_motif_length, genomes_in, genomes_not_in, dict_org_code):
    """
    Definition
    """
    start_time = time.time()
    num_thread = 4
    num_file = 4
    cf.eprint('## Search for', len(genomes_in), "included genomes and",
              len(genomes_not_in), 'excluded genomes with', num_thread,
              'thread(s) ##')

    cf.unzip_files(UNCOMPRESSED_GEN_DIR, genomes_in + genomes_not_in,
                   dict_org_code)

    if len(genomes_in) != 1:
        sorted_genomes = cf.sort_genomes(genomes_in, fasta_path, dict_org_code,
                                         False)
    else:
        sorted_genomes = genomes_in

    if len(genomes_not_in) >= 1:
        sorted_genomes_notin = cf.sort_genomes_desc(genomes_not_in, fasta_path,
                                                    dict_org_code, True)
    else:
        sorted_genomes_notin = []

    cf.eprint('-- RESEARCH --')

    # Research in first genome
    dic_seq = construct_in(fasta_path, sorted_genomes[0],
                           dict_org_code[sorted_genomes[0]][0], pam,
                           non_pam_motif_length)
    cf.eprint(str(len(dic_seq)) + ' hits in first included genome ' +
              sorted_genomes[0])
    list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file)

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
            dic_seq = cf.add_notin_parallel(num_thread, list_fasta,
                                            dict_org_code[genome][0], dic_seq)
            if not dic_seq:
                print("Program terminated&No hits remain after exclude genome "
                      + genome)
                end_time = time.time()
                total_time = end_time - start_time
                cf.eprint('TIME', total_time)
                cf.delete_used_files()
                sys.exit(1)
            cf.eprint(str(len(dic_seq)) + " hits remain after exclude genome "
                      + genome)
            list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file)

        elif i[1] == 'in':
            dic_seq = add_in_parallel(num_thread, list_fasta,
                                      dict_org_code[genome][0], dic_seq,
                                      genome, len(pam) + non_pam_motif_length)
            if not dic_seq:
                print("Program terminated&No hits remain after include genome "
                      + genome)
                end_time = time.time()
                total_time = end_time - start_time
                cf.eprint('TIME', total_time)
                cf.delete_used_files()
                sys.exit(1)
            cf.eprint(str(len(dic_seq)) + " hits remain after include genome "
                      + genome)
            list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file)

    cf.delete_used_files()

    hit_list = cf.construct_hitlist(dic_seq)
    hit_list = sort_hits(hit_list)

    # Put results in local file for access via the interface.
    cf.write_to_file(genomes_in, genomes_not_in, hit_list[:10000], pam,
                     non_pam_motif_length)

    # Output formatting for printing to interface
    cf.output_interface(hit_list[:100])


def set_globel_env(param):
    """
    Definition
    """
    global TASK_KEY
    global UNCOMPRESSED_GEN_DIR
    global WORKDIR
    global ASYNC  # A switch to properly configure sbatch delagation
    TASK_KEY = str(uuid.uuid1())
    cf.TASK_KEY = TASK_KEY
    UNCOMPRESSED_GEN_DIR = param.rfg
    ASYNC = param.async


def main():
    """
    Definition
    """
    global WORKDIR

    start_time = time.time()

    parameters = args_gestion()

    set_globel_env(parameters)

    if ASYNC:
        WORKDIR = os.getcwd()
        cf.WORKDIR = WORKDIR
    else:
        WORKDIR = cf.setup_work_space(parameters)
        cf.WORKDIR = WORKDIR

    fasta_path = WORKDIR + '/reference_genomes/fasta'
    # Keys: organism, values: genomic reference (ncbi)
    dict_organism_code = cf.read_json_dic(UNCOMPRESSED_GEN_DIR +
                                        '/genome_ref_taxid.json')
    # organisms_selected,organisms_excluded,pam,non_pam_motif_length=args_gestion(dict_organism_code)
    organisms_selected, organisms_excluded, pam, non_pam_motif_length = setup_application(parameters, dict_organism_code)
    print(','.join(organisms_excluded))
    cf.eprint('SELECTED', organisms_selected)
    cf.eprint('EXCLUDED', organisms_excluded)
    print(TASK_KEY)

    cf.eprint('---- CSTB complete genomes ----')
    cf.eprint('Parallelisation with distance matrix')
    construction(fasta_path, pam, non_pam_motif_length, organisms_selected,
                 organisms_excluded, dict_organism_code)
    end_time = time.time()
    total_time = end_time - start_time
    cf.eprint('TIME', total_time)


if __name__ == '__main__':
    main()

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
import common_functions as cf
import pre_treatment as pt
import scan_tree as st
sys.path.append("/Users/slematre/Documents/pyCouch/src")
import pycouch.wrapper as couchDB

# Global variable
TASK_KEY = str(uuid.uuid1())
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
    parser.add_argument("-taxid", type=lambda x: pt.valid_taxid(parser, x),
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
    global ASYNC  # A switch to properly configure sbatch delagation
    TASK_KEY = str(uuid.uuid1())
    cf.TASK_KEY = TASK_KEY
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


def display_no_hits(genome, start_time, workdir, in_exc):
    """
    Print the name of the genome after which there are no more hits and
    exit the program
    """
    print("Program terminated&No hits remain after {} genome {}"
          .format(in_exc, genome))
    end_time = time.time()
    total_time = end_time - start_time
    cf.eprint('TIME', total_time)


def display_hits(dic_seq, genomes_in, genomes_not_in, pam, non_pam_motif_length, workdir):
    """
    Sort hits and write output for interface
    """
    hit_list = cf.construct_hitlist(dic_seq)
    hit_list = sort_hits(hit_list)

    # Put results in local file for access via the interface.
    cf.write_to_file(genomes_in, genomes_not_in, hit_list[:10000], pam,
                     non_pam_motif_length, workdir)

    # Output formatting for printing to interface
    cf.output_interface(hit_list[:100], workdir)


def add_tmp_genome(param):
    """
    Enter a temporary genome in the database
    """
    shutil.copyfile(param.rfg + "/genome_ref_taxid.json",
                    param.rfg + "/genome_ref_taxid_tmp.json")
    shutil.copyfile(param.rfg + "/distance_dic.json",
                    param.rfg + "/distance_dic_tmp.json")
    dic_taxid, ref_new, name = pt.set_dic_taxid(param.file, param.gcf,
                                                param.asm, param.taxid,
                                                param.rfg)
    pt.index_bowtie_blast(ref_new)
    # the fasta file was copied in the tmp directory ./reference_genomes
    pt.construct_in("reference_genomes/fasta", name + " " + param.gcf,
                    ref_new, param.rfg)
    pt.compress(param.rfg, ref_new)
    pt.distance_matrix(dic_taxid, ref_new, param.rfg)
    # Delete temporary directory
    shutil.rmtree("reference_genomes")
    return name, ref_new


def remove_tmp_genome(param, name, ref_new):
    """
    Remove the temporary genome from json files and from the BDD
    """
    # Remvoe the dict taxid containing the tmp genome and replace it by
    # the original file. Same thing fot the distance file
    os.remove(param.rfg + "/genome_ref_taxid.json")
    os.system("mv {}/genome_ref_taxid_tmp.json {}/genome_ref_taxid.json"
              .format(param.rfg, param.rfg))
    os.remove(param.rfg + "/distance_dic.json")
    os.system("mv {}/distance_dic_tmp.json {}/distance_dic.json"
              .format(param.rfg, param.rfg))
    # Remove compressed files from the BDD
    os.remove(param.rfg + "/fasta/" + ref_new + ".tar.gz")
    os.remove(param.rfg + "/index2/" + ref_new + ".tar.gz")
    organism = name + " " + param.gcf
    os.remove(param.rfg + "/pickle/" + organism.replace("/", "_") +
              "." + ref_new + ".p")


def traverse_tree(genomes_in, dic_org_code, uncompressed_gen_dir):
    """
    Definition
    """
    # Dic to assocaite taxid with its organism name
    dic_taxid_org = {}
    list_taxid_in = []
    for sp in genomes_in:
        taxid = dic_org_code[sp][1]
        dic_taxid_org[taxid] = sp
        list_taxid_in.append(taxid)
    dic_node_in, list_taxid_in = st.find_node_complete(dic_org_code, list_taxid_in, uncompressed_gen_dir)
    # Keep organism name not under node
    genomes_in = [dic_taxid_org[taxid] for taxid in list_taxid_in]
    # cf.eprint("***** Number of hits : {}  *****".format(len(dic_seq)))
    return dic_node_in, genomes_in


def principal_search(list_order, list_fasta, dict_org_code, dic_seq, pam, non_pam_motif_length, workdir, start_time, num_thread=4, num_file=4):
    """
    The rest of the Search
    """
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
            break

        cf.eprint(str(len(dic_seq)) + " hits remain after " + in_exc +
                  " genome " + genome)
        list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file, workdir)
    return dic_seq


def genome_ordered_research(genomes_in, fasta_path, dict_org_code, uncompressed_gen_dir):
    """
    Sort genome for the research by distance between organism
    """
    if len(genomes_in) != 1:
        sorted_genomes = cf.sort_genomes(genomes_in, fasta_path, dict_org_code,
                                         False)
    else:
        sorted_genomes = genomes_in

    if sorted_genomes:
        # Order for research
        dist_file = uncompressed_gen_dir + "/distance_dic.json"
        dist_dic = cf.read_json_dic(dist_file)
        cf.eprint('\n\n-- Determinate order for research -- ')
        list_order = cf.order_for_research(sorted_genomes[1:],
                                           [], sorted_genomes[0],
                                           dict_org_code, dist_dic, [])
        return list_order, sorted_genomes[0]
    else:
        return [], None


def construction(fasta_path, pam, non_pam_motif_length, genomes_in,
                 genomes_not_in, dict_org_code, workdir, uncompressed_gen_dir):
    """
    Principal algorithm to launch to find sgRNA common
    """
    start_time = time.time()
    num_thread = 4
    num_file = 4
    cf.eprint('\n\n## Search for', len(genomes_in), "included genomes and",
              len(genomes_not_in), 'excluded genomes with', num_thread,
              'thread(s) ##')
    # Retrieve the dic from the pre-processing nodes and list of genomes not
    # under nodes pre-processed
    dic_seq, genomes_in = traverse_tree(genomes_in, dict_org_code, uncompressed_gen_dir)

    cf.unzip_files(uncompressed_gen_dir, genomes_in,
                   dict_org_code, workdir)
    list_order, genome_ref = genome_ordered_research(genomes_in,
                                                     fasta_path, dict_org_code,
                                                     uncompressed_gen_dir)
    cf.eprint('\n\n-- RESEARCH --')
    # dic_seq = {}
    if not dic_seq:
        # Research in first genome
        if not os.path.isdir(uncompressed_gen_dir + "/pickle"):
            os.mkdir(uncompressed_gen_dir + "/pickle")
        name_file = genome_ref.replace("/", "_")
        pickle_file = (uncompressed_gen_dir + "/pickle/" + name_file +
                       "." + dict_org_code[genome_ref][0] + ".p")
        # Load the pickle file
        if os.path.isfile(pickle_file):
            cf.eprint("Load pickle file")
            dic_seq = pickle.load(open(pickle_file, "rb"))
        # Search after sgRNA sequences in the genome by regex
        else:
            dic_seq = pt.construct_in(fasta_path, genome_ref,
                                      dict_org_code[genome_ref][0],
                                      uncompressed_gen_dir)
        cf.eprint(str(len(dic_seq)) + ' hits in first included genome ' +
                  genome_ref)
    else:
        list_order.insert(0, (genome_ref, "in"))

    # If they are genomes not under nodes, bowtie
    if genome_ref and list_order:
        cf.eprint("Research")
        list_fasta = cf.write_to_fasta_parallel(dic_seq, num_file, workdir)
        dic_seq = principal_search(list_order, list_fasta, dict_org_code, dic_seq,
                                   pam, non_pam_motif_length, workdir, start_time)
    dic_seq = couchdb_search(dic_seq, genomes_not_in)
    cf.delete_used_files(workdir)
    cf.eprint(str(len(dic_seq)) + ' hits remaining')
    return dic_seq


def couchdb_search(dic_seq, genomes_not_in):
    """
    Definition
    """
    couchDB.setServerUrl('http://localhost:1234')
    if couchDB.couchPing():
        sgrna = list(dic_seq.keys())
        start = time.time()
        results = couchDB.bulkRequestByKey(sgrna, "crispr_dvl")
        print("Time for the request : " + str(time.time() - start))
        # For each request
        for rslt in results["results"]:
            list_name = []
            request = rslt["id"]
            doc = rslt["docs"][0]
            # The request if ok
            if list(doc.keys())[0] == 'ok':
                # For each key in request result, so in json file
                for json_dic in doc["ok"]:
                    if json_dic not in ["_id", "_rev"]:
                        # name of organism
                        list_name.append(json_dic)
                if set(list_name).intersection(set(genomes_not_in)):
                    del dic_seq[request]
        return dic_seq
    else:
        sys.exit("Error : can't connect to the database")


def main():
    """
    Main function. Create the organism - code dictionnary and launch
    comparisons between selected genomes
    """
    start_time = time.time()

    parameters = args_gestion()
    UNCOMPRESSED_GEN_DIR = parameters.rfg
    workdir = set_global_env(parameters)
    fasta_path = workdir + '/reference_genomes/fasta'
    # Add a temporary genome in the database
    if parameters.add:
        cf.eprint('---- Add the temporary new genome ----')
        name, ref_new = add_tmp_genome(parameters)

    # Keys: organism, values: genomic reference (ncbi)
    dict_organism_code = cf.read_json_dic(UNCOMPRESSED_GEN_DIR +
                                          '/genome_ref_taxid.json')
    organisms_selected, organisms_excluded, non_pam_motif_length = setup_application(parameters, dict_organism_code)
    cf.eprint('SELECTED', organisms_selected)
    cf.eprint('EXCLUDED', organisms_excluded)
    print(TASK_KEY)

    cf.eprint('\n\n---- CSTB complete genomes ----')
    cf.eprint('Parallelisation with distance matrix')
    dic_seq = construction(fasta_path, parameters.pam, non_pam_motif_length, organisms_selected,
                           organisms_excluded, dict_organism_code, workdir, UNCOMPRESSED_GEN_DIR)
    if not dic_seq:
        sys.exit(1)
    else:
        display_hits(dic_seq, organisms_selected, organisms_excluded,
                     parameters.pam, non_pam_motif_length, workdir)

    # Remove the temporary genome from the database
    if parameters.add:
        cf.eprint('\n\n---- Remove the temporary new genome ----')
        remove_tmp_genome(parameters, name, ref_new)

    end_time = time.time()
    total_time = end_time - start_time
    cf.eprint('\n\n{}{}\n{}* TIME : {} *\n{}{}'.format(" "*10, "*"*29, " "*10,
                                                       total_time, " "*10, "*"*29))


if __name__ == '__main__':
    main()

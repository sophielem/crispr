#!/usr/bin/env python3
"""
Display results
"""

import sys
import json


class Hit():
    """
    Object Hit containing the CRISPR sequence, indexes where it is found
    and its associated score
    """
    def __init__(self, sequence, in_dic):
        self.sequence = sequence
        self.genomes_dict = in_dic
        self.score = 0


def eprint(*args, **kwargs):
    """
    For printing only on the terminal
    """
    print(*args, file=sys.stderr, **kwargs)


def construct_hitlist(dict_seq):
    '''
    Will construct an object Hit for each sequence in dictionnary,
    and store all Hits in a list. These function only fill the attributes
    sequence and genomes_dict of the object
    '''
    eprint('\n\n-- Construct final list --')
    hits_list = []
    count = 0
    for seq in dict_seq:
        count += 1
        new_hit = Hit(seq, dict_seq[seq])
        hits_list.append(new_hit)
    return hits_list


def sort_hits(hitlist):
    """
    Scoring of the hits found, where positive scores mean stronger
    sgrna constructs.
    Complete Hit objects with score and sort the hits with this scores.
    """
    eprint('-- Sort hits --')
    for hit in hitlist:
        score = 0
        for genome in hit.genomes_dict:
            for ref in hit.genomes_dict[genome]:
                score += len(hit.genomes_dict[genome][ref])
        hit.score = score
    sorted_hitlist = sorted(hitlist, key=lambda hit: hit.score, reverse=True)
    return sorted_hitlist


def write_to_file(genomes_in, genomes_not_in, hit_list, pam, non_pam_motif_length, workdir):
    """
    Write results in a file.
    The file is a tabulated file, with first column=sgrna sequence,
    then one column for each included organisms with list of positions
    of the sequence in each.
    """
    eprint('\n\n-- Write results to file --')
    rep_rslt_file = workdir + '/results_allgenome.txt'
    output = open(rep_rslt_file, 'w')
    gen_i = ','.join(genomes_in)
    gen_ni = ','.join(genomes_not_in)
    if gen_ni == '':
        gen_ni = 'None'
    output.write('#ALL GENOMES\n#Genomes included :' + gen_i +
                 ' ; Genomes excluded :' + gen_ni + '\n'+'#Parameters: pam:' +
                 pam + ' ; sgrna size:' + str(non_pam_motif_length) + '\n')
    output.write('sgrna sequence')
    for genome_i in genomes_in:
        output.write('\t' + genome_i)
    output.write('\n')
    for hit in hit_list:
        output.write(hit.sequence+'\t')
        for gi in genomes_in:
            for ref in hit.genomes_dict[gi]:
                to_write = ref + ':' + ','.join(hit.genomes_dict[gi][ref]) + ';'
        to_write = to_write.strip(';')
        output.write(to_write + '\n')
    output.close()


def create_list_ref(dic_ref):
    """
    Definition
    """
    list_ref = []
    for ref in dic_ref:
        dic = {"ref": ref, "coords": dic_ref[ref]}
        list_ref.append(dic)
    return list_ref


def create_list_occurences(dic_occurences):
    """
    Definition
    """
    list_occurences = []
    for genome in dic_occurences:
        dic_genome = {'org': genome, 'all_ref': create_list_ref(dic_occurences[genome])}
        list_occurences.append(dic_genome)
    return list_occurences


def output_interface(hit_list, workdir):
    """
    Reformat the results to create a json file.
    It will be parsed in javascript to display it in interface.
    """
    eprint('\n\n-- Construct results for graphical interface --')
    json_result_file = workdir + '/results.json'
    # print(json_result_file)
    list_dic = []

    for hit in hit_list:
        dic = {'sequence': hit.sequence, 'occurences': create_list_occurences(hit.genomes_dict)}
        list_dic.append(dic)

    list_dic.reverse()
    with open(json_result_file, 'w') as filout:
        json.dump(list_dic, filout, indent=4)


def display_hits(dic_seq, genomes_in, genomes_not_in, pam, non_pam_motif_length, workdir):
    """
    Sort hits and write output for interface
    """
    hit_list = construct_hitlist(dic_seq)
    hit_list = sort_hits(hit_list)
    
    # Put results in local file for access via the interface.
    write_to_file(genomes_in, genomes_not_in, hit_list[:10000], pam,
                  non_pam_motif_length, workdir)

    # Output formatting for printing to interface
    output_interface(hit_list[:100], workdir)

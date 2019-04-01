#!/usr/bin/env python3
"""
Display results
"""

import sys
import json
import itertools


class Hit():
    """
    Object Hit containing the CRISPR sequence, indexes where it is found
    and its associated score
    """
    def __init__(self, sequence, score):
        self.sequence = sequence
        self.genomes_dict = {}
        self.score = score

    def set_genomes_dict(self, dic_seq):
        """
        Set the dictionary containing coordinates
        """
        self.genomes_dict = dic_seq

    def write(self, genomes_in):
        to_write = ""
        for gi in genomes_in:
            for ref in self.genomes_dict[gi]:
                to_write += ref + ':' + ','.join(self.genomes_dict[gi][ref]) + ';'
        return to_write.strip(';')

    def list_ref(self, org_name):
        list_ref = [{"ref": ref, "coords": self.genomes_dict[org_name][ref]} for ref in self.genomes_dict[org_name]]
        return list_ref

    def list_occ(self):
        list_occurences = [{'org': genome, 'all_ref': self.list_ref(genome)} for genome in self.genomes_dict]
        return list_occurences


def eprint(*args, **kwargs):
    """
    For printing only on the terminal
    """
    print(*args, file=sys.stderr, **kwargs)


def write_to_file(genomes_in, genomes_not_in, dic_hits, pam, non_pam_motif_length, workdir, nb_top, hit_obj):
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
    i = 0
    for sgrna in dic_hits:
        if i == nb_top : break
        i += 1
        hit = dic_hits[sgrna]
        output.write(sgrna+'\t')
        to_write = hit.write(genomes_in)
        output.write(to_write + '\n')
    output.close()


def output_interface(dic_hits, workdir, nb_top):
    """
    Reformat the results to create a json file.
    It will be parsed in javascript to display it in interface.
    """
    # eprint('\n\n-- Construct results for graphical interface --')
    json_result_file = workdir + '/results.json'
    list_dic = [{'sequence': sgrna, 'occurences':dic_hits[sgrna].list_occ()} for i, sgrna in itertools.takewhile(lambda x: x[0] < nb_top, enumerate(dic_hits)) ]
    list_dic.reverse()
    with open(json_result_file, 'w') as filout:
        json.dump(list_dic, filout, indent=4)


def display_hits(dic_hits, genomes_in, genomes_not_in, pam, non_pam_motif_length, workdir, nb_top, hit_obj):
    """
    write output for interface
    """
    # Put results in local file for access via the interface.
    write_to_file(genomes_in, genomes_not_in, dic_hits, pam,
                  non_pam_motif_length, workdir, nb_top, hit_obj)

    # Output formatting for printing to interface
    output_interface(dic_hits, workdir, 100)

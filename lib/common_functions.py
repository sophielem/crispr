"""
    Definition
"""

from __future__ import print_function
import os
import sys
import json
import re
import subprocess
import tarfile
from queue import Queue
from threading import Thread
from glob import glob
from Bio import SeqIO


def humanize_bytes(bytes, precision=1):
    """Return a humanized string representation of a number of bytes.

    Assumes `from __future__ import division`.

    >>> humanize_bytes(1)
    '1 byte'
    >>> humanize_bytes(1024)
    '1.0 kB'
    >>> humanize_bytes(1024*123)
    '123.0 kB'
    >>> humanize_bytes(1024*12342)
    '12.1 MB'
    >>> humanize_bytes(1024*12342,2)
    '12.05 MB'
    >>> humanize_bytes(1024*1234,2)
    '1.21 MB'
    >>> humanize_bytes(1024*1234*1111,2)
    '1.31 GB'
    >>> humanize_bytes(1024*1234*1111,1)
    '1.3 GB'
    """
    abbrevs = (
        (1 << 50, 'PB'),
        (1 << 40, 'TB'),
        (1 << 30, 'GB'),
        (1 << 20, 'MB'),
        (1 << 10, 'kB'),
        (1, 'bytes')
    )
    if bytes == 1:
        return '1 byte'
    for factor, suffix in abbrevs:
        if bytes >= factor:
            break
    return '%.*f %s' % (precision, bytes / factor, suffix)


def _statdir(location):
    """
    Calcul the size of a given location
    """
    data = {}
    for sub_dir in glob(location + "/*/"):
        data[sub_dir] = 0
        for dirpath, dirnames, filenames in os.walk(sub_dir):
            for file in filenames:
                fpath = os.path.join(dirpath, file)
                data[sub_dir] += os.path.getsize(fpath)
        data[sub_dir] = humanize_bytes(data[sub_dir])
    return data


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


def reverse_complement(sequence):
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


def sort_genomes(list_genomes, fasta_path, dict_org_code, descending):
    """
    Sort genomes by ascending size
    """
    tmp_list = []
    for genome in list_genomes:
        len_genome = 0
        for seq_record in SeqIO.parse(fasta_path + '/' +
                                      dict_org_code[genome][0] + '/' +
                                      dict_org_code[genome][0] +
                                      '_genomic.fna', 'fasta'):
            len_genome += len(seq_record)
        tmp_list.append((len_genome, genome))
    # Sort by ascending or descending size
    genomes_sorted = [i[1] for i in sorted(tmp_list, key=lambda genome: genome[0], reverse=descending)]
    return genomes_sorted
    

def unzip_files(ref_gen_dir, list_genomes, dict_org_code):
    """
    Unzip required files for reserach
    """
    eprint('-- Unzip selected genomes --')
    eprint('location ' + WORKDIR)
    for genome in list_genomes:
        ref = dict_org_code[genome][0]
        extract_unzip_file(ref_gen_dir, ref, "/fasta/")
        extract_unzip_file(ref_gen_dir, ref, "/index2/")

    stats = _statdir(WORKDIR)
    eprint("Temporary file size\n" + str(stats))


def extract_unzip_file(ref_gen_dir, ref, folder):
    """
    Extract files from a tar.gz compression
    """
    tar = tarfile.open(ref_gen_dir + folder + ref + '.tar.gz')
    tar.extractall(WORKDIR + "/")
    tar.close()


def write_to_fasta_parallel(dic_seq, num_file):
    """
    Write sequences in fasta file, and separate its
    in several files if specified.
    """
    list_seq = list(dic_seq.keys())
    sep = len(dic_seq) // num_file
    list_dic_fasta = []
    for num in range(num_file):
        out = open(WORKDIR + '/sgrna' + str(num) + '.fa', 'w')
        list_dic_fasta.append({'num': num,
                               'input_fasta': WORKDIR + '/sgrna' + str(num) + '.fa',
                               'results': None
                              })
        i = 0
        for seq in list_seq:
            while i < sep:
                remove_seq = list_seq.pop()
                out.write('>' + remove_seq + '\n' + remove_seq + '\n')
                i += 1
        if num == num_file-1:
            if list_seq:
                for seq in list_seq:
                    out.write('>' + seq + '\n' + seq + '\n')
        out.close()
    return list_dic_fasta


def run_bowtie(organism_code, fasta_file, num):
    """
    Definition
    """
    result_file = WORKDIR + '/results_bowtie' + num + '.sam'
    bowtie_tab = ['bowtie2', '-x ' + WORKDIR + '/reference_genomes/index2/' +
                  organism_code + '/' + organism_code + ' -f ' + fasta_file +
                  ' -S ' + result_file + ' -L 13 -a --quiet ']
    subprocess.call(bowtie_tab)
    return result_file


def treat_bowtie_notin(result_file, e, dic_seq):
    """
    Definition
    """
    res = open(result_file, 'r')
    dic_result = {}
    for l in res:
        if l[0] != '@':
            if l.split('\t')[2] == '*':
                seq = l.split('\t')[0]
                dic_result[seq] = dic_seq[seq]
    e['results'] = dic_result
    res.close()


def add_notin_parallel(num_thread, list_fasta, organism_code, dic_seq):
    """
    Launch bowtie alignments for excluded genomes and treat the results,
    with parallelization (ie if 4 threads are selected,
    then 4 bowtie will be launch at the same time, with 4 subfiles
    of the beginning file.
    For excluded genomes, only the sequence NOT matching with genome
    will be conserved.
    """
    def worker():
        """
        Definition
        """
        while True:
            e = q.get()
            fasta_file = e['input_fasta']
            num_str = str(e['num'])
            result_file = run_bowtie(organism_code, fasta_file, num_str)
            treat_bowtie_notin(result_file, e, dic_seq)
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


def delete_used_files():
    """
    Delete unzipped files that have been used for research
    """
    eprint("Deleting stuff")
    eprint('-- Delete selected genomes --')
    out = open(WORKDIR + '/delete.sh', 'w')
    out.write('rm -r ' + WORKDIR + '/*.sam\nrm -r ' + WORKDIR +
              '/reference_genomes\nrm ' + WORKDIR + '/delete.sh\n')
    out.close()
    os.system('bash ' + WORKDIR + '/delete.sh')


def construct_hitlist(dict_seq):
    '''
    Will construct an object Hit for each sequence in dictionnary,
    and store all Hits in a list. These function only fill the attributes
    sequence and genomes_dict of the object
    '''
    eprint('-- Construct final list --')
    hits_list = []
    count = 0
    for seq in dict_seq:
        count += 1
        new_hit = Hit(seq, dict_seq[seq])
        hits_list.append(new_hit)
    return hits_list


def write_to_file(genomes_in, genomes_not_in, hit_list, pam, non_pam_motif_length):
    """
    Write results in a file.
    The file is a tabulated file, with first column=sgrna sequence,
    then one column for each included organisms with list of positions
    of the sequence in each.
    """
    eprint('-- Write results to file --')
    rep_rslt_file = WORKDIR + '/results_allgenome.txt'
    output = open(rep_rslt_file, 'w')
    gi = ','.join(genomes_in)
    gni = ','.join(genomes_not_in)
    if gni == '':
        gni = 'None'
    output.write('#ALL GENOMES\n#Genomes included :' + gi +
                 ' ; Genomes excluded :' + gni + '\n'+'#Parameters: pam:' +
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


def setup_work_space(parameters):
    """
    Create the work folder where all temporary or
    results files will be stored
    """
    work_folder = parameters.cah + '/' + TASK_KEY
    os.mkdir(work_folder)
    return work_folder


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


def output_interface(hit_list):
    """
    Reformat the results to create a json file.
    It will be parsed in javascript to display it in interface.
    """
    eprint('-- Construct results for graphical interface --')
    json_result_file = WORKDIR + '/results.json'
    # print(json_result_file)
    list_dic = []

    for hit in hit_list:
        dic = {'sequence': hit.sequence, 'occurences': create_list_occurences(hit.genomes_dict)}
        list_dic.append(dic)

    list_dic.reverse()
    with open(json_result_file, 'w') as filout:
        json.dump(list_dic, filout, indent=4)


def read_json_dic(path):
    """
    Load the dictionnary that contains genomes name associated with
    their references
    """
    with open(path) as json_data:
        dict_json = json.load(json_data)
    return dict_json


def order_for_research(list_in, list_notin, genome, dict_org_code, dist_dic, list_order):
    """
    Determinate the order of research, based on the distance matrix
    created for all genomes.
    The best comparisons to do are comparisons between similar genomes
    for the exclusion and between distant genomes for inclusion
    """
    ref1 = dict_org_code[genome][0]
    if list_in and list_notin:
        in_compare = -1
        for gi in list_in:
            ref2 = dict_org_code[gi][0]
            dist = dist_dic[ref1][ref2]
            if dist > in_compare:
                in_compare = dist
                in_compare_genome = gi
        notin_compare = 11
        for gni in list_notin:
            ref2 = dict_org_code[gni][0]
            dist = dist_dic[ref1][ref2]
            if dist < notin_compare:
                notin_compare = dist
                notin_compare_genome = gni

        if in_compare > notin_compare:
            new_genome = in_compare_genome
            list_order.append((new_genome, 'in'))
            list_in.remove(new_genome)
        else:
            new_genome = notin_compare_genome
            list_order.append((new_genome, 'notin'))
            list_notin.remove(new_genome)

    elif list_in:
        in_compare = -1
        for gi in list_in:
            ref2 = dict_org_code[gi][0]
            dist = dist_dic[ref1][ref2]
            if dist > in_compare:
                in_compare = dist
                in_compare_genome = gi
        new_genome = in_compare_genome
        list_order.append((new_genome, 'in'))
        list_in.remove(new_genome)

    elif list_notin:
        notin_compare = 11
        for gni in list_notin:
            ref2 = dict_org_code[gni][0]
            dist = dist_dic[ref1][ref2]
            if dist < notin_compare:
                notin_compare = dist
                notin_compare_genome = gni
        new_genome = notin_compare_genome
        list_order.append((new_genome, 'notin'))
        list_notin.remove(new_genome)
    else:
        return list_order

    return order_for_research(list_in, list_notin, new_genome, dict_org_code, dist_dic, list_order)

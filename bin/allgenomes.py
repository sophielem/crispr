#!/usr/bin/env python3

# -*-coding:Utf-8 -*-

from __future__ import print_function
import uuid
import time
import argparse
import sys
import os
from Bio import SeqIO
from queue import Queue
from threading import Thread
import common_functions as cf

#Global variable
TASK_KEY = str(uuid.uuid1())
WORKDIR = None
UNCOMPRESSED_GEN_DIR = None
ASYNC = False

def args_gestion():
    '''Take and treat arguments that user gives in command line'''
    ##Argparsing
    parser=argparse.ArgumentParser(description="Allgenomes program.")
    parser.add_argument("-async", action='store_true')
    parser.add_argument("-gi",metavar="<str>",help="The organisms to search inclusion in.",required=True)
    parser.add_argument("-gni",metavar="<str>",help="The organisms to search exclusion from",required=True)
    parser.add_argument("-sl",metavar="<int>",help="The length of the sgRNA, excluding PAM")
    parser.add_argument("-pam",metavar="<str>",help="The PAM motif",required=True)
    parser.add_argument("-rfg",metavar="<str>",help="The path to the reference genome folder")
    parser.add_argument("-cah",metavar="<str>",help="The path to cache folder for the webservice")
    args=parser.parse_args()

    return args


def setupApplication(parameters, dict_organism_code):
    '''Processes the arguments to be usable for research'''
    organisms_selected=parameters.gi.split('&')
    organisms_excluded=parameters.gni.split('&')
    organisms_selected=[i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded=[i for i in organisms_excluded if i in dict_organism_code]
    non_PAM_motif_length=int(parameters.sl)

    return organisms_selected, organisms_excluded, parameters.pam, non_PAM_motif_length

def construct_in(fasta_path, organism, organism_code, PAM, non_PAM_motif_length):
    '''Construct the sequences for first organism, with python regular expression research'''
    fasta_file=fasta_path + '/' + organism_code + '/' + organism_code +'_genomic.fna'

    sgRNA=''
    for i in range(non_PAM_motif_length):
        sgRNA+='N'
    sgRNA+=PAM

    seq_dict={}

    for genome_seqrecord in SeqIO.parse(fasta_file,'fasta'):
        genome_seq=genome_seqrecord.seq
        ref=genome_seqrecord.id
        seq_list_forward=cf.find_sgRNA_seq(str(genome_seq),cf.reverse_complement(sgRNA))
        seq_list_reverse=cf.find_sgRNA_seq(str(genome_seq),sgRNA)

        for indice in seq_list_forward:
            end=indice+len(PAM)+non_PAM_motif_length
            seq=genome_seq[indice:end].reverse_complement()
            seq=str(seq)
            if seq not in seq_dict:
                seq_dict[seq]={organism:{}}
            if ref not in seq_dict[seq][organism]:
                seq_dict[seq][organism][ref]=[]
            seq_dict[seq][organism][ref].append('+('+str(indice+1)+','+str(end)+')')

        for indice in seq_list_reverse:
            end=indice+len(PAM)+non_PAM_motif_length
            seq=genome_seq[indice:end]
            seq=str(seq)
            if seq not in seq_dict:
                seq_dict[seq]={organism:{}}
            if ref not in seq_dict[seq][organism]:
                seq_dict[seq][organism][ref]=[]
            seq_dict[seq][organism][ref].append('-('+str(indice+1)+','+str(end)+')')
    return seq_dict


def add_in_parallel(num_thread,list_fasta,organism_code,dic_seq,genome,len_sgrna):
    '''Launch bowtie alignments for included genomes and treat the results, with parallelization (ie if 4 threads are selected, then 4 bowtie will be launch at the same time, with 4 subfiles of the beginning file.
    For included genomes, only the sequences matching exactly (no mismatch0 with genome will be conserved.
    '''

    def worker():
        while True:
            e = q.get()
            fasta_file=e['input_fasta']
            num_str=str(e['num'])
            resultFile = cf.run_bowtie(organism_code, fasta_file, num_str)
            res=open(resultFile, 'r')
            dic_result={}
            for l in res:
                # not a comment
                if l[0]!='@':
                    # can align these 2 reads
                    if l.split('\t')[2]!='*':
                        l_split=l.split('\t')
                        #CIGAR string representation of alignment
                        cigar=l_split[5]
                        if cigar=='23M':
                            # reads quality
                            mm=l_split[-2]
                            # Name of reference sequence where alignment occurs
                            ref=l_split[2]
                            if mm.split(':')[-1]=='23':
                                seq=l_split[0]
                                if seq not in dic_result:
                                    dic_result[seq]=dic_seq[seq]
                                if genome not in dic_result[seq]:
                                    dic_result[seq][genome]={}
                                seq_align=l_split[9]
                                if seq != seq_align:
                                    strand='+'
                                else:
                                    strand='-'
                                start=l_split[3]
                                end=int(start)+len_sgrna-1
                                if ref not in dic_result[seq][genome]:
                                    dic_result[seq][genome][ref]=[]
                                coord=strand+'('+start+','+str(end)+')'
                                dic_result[seq][genome][ref].append(coord)
            e['results']=dic_result
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
    total_results={}
    for e in list_fasta:
        total_results.update(e['results'])

    return total_results

def sort_hits(hitlist):
    '''
    Scoring of the hits found, where positive scores mean stronger sgRNA constructs.
    Complete Hit objects with score and sort the hits with this scores.
    '''
    cf.eprint('-- Sort hits --')
    for hit in hitlist:
        score=0
        for genome in hit.genomes_Dict:
            for ref in hit.genomes_Dict[genome]:
                score+=len(hit.genomes_Dict[genome][ref])
        hit.score=score
    sorted_hitlist=sorted(hitlist,key=lambda hit:hit.score,reverse=True)
    return(sorted_hitlist)

def construction(fasta_path,PAM,non_PAM_motif_length,genomes_IN,genomes_NOT_IN,dict_org_code):
    start_time=time.time()
    start = time.time()
    num_thread=4
    num_file=4
    cf.eprint('## Search for',len(genomes_IN),"included genomes and",len(genomes_NOT_IN),'excluded genomes with',num_thread,'thread(s) ##')

    cf.unzip_files(UNCOMPRESSED_GEN_DIR, genomes_IN+genomes_NOT_IN, dict_org_code)

    if len(genomes_IN)!=1:
        sorted_genomes=cf.sort_genomes(genomes_IN,fasta_path,dict_org_code)
    else:
        sorted_genomes=genomes_IN

    if len(genomes_NOT_IN)>=1:
        sorted_genomes_notin=cf.sort_genomes_desc(genomes_NOT_IN,fasta_path,dict_org_code)
    else:
        sorted_genomes_notin=[]

    cf.eprint('-- RESEARCH --')

    #Research in first genome
    dic_seq=construct_in(fasta_path,sorted_genomes[0],dict_org_code[sorted_genomes[0]][0],PAM,non_PAM_motif_length)
    cf.eprint(str(len(dic_seq))+' hits in first included genome '+sorted_genomes[0])
    list_fasta=cf.write_to_fasta_parallel(dic_seq,num_file)

    #Order for research
    distFile = UNCOMPRESSED_GEN_DIR + "/distance_dic.json"
    dist_dic=cf.readJsonDic(distFile)
    cf.eprint('-- Determinate order for research -- ')
    list_order=cf.order_for_research(sorted_genomes[1:],sorted_genomes_notin,sorted_genomes[0],dict_org_code,dist_dic,[])

    #Execute the rest of the search
    for i in list_order:
        genome=i[0]
        if i[1]=='notin':
            dic_seq=cf.add_notin_parallel(num_thread,list_fasta,dict_org_code[genome][0],dic_seq)
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after exclude genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                cf.eprint('TIME',total_time)
                cf.delete_used_files()
                sys.exit(1)
            cf.eprint(str(len(dic_seq))+' hits remain after exclude genome '+genome)
            list_fasta=cf.write_to_fasta_parallel(dic_seq,num_file)

        elif i[1]=='in':
            dic_seq=add_in_parallel(num_thread,list_fasta,dict_org_code[genome][0],dic_seq,genome,len(PAM)+non_PAM_motif_length)
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after include genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                cf.eprint('TIME',total_time)
                cf.delete_used_files()
                sys.exit(1)
            cf.eprint(str(len(dic_seq))+' hits remain after include genome '+genome)
            list_fasta=cf.write_to_fasta_parallel(dic_seq,num_file)

    cf.delete_used_files()
    print(len(dic_seq))

    hit_list=cf.construct_hitlist(dic_seq)

    hit_list=sort_hits(hit_list)

    ##Put results in local file for access via the interface.
    cf.write_to_file(genomes_IN,genomes_NOT_IN,hit_list[:10000],PAM,non_PAM_motif_length)

    ##Output formatting for printing to interface
    cf.output_interface(hit_list[:100])


def setGlobalEnvironment(param):
    global TASK_KEY
    global UNCOMPRESSED_GEN_DIR
    global WORKDIR
    global ASYNC # A switch to properly configure sbatch delagation
    TASK_KEY=str(uuid.uuid1())
    cf.TASK_KEY=TASK_KEY
    UNCOMPRESSED_GEN_DIR=param.rfg
    ASYNC = param.async

def main():
    global WORKDIR

    start_time=time.time()

    parameters = args_gestion()

    setGlobalEnvironment(parameters)

    if ASYNC:
        WORKDIR = os.getcwd()
        cf.WORKDIR=WORKDIR
    else:
        WORKDIR=cf.setupWorkSpace(parameters)
        cf.WORKDIR=WORKDIR


    fasta_path  = WORKDIR + '/reference_genomes/fasta'

    dict_organism_code = cf.readJsonDic(UNCOMPRESSED_GEN_DIR + '/genome_ref_taxid.json')  ##Keys: organism, values: genomic reference (ncbi)
   # organisms_selected,organisms_excluded,PAM,non_PAM_motif_length=args_gestion(dict_organism_code)
    organisms_selected, organisms_excluded, PAM, non_PAM_motif_length = setupApplication(parameters, dict_organism_code)
    print(','.join(organisms_excluded))
    cf.eprint('SELECTED',organisms_selected)
    cf.eprint('EXCLUDED',organisms_excluded)
    print(TASK_KEY)

    cf.eprint('---- CSTB complete genomes ----')
    cf.eprint('Parallelisation with distance matrix')
    construction(fasta_path,PAM,non_PAM_motif_length,organisms_selected,organisms_excluded,dict_organism_code)
    end_time=time.time()
    total_time=end_time-start_time
    cf.eprint('TIME',total_time)

if __name__=='__main__':
    main()

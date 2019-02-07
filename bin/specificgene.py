#!/usr/bin/python

# -*-coding:Utf-8 -*-



from __future__ import print_function
import common_functions as cf
import time,argparse,uuid,os,subprocess,sys
from Bio.Blast import NCBIXML
from Queue import Queue
from threading import Thread
ASYNC = False

### METHODS FOR THE ARGUMENTS ###

def args_gestion():
    '''Take and treat arguments that user gives in command line'''
    parser = argparse.ArgumentParser(__file__, description="Specific gene program.")
    parser.add_argument("-async", action='store_true')
    parser.add_argument("-seq", help="String query sequence", required=True)
    parser.add_argument("-gi", help="List with the name(s) of genome(s)", required=True)
    parser.add_argument("-n", help="The research will be done on the n first bases of the gene", required=True,type=int)
    parser.add_argument("-gni", help='List with the name(s) of not in genome(s)', required=True)
    parser.add_argument('-ip',
                        help='identity percentage min for the research of homologous genes using blastn (default:70)', required=True,type=int)
    parser.add_argument('-pam',help='pam motif for sgrna', required=True)
    parser.add_argument('-sl',help='sgrna length (without pam motif)', required=True, type=int)
    parser.add_argument("-rfg",metavar="<str>",help="The path to the reference genome folder")
    parser.add_argument("-cah",metavar="<str>",help="The path to cache folder for the webservice")
    args = parser.parse_args()

    return args

def setup_application(parameters,dict_organism_code):
    organisms_selected=parameters.gi.split('&')
    organisms_excluded=parameters.gni.split('&')
    organisms_selected=[i for i in organisms_selected if i in dict_organism_code]
    organisms_excluded=[i for i in organisms_excluded if i in dict_organism_code]
    non_pam_motif_length=int(parameters.sl)

    return parameters.seq,organisms_selected,organisms_excluded,non_pam_motif_length,parameters.n,parameters.ip,parameters.pam


### METHODS FOR THE RESEARCH ###

def blast_to_find_all_genes(query_seq, genomes, identity_percentage_min,dic_genome):
    '''Launch blastn between query sequence and reference genome and all in genomes (separately).
    For origin genome, will only save the genes position, not the subject sequence (the query sequence is conserved for origin genome)
    For other genomes IN, will take the subject sequence found by blast as sequence for the organism and the sgrnas will be searched on it. It will also return gene's genomic position
    The hits must have a %id greater than the %id given by user
    If there is several blast hits, we only take the first (as defined by Blast)
    '''

    dic_genes={}
    size_gene=len(query_seq)
    file_query=WORKDIR+'/blast_request.fna'
    with open(file_query,'w') as f:
        f.write(query_seq)
    f.close()
    for i in range(len(genomes)):
        cf.eprint("Searching for sequence in",genomes[i])
        ref=dic_genome[genomes[i]][0]
        db_path = WORKDIR + "/reference_genomes/fasta/" + ref + '/' + ref + "_genomic.fna"
        blast_command = "blastn -db " + db_path + " -query " + file_query + " -outfmt 5"
        blast_output = os.popen(blast_command, 'r')
        blast_records = NCBIXML.parse(blast_output)
        for blast_record in blast_records:
            if blast_record.alignments == []:
                cf.eprint("No blast hits for",genomes[i])
                print("Program terminated&No blast hits for "+genomes[i]+". Your sequence probably does not have homologs in this genome.")
                exit(1)
            else:
                first_alignment=blast_record.alignments[0]
                first_hsp=first_alignment.hsps[0]
                identity_percentage=first_hsp.identities/len(first_hsp.query)
                if identity_percentage < (identity_percentage_min/100):
                    cf.eprint("Blast hit(s) do not meet stringency criterium: ID percentage. Program termination.")
                    print("Program terminated& Blast hit in "+genomes[i]+" do not meet stringency criterium ID percentage. Maybe try with a lower ID percentage.")
                    exit(1)
                else:
                    if i==0: #first genome, we want the query seq to search sgrna
                        if first_hsp.sbjct_start > first_hsp.sbjct_end: #gene is on reverse strand, + treatment to have position of the query seq and not the position of the sbjct hit (both can be not exactly the same, and we want to search sgrna in our query seq)
                            gene_start = first_hsp.sbjct_end - (size_gene - first_hsp.query_end) - 1
                            gene_end = first_hsp.sbjct_start + first_hsp.query_start - 2

                        else: #gene on forward strand, treatment to have position of the query seq
                            gene_start = first_hsp.sbjct_start - first_hsp.query_start
                            gene_end = first_hsp.sbjct_end + (size_gene - first_hsp.query_end) - 1

                        dic_genes[genomes[i]]=(gene_start,gene_end)

                    else: #other in genomes, now we take the hit subject as seq for sgrna searching
                        if first_hsp.sbjct_start > first_hsp.sbjct_end: #reverse strand
                            gene_start=first_hsp.sbjct_end-1 #blast is on 1-based and our program need 0-based
                            gene_end=first_hsp.sbjct_start-1
                            strand='-'
                        else:
                            gene_start=first_hsp.sbjct_start-1
                            gene_end=first_hsp.sbjct_end-1
                            strand='+'
                        dic_genes[genomes[i]]=(gene_start,gene_end)
    return dic_genes


def sgrna_on_gene(start_gene, end_gene, start_sgrna, end_sgrna):
    '''Check if a sgrna is in a gene (with both genomics positions)'''
    on_gene = False
    if start_gene <= start_sgrna and end_gene >= end_sgrna:
        on_gene = True
    return on_gene


def score_and_sort(hitlist):
    '''Score and sort the results of an analysis.
    For now, the score for 1 sgrna is the sum of all mismatchs + 1 if it's present in an excluded organism (no matter if it's exact match or not)
    Fill the attribute score of Hit objects.
    So best sgrna is the one that have lowest score.
    '''
    number_on_gene=0
    for hit in hitlist:
        dic_on_gene={}
        count_on_gene=0
        for genome in hit.genomes_dict:
            dic_on_gene[genome]=0
            count_genomes_on_gene=0
            for ref in hit.genomes_dict[genome]:
                list_coords=hit.genomes_dict[genome][ref]
                for coord in list_coords:
                    if coord.split(':')[-1]=='OnGene':
                        count_on_gene+=1
                        dic_on_gene[genome]+=1
        hit.score=count_on_gene
        count=0
        for genome in dic_on_gene:
            if dic_on_gene[genome]>=1:
                count+=1
        if count==len(hit.genomes_dict):
            number_on_gene+=1

    print(number_on_gene)
    sorted_hitlist=sorted(hitlist,key=lambda hit:hit.score,reverse=True)
    return sorted_hitlist

def search_sgrna_in_query_seq(query_seq,non_pam_motif_length, pam):

    sgrna=''
    for i in range(non_pam_motif_length):
        sgrna+='N'
    sgrna+=pam

    dic_seq={}

    seq_list_forward=cf.find_sgrna_seq(query_seq,cf.reverse_complement(sgrna))

    for indice in seq_list_forward:

        end=indice+len(pam)+non_pam_motif_length
        seq=cf.reverse_complement(query_seq[indice:end])
        dic_seq[seq]={}

    return dic_seq

def treat_bowtie_in(result_file,e,dic_seq,genome,len_sgrna,gene_coords):
    res=open(result_file, 'r')
    dic_result={}
    for l in res:
        if l[0]!='@':
            if l.split('\t')[2]!='*':
                l_split=l.split('\t')
                cigar=l_split[5]
                if cigar=='23M':
                    mm=l_split[-2]
                    if mm.split(':')[-1]=='23':
                        ref=l_split[2]
                        seq=l_split[0]
                        if seq not in dic_result:
                            dic_result[seq]=dic_seq[seq]
                        dic_result[seq][genome]={}
                        if ref not in dic_result[seq][genome]:
                            dic_result[seq][genome][ref]=[]
                        seq_align=l_split[9]
                        if seq != seq_align:
                            strand='+'
                        else:
                            strand='-'
                        start=l_split[3]
                        end=int(start)+len_sgrna-1
                        if sgrna_on_gene(gene_coords[0],gene_coords[1],int(start),int(end)):
                            info_gene='OnGene'
                        else:
                            info_gene='OffGene'
                        coord=strand+'('+start+','+str(end)+'):'+info_gene
                        dic_result[seq][genome][ref].append(coord)
    e['results']=dic_result
    res.close()

def add_in_parallel(num_thread,list_fasta,organism_code,dic_seq,genome,len_sgrna,gene_coords):
    '''Launch bowtie alignments for included genomes and treat the results, with parallelization (ie if 4 threads are selected, then 4 bowtie will be launch at the same time, with 4 subfiles of the beginning file.
    For included genomes, only the sequences matching exactly (no mismatch0 with genome will be conserved.
    '''
    def worker():
        while True:
            e = q.get()
            fasta_file=e['input_fasta']
            num_str=str(e['num'])
            result_file =cf.run_bowtie(organism_code, fasta_file, num_str)
            treat_bowtie_in(result_file,e,dic_seq,genome,len_sgrna,gene_coords)
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

def set_globel_env(param):
    global TASK_KEY
    global REF_GEN_DIR
    global WORKDIR
    global ASYNC

    TASK_KEY=str(uuid.uuid1())
    cf.TASK_KEY=TASK_KEY
    REF_GEN_DIR=param.rfg
    cf.REF_GEN_DIR=REF_GEN_DIR
    ASYNC = param.async


def do_search(query_seq, n, genome_list, dict_org_code, not_in_genome_list,
              identity_percentage_min,pam,sgrna_length):
    '''Launch the research with all parameters given by user. Principal function, it will call other functions to do the research'''

    start_time=time.time()
    num_threads=2
    fasta_path=WORKDIR+'/reference_genomes/fasta'
    len_all_sgrna=sgrna_length+len(pam)


    cf.unzip_files(genome_list+not_in_genome_list,dict_org_code) #unzip folder required for search

    #find homologous genes coordinates in each selected genome
    dic_genes=blast_to_find_all_genes(query_seq, genome_list, identity_percentage_min,dict_org_code)

    cf.eprint(dic_genes)

    #search sgrna sequences in given gene
    dic_seq=search_sgrna_in_query_seq(query_seq,sgrna_length,pam)

    #sort included and excluded genomes to optimize performances
    sorted_genomes_in=cf.sort_genomes(genome_list,fasta_path,dict_org_code)
    if not_in_genome_list:
        sorted_genomes_notin=cf.sort_genomes_desc(not_in_genome_list,fasta_path,dict_org_code)
    else:
        sorted_genomes_notin=[]

    #look for sequence previously found in first genome
    first_genome=sorted_genomes_in[0]
    list_fasta=cf.write_to_fasta_parallel(dic_seq,num_threads)
    dic_seq=add_in_parallel(num_threads,list_fasta,dict_org_code[first_genome][0],dic_seq,first_genome,len_all_sgrna,dic_genes[first_genome])
    list_fasta=cf.write_to_fasta_parallel(dic_seq,num_threads)

    #Determinate the order for further research, to optimize performance
    dist_file = REF_GEN_DIR + "/distance_dic.json"
    distDic=cf.read_json_dic(dist_file)
    list_order=cf.order_for_research(sorted_genomes_in[1:],sorted_genomes_notin,sorted_genomes_in[0],dict_org_code,distDic,[])

    #Execute the rest of the search, call different functions if genome is included or excluded.
    for i in list_order:
        genome=i[0]
        if i[1]=='notin':
            dic_seq = cf.bowtie_multithr(num_thread, list_fasta,
                                       dict_org_code[genome][0], dic_seq,
                                       genome, len(pam) + non_pam_motif_length, False)

            if len(dic_seq)==0:
                print("Program terminated&No hits remain after exclude genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                cf.eprint('TIME',total_time)
                cf.delete_used_files()
                sys.exit(1)
            cf.eprint(str(len(dic_seq))+' hits remain after exclude genome '+genome)
            list_fasta=cf.write_to_fasta_parallel(dic_seq,num_threads)

        elif i[1]=='in':
            dic_seq=add_in_parallel(num_threads,list_fasta,dict_org_code[genome][0],dic_seq,genome,len_all_sgrna,dic_genes[genome])
            if len(dic_seq)==0:
                print("Program terminated&No hits remain after include genome "+genome)
                end_time=time.time()
                total_time=end_time-start_time
                cf.eprint('TIME',total_time)
                cf.delete_used_files()
                sys.exit(1)
            cf.eprint(str(len(dic_seq))+' hits remain after include genome '+genome)
            list_fasta=cf.write_to_fasta_parallel(dic_seq,num_threads)


    cf.delete_used_files() # delete unzip folders previously unzipped

    print(len(dic_seq)) #give the number of hits to inferface

    #Sort and format results
    hit_list=cf.construct_hitlist(dic_seq)
    hit_list=score_and_sort(hit_list)
    cf.write_to_file(genome_list,not_in_genome_list,hit_list[:10000],pam,sgrna_length)
    cf.output_interface(hit_list[:100])


def main():
    global WORKDIR
    global ASYNC

    start_time=time.time()
    parameters = args_gestion()


   # cf.eprint('--->' + str(ASYNC))
    set_globel_env(parameters)

    if ASYNC:
        WORKDIR = os.getcwd()
        cf.WORKDIR=WORKDIR
    else:
        WORKDIR=cf.setup_work_space(parameters)
        cf.WORKDIR=WORKDIR

    cf.eprint('TASK_KEY',TASK_KEY)
    cf.eprint('REF_GEN_DIR',REF_GEN_DIR)
    cf.eprint('WORKDIR',WORKDIR)

    dict_organism_code = cf.read_json_dic(REF_GEN_DIR + '/genome_ref_taxid.json')
    query_seq,organisms_selected,organisms_excluded,non_pam_motif_length,n,identity_percentage_min,pam=setup_application(parameters,dict_organism_code)



    print(','.join(organisms_excluded)) # to give the not in genomes to interface
    print(TASK_KEY) #to give the workdir key to interface
    do_search(query_seq, n, organisms_selected, dict_organism_code, organisms_excluded, identity_percentage_min, pam, non_pam_motif_length)
    end_time=time.time()
    total_time=end_time-start_time
    cf.eprint('TIME',total_time)
    #eprint('CUMULATIVE_LENGTH',cumulative_length(genome_list,dic_genome))

if __name__ == "__main__":
    main()

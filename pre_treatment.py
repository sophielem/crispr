import os
from Bio import SeqIO
# from ete3 import NCBITaxa
import pickle
import argparse
from common_functions import reverse_complement
import re, json, sys
import tarfile
import shutil

class Lineage:
    def __init__(self):
        self.species="No specie"
        self.genus="No genus"
        self.family="No family"
        self.order="No order"
        self.classe='No class'
        self.phylum="No phylum"


def valid_file(parser, filename):
    """
    Check if the file exists
    """
    if not os.path.isfile(filename):
        parser.error("The file {} does not exist !".format(filename))
    else:
        return filename

def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Pre-treatment to import new genomes")
    parser.add_argument("-async", action='store_true')
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the reference genome folder")
    parser.add_argument("-file", metavar="FILE", type=lambda x: valid_file(parser, x),
                        help="The path to the fasta file", required=True)
    parser.add_argument("-gcf", metavar="<str>",
                        help="The GCF assembly ID ",
                        required=True)
    parser.add_argument("-asm", metavar="<str>",
                        help="The ASM assembly ID",
                        required=True)
    parser.add_argument("-taxid", metavar="<str>",
                        help="The taxon ID",
                        required=True)
    args = parser.parse_args()
    return args


def setup():
    os.system('mkdir reference_genomes')
    os.system('mkdir reference_genomes/fasta')
    os.system('mkdir reference_genomes/index2')


def fname(param):
    """
    Create dictionnary for json file and gzip the fasta file
    """
    # MOCHE
    seq_record = SeqIO.parse(param.file, "fasta")
    name = next(seq_record).description
    name = name.split(",")[0]
    ref = param.gcf + "_" + param.asm

    with open(param.rfg + "/genome_ref_taxid.json", "r") as json_data:
        dic_ref = json.load(json_data)

    # Check if the reference is not in the dic_ref, so in the database
    for organism in dic_ref:
        if dic_ref[organism][0] == ref:
            sys.exit("ERROR : This genome is already in the database")
        elif dic_ref[organism][1] == param.taxid:
            sys.exit("ERROR : This taxon ID is already in the database")

    dic_taxid={}
    dic_taxid[ref] = param.taxid
    dic_ref[name + ' ' + param.gcf] = [ref, param.taxid]

    os.mkdir("reference_genomes")
    os.mkdir("reference_genomes/fasta/")
    os.mkdir("reference_genomes/fasta/" + ref)
    shutil.copyfile(param.file, "reference_genomes/fasta/" + ref + "/" + ref + "_genomic.fna")
    with tarfile.open(param.rfg + "/fasta/" + ref + ".tar.gz", "w:gz") as tar:
        tar.add("reference_genomes/fasta/" + ref + "/" + ref + "_genomic.fna")
    shutil.rmtree("reference_genomes")

    json.dump(dic_ref,open(param.rfg + "/genome_ref_taxid.json",'w'),indent=4)
    return dic_taxid


def dic_download(ref_bacteria):
    print('DOWNLOAD')
    out=open('to_download.sh','w')
    taxfile=open('scripts/taxfile.txt','w')
    f=open(ref_bacteria,'r')
    cmd=''
    dic_taxid={}
    dic_ref={}
    for l in f:
        if l[0]!='#':
            l_split=l.rstrip().split('\t')
            name=l_split[7]
            name=name.replace("'",'')
            short_ref=l_split[0]
            taxfile.write(name+'\n')
            ftp_link=l_split[19]
            ref=ftp_link.split('/')[-1]
            taxid=l_split[5]
            dic_taxid[ref]=taxid
            dic_ref[name+' '+short_ref]=[ref,taxid]
            ftp_suffix=ftp_link.split('/')[-1]+'_genomic.fna'
            genome_link='https://'+ftp_link.split('//')[1]+'/'+ftp_suffix
            if ref+'.tar.gz' not in os.listdir('reference_genomes/fasta'):
                cmd+='mkdir reference_genomes/fasta/'+ref+'\n'
                cmd+='curl --remote-name --remote-time '+genome_link+'.gz\n'
                cmd+='mv '+ftp_suffix+'.gz reference_genomes/fasta/'+ref+'/\n'
                cmd+='gunzip reference_genomes/fasta/'+ref+'/'+ftp_suffix+'.gz\n'

    out.write(cmd)
    out.close()
    taxfile.close()
    f.close()
    json.dump(dic_ref,open('reference_genomes/genome_ref_taxid.json','w'),indent=4)
    #os.system('bash to_download.sh')
    os.system('rm to_download.sh')
    return dic_taxid


def index_bowtie_blast(list_ref):
    print('INDEX')
    out=open('index.sh','w')
    for ref in list_ref:
        cmd=''
        if ref not in os.listdir('reference_genomes/index2'):
            #cmd='gunzip reference_genomes/fasta/'+ftp_suffix+'.gz\n'
            cmd+='mkdir reference_genomes/index2/'+ref+'\n'
            cmd+='bowtie2-build reference_genomes/fasta/'+ref+'/'+ref+'_genomic.fna reference_genomes/index2/'+ref+'/'+ref+'\n'
            cmd+='makeblastdb -in reference_genomes/fasta/'+ref+'/'+ref+'_genomic.fna -dbtype nucl\n'
        out.write(cmd+'\n')
    out.close()
    os.system('bash index.sh')
    os.system('rm index.sh')

def store_positions(seq_list_forward,seq_list_reverse,organism_code,pam,non_pam_motif_length):
    '''
    Same function as construct_seq_dict except the dictionnary will not be like value=list of coordinates. It add the information about organism with dictionnary values like : dictionnary with key=organism and value=list of coordinates in this organism
    '''
    fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=genome_seqrecord.seq
    seq_dict={}
    for indice in seq_list_forward:
        end=indice+len(pam)+non_pam_motif_length
        seq=genome_seq[indice:end].reverse_complement()
        seq=str(seq)
        if seq not in seq_dict:
            seq_dict[seq]=[]
        seq_dict[seq].append('+('+str(indice+1)+','+str(end)+')')

    for indice in seq_list_reverse:
        end=indice+len(pam)+non_pam_motif_length
        seq=genome_seq[indice:end]
        seq=str(seq)
        if seq not in seq_dict:
            seq_dict[seq]=[]
        seq_dict[seq].append('-('+str(indice+1)+','+str(end)+')')

    pickle.dump(seq_dict,open("reference_genomes/pre_calculate/"+organism_code+"_dicpos.pic","wb"))


def build_expression(seq):
    result = ''
    iupac_code={'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]', 'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}
    for c in seq:
        if c in iupac_code:
            result = result + iupac_code[c]
        else:
            result = result + c
    return result

def find_pam(seq,motif):
    list_seq=[]
    reg_exp = build_expression(motif)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices

def find_sgrna(organism_code,pam,non_pam_motif_length):
    fasta_file='reference_genomes/fasta/' + organism_code +'_genomic.fna'
    genome_seqrecord=next(SeqIO.parse(fasta_file, 'fasta'))
    genome_seq=str(genome_seqrecord.seq)
    sgrna=''
    for i in range(non_pam_motif_length):
        sgrna+='N'
    sgrna+=pam
    seq_list_forward=find_pam(genome_seq,reverse_complement(sgrna))
    seq_list_reverse=find_pam(genome_seq,sgrna)

    return seq_list_forward,seq_list_reverse

def pre_calculate(list_ref):
    print('PRECALCULATE')
    count=0
    for ref in list_ref:
        count+=1
        print(count)
        if ref+'_dicpos.pic' not in os.listdir('reference_genomes/pre_calculate'):
            list_forward,list_reverse=find_sgrna(ref,'NGG',20)
            store_positions(list_forward,list_reverse,ref,'NGG',20)

def create_lineage_objects(dic_tax):
    ncbi=NCBITaxa()
    dic_lineage={}
    count=0
    for ref in dic_tax:
        lineage_object=Lineage()
        tax_ref=dic_tax[ref]
        lineage=ncbi.get_lineage(tax_ref)
        names=ncbi.get_taxid_translator(lineage)
        ranks=ncbi.get_rank(lineage)
        for i in ranks:
            if ranks[i]=='species':
                lineage_object.species=names[i]
            elif ranks[i]=='genus':
                lineage_object.genus=names[i]
            elif ranks[i]=='family':
                lineage_object.family=names[i]
            elif ranks[i]=='order':
                lineage_object.order=names[i]
            elif ranks[i]=='class':
                lineage_object.classe=names[i]
            elif ranks[i]=='phylum':
                lineage_object.phylum=names[i]
        dic_lineage[ref]=(lineage_object,count)
        count+=1
    return dic_lineage

def distance_dic(dic_lineage):
    dic={}
    for ref1 in dic_lineage:
        dic[ref1]={}
        for ref2 in dic_lineage:
            if ref1==ref2:
                dic[ref1][ref2]=0
            elif dic_lineage[ref1][0].species==dic_lineage[ref2][0].species:
                dic[ref1][ref2]=1
            elif dic_lineage[ref1][0].genus==dic_lineage[ref2][0].genus:
                dic[ref1][ref2]=2
            elif dic_lineage[ref1][0].family==dic_lineage[ref2][0].family:
                dic[ref1][ref2]=3
            elif dic_lineage[ref1][0].order==dic_lineage[ref2][0].order:
                dic[ref1][ref2]=4
            elif dic_lineage[ref1][0].classe==dic_lineage[ref2][0].classe:
                dic[ref1][ref2]=5
            elif dic_lineage[ref1][0].phylum==dic_lineage[ref2][0].phylum:
                dic[ref1][ref2]=6
            else:
                dic[ref1][ref2]=7
    return dic

def distance_matrix(dic_taxid):
    print('DISTANCE MATRIX')
    dic_lineage=create_lineage_objects(dic_taxid)
    dist_dic=distance_dic(dic_lineage)
    json.dump(dist_dic, open( "reference_genomes/distance_dic.json", "w" ),indent=4 )

def json_tree():
    print('JSON TREE')
    os.system('python3 scripts/tax2json.py')

def genome_file_for_list():
    print('GENOME FILE FOR LIST')
    dic=json.load(open('reference_genomes/genome_ref_taxid.json','r'))
    list_genomes=list(dic.keys())
    list_genomes=sorted(list_genomes)
    out=open('scripts/static/sortedgenomes.txt','w')
    for i in list_genomes:
        name=' '.join(i.split(' ')[:-1])
        out.write('<OPTION>'+name+'\n')
    out.close()

def compress():
    out=open('compress.sh','w')
    for i in os.listdir('reference_genomes/fasta/'):
        if i != '.DS_Store':
            if i.split('.')[-1]!='gz':
                out.write('tar -zcf reference_genomes/fasta/'+i+'.tar.gz reference_genomes/fasta/'+i+'\n')
                out.write('rm -r reference_genomes/fasta/'+i+'\n')
    for j in os.listdir('reference_genomes/index2/'):
        if j != '.DS_Store':
            if j.split('.')[-1]!='gz':
                out.write('tar -zcf reference_genomes/index2/'+j+'.tar.gz reference_genomes/index2/'+j+'\n')
                out.write('rm -r reference_genomes/index2/'+j+'\n')

    out.close()
    os.system('bash compress.sh')
    os.system('rm compress.sh')


if __name__ == '__main__':
    param = args_gestion()
    fname(param)
    #setup()
    # dic_taxid=dic_download(param.file)
    #index_bowtie_blast(dic_taxid)
    #distance_matrix(dic_taxid)
    # json_tree()
    #genome_file_for_list()
    #compress()

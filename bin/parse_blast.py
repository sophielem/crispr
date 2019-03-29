#!/usr/bin/env python3
"""
Definition
"""

import argparse
import pickle
import re
import xml.etree.ElementTree as ET


class BlastHit(object):
    """docstring for BlastHit."""
    def __init__(self, start, end, len):
        self.start = start
        self.end = end
        self.len = len

    def __init__(self, hsp):
        self.start = int(hsp.find("Hsp_hit-from").text) - 1
        self.end = int(hsp.find("Hsp_hit-to").text) - 1
        self.len = int(hsp.find("Hsp_align-len").text)

    def __repr__(self):
        return "\nStart : {}\nEnd : {}\nLength aln : {}".format(self.start, self.end, self.len)


class BlastReport(object):
    """docstring for BlastHit."""
    def __init__(self, f_name, id_min, genomes_in):
        tree = ET.parse(f_name)
        self.root = tree.getroot()
        self.homolog_gene = self._parse_(id_min, genomes_in)

    def __getitem__(self, k):
        if k in self.homolog_gene:
            return self.homolog_gene[k]
        return None

    def _parse_(self, id_min, genomes_in):
        data = {}
        len_query = int(self.root.find("BlastOutput_query-len").text)
        for hit in self.root.iter(tag="Hit"):
            org_name = hit.find("Hit_def").text
            ref = org_name.split(" ")[0]
            org_name = org_name.replace(ref, "").strip()
            subdata = {ref : []}
            real_name = "".join(list(filter(lambda x: check_in_gi(x, org_name), genomes_in)))
            if real_name:
                for hsp in hit.iter(tag="Hsp"):
                    if (int(hsp.find("Hsp_identity").text)/ len_query) * 100 > id_min:
                        subdata[ref].append(BlastHit(hsp))
                if real_name not in data:
                    data[real_name] = {}
                data[real_name].update(subdata)
        return data

    def is_hit(self):
        msg = self.root.find("./BlastOutput_iterations/Iteration/Iteration_message")
        if not msg:
            return True
        else:
            msg = msg.text
            print(msg)
            return False

    def org_names(self):
        return list(self.homolog_gene.keys())

    def ref_names(self, org_name):
        if org_name in self.org_names():
            return list(self[org_name].keys())
        return None

    def gene_coords(self, org_name, ref_name):
        if ref_name in self.ref_names(org_name):
            return list(self[org_name][ref_name])
        return None


def check_in_gi(gi_name, blast_name):
    # Remove the GCF code to the end
    gi_name = " ".join(gi_name.split(" ")[:-1]).strip()
    return re.search(gi_name + "[ ,$]", blast_name)


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Parse the output of blastn")
    parser.add_argument("-blast", metavar="<str>",
                        help="The path to the blastn output file",
                        required=True)
    parser.add_argument('-ip', type=int,
                        help='identity percentage min for the research of homologous genes using blastn (default:70)',
                        nargs="?", const=70)
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    OUTPUT_BLAST = BlastReport(PARAM.blast, PARAM.ip, GENOMES_IN)
    if OUTPUT_BLAST.is_hit():
        pickle.dump(OUTPUT_BLAST, open("output_blast.p", "wb"), protocol=3)
    else:
        print("Progam terminated&No homologous gene found")

# OUTPUT_BLAST = BlastReport("data/xml.xml", 70, ["Bacillus pseudofirmus OF48", "Bacillus pseudofirmus OF4", "NOTIN"])
# print(OUTPUT_BLAST.homolog_gene)

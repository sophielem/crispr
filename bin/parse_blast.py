#!/usr/bin/env python3
"""
Definition
"""

import argparse
import pickle
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
    def __init__(self, f_name):
        tree = ET.parse(f_name)
        self.root = tree.getroot()
        self.homolog_gene = self._parse_()

    def _parse_(self):
        data = {}
        for hit in self.root.iter(tag="Hit"):
            org_name = hit.find("Hit_def").text
            ref = org_name.split(" ")[0]
            org_name = org_name.replace(ref, "").strip()
            subdata = {ref : []}

            if org_name not in data:
                data[org_name] = {}
            for hsp in hit.iter(tag="Hsp"):
                subdata[ref].append(BlastHit(hsp))
            data[org_name].update(subdata)
        return data

    def __getitem__(self, k):
        if k in self.homolog_gene:
            return self.homolog_gene[k]
        return None

    def is_hit(self):
        msg = self.root.find("./BlastOutput_iterations/Iteration/Iteration_message")
        if not msg:
            return True
        else:
            msg = msg.text
            print(msg)
            return False

    def coord_gene_iter(self, k, c):
        res = self[k][c]
        for i in res:
            yield i

    def ref_iter(self, k):
        for ref in self[k]:
            yield ref


    def homol_gene_iter(self):
        for org in self.homolog_gene:
            yield org


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Parse the output of blastn")
    parser.add_argument("-blast", metavar="<str>",
                        help="The path to the blastn output file",
                        required=True)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    PARAM = args_gestion()
    OUTPUT_BLAST = BlastReport(PARAM.blast)
    if OUTPUT_BLAST.is_hit():
        pickle.dump(OUTPUT_BLAST, open("output_blast.p", "wb"), protocol=3)
    # for org in OUTPUT_BLAST.homol_gene_iter():
    #     for ref in OUTPUT_BLAST.ref_iter(org):
    #         for coord in OUTPUT_BLAST.coord_gene_iter(org, ref):
    #             print(coord.start)

#!/usr/bin/env python3

import argparse
import word_detect
import wordIntegerIndexing as decoding


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Create pickleand index metafile")
    parser.add_argument("-file", metavar="<str>",
                        help="The fasta file to parse",
                        required=True)
    parser.add_argument("-out", metavar="<str>",
                        help="The output path without the extension",
                        required=True)
    parser.add_argument("-rfg", metavar="<str>",
                        help="The path to the database for index and pickle file",
                        nargs="?", const="")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    PARAM = args_gestion()
    PATH = PARAM.rfg + "/genome_pickle/" if PARAM.rfg else ""
    DIC_PICKLE = word_detect.construct_in(PARAM.file, PATH + PARAM.out + ".p")
    if DIC_PICKLE:
        PATH = PARAM.rfg + "/genome_index/" if PARAM.rfg else ""
        decoding.indexAndOccurencePickle(PARAM.out + ".p", PATH + PARAM.out + ".index")
    else:
        print("Program terminated&No sgRNA sequences on the query")

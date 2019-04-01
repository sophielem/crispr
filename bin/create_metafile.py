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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    PARAM = args_gestion()
    word_detect.construct_in(PARAM.out + ".p", PARAM.file)
    decoding.indexPickle(PARAM.out + ".p", PARAM.out + ".index")

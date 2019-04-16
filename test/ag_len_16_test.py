#!/usr/bin/env python3
"""
Test for All genome script with length of sgRNA of 16. Use the output of setCompare from the folder
test/data/ag_len_16
The total number of sgRNA is 137
"""

import os
import sys
import json
import difflib


if __name__ == '__main__':
    gi = "Buchnera aphidicola (Cinara tujafilina) GCF_000217635.1&Aliivibrio wodanis GCF_000953695.1&Photobacterium gaetbulicola Gung47 GCF_000940995.1"
        os.system("python -u bin/post_processing.py -f test/data/ag_len_16/set_index.txt -sl 16 -pam \"NGG\" -gi \"" + gi + "\" -gni \"\" -r " + sys.argv[1] + "  -c 2000 --no-proxy")

    res = json.load(open("test/results.json", "r"))
    # Check if the total number of sgRNA is correct
    if len(res) != 100:
        sys.exit("Problem with length 16 test")

    # Check if found sequences are corect
    ref = json.load(open("test/data/ag_len_16/results.json", "r"))

    diff = difflib.ndiff(res.readlines(), ref.readlines())
    delta = ''.join(x[2:] for x in diff if x.startswith('- ') or x.startswith("+ "))

    if delta:
        sys.exit("Problem with length 16 test")

#!/usr/bin/env python3

import os
import sys
import json
import difflib


if __name__ == '__main__':
    gi = "Buchnera aphidicola (Cinara tujafilina) GCF_000217635.1&Aliivibrio wodanis GCF_000953695.1"
    gni = "Grimontia hollisae GCF_001558255.1"
    os.system("python -u bin/post_processing.py -f test/data/ag_exclude/set_index.txt -sl 20 -pam \"NGG\" -gi \"" + gi + "\" -gni \"" + gni + "\" -r " + sys.argv[1] + "  -c 2000 --no-proxy")
    res = json.load(open("test/results.json", "r"))

    if len(res) != 12:
        sys.exit("Problem with exclude test")

    ref = json.load(open("test/data/ag_exclude/results.json", "r"))

    diff = difflib.ndiff(res.readlines(), ref.readlines())
    delta = ''.join(x[2:] for x in diff if x.startswith('- ') or x.startswith("+ "))

    if delta:
        sys.exit("Problem with exclude test")

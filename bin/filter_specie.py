"""
    Filter out passed specie which are not found in ref dictionnary
    Usage:
        filter_specie.py --ref=<specieDictionary.json> --query=<specie1&specie2>

    Options:
        -h --help                               Show this screen

"""

import json
from docopt import docopt



if __name__ == "__main__":
    ARGUMENTS = docopt(__doc__, version='filter_specie 1.0')

    #if ARGUMENTS['--ref']:
    #    print (ARGUMENTS['--ref'])


    #if ARGUMENTS['--query']:
    #    print (ARGUMENTS['--query'])


    with open(ARGUMENTS['--ref'], "r") as fp:
        data = json.load(fp)
        #print(data)

        print ( '&'.join([ sp for sp in ARGUMENTS['--query'].split('&') if sp in data ]) ) 
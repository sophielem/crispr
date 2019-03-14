"""
    Translate a dictionary of constant-length CRISPR motifs into an ordered set of integer, using base4 encoding.
    Usage:
        wordIntegerIndexing.py <pickledDictionary> [--out=<outFile>]

    Options:
        -h --help                               Show this screen
        -o <outFile>, --out <outFile>           Name of the output file [default : ./pickledDictionary.index]

"""

from docopt import docopt

import os, pickle
import math

def indexPickle(filePath, targetFile):
    pData = pickle.load( open( filePath, "rb" ) )
    wordList = list(pData.keys())
    data = sorted([ weightWord(w, ["A", "T", "C", "G"], len(wordList[0])) for w in wordList ])

    with open(targetFile, "w") as fp:
        fp.write(str(len(data)) + "\n")
        for l in data:
            fp.write(str(l) + "\n")

    return len(data)


def weightWord(word, alphabet, length=None) :
    rank = 0
    if length:
        if length != len(word):
            raise ValueError ("Irregular word length " + str(len(word)) + " (expected " + str(length) + ")")
    for n,c in enumerate(word[::-1]):
        wei = alphabet.index(c)
        base = len(alphabet)
        rank += wei * pow(base,n)
    return rank


def decode(rank, alphabet, length=20):
    """
    Definition
    """
    word = ""
    base = len(alphabet)
    for n in range(length - 1, -1, -1):
        index = math.trunc(rank / pow(base, n))
        word += alphabet[index]
        rank = rank % pow(base, n)
    return word


if __name__ == "__main__":
    arguments = docopt(__doc__, version='wordIntegerIndexing 1.0')

    targetFile = '.'.join(os.path.basename(arguments['<pickledDictionary>']).split('.')[0:-1]) + '.index'
    if arguments['--out']:
        targetFile = arguments['--out']

    total = indexPickle(arguments['<pickledDictionary>'], targetFile)
    print ("Successfully indexed", total, "words\nfrom:", arguments['<pickledDictionary>'], "\ninto:", targetFile)

    rank = 45
    mot = decode(rank, ["A", "T", "C", "G"], 23)
    print(mot)

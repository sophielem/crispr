# CRISPR : search common sgrna sequences excluding genomes
This script is a part of a pipeline which allows to find common sgrna from a set of genomes of bacteria which cannot match with a set of other genomes.

With this, it is possible to select an heterogeneous population of bacteria and to kill a second population.

## Implemented functions
This python script decode the plain text given in argument with the alphabet ["A", "T", "C", "G"].<br>
Then, it does a research in the crispr database with the decoded sequences and retrieves all organism containing the sequence for each sequence with their coordinates in each organism. Only genomes selected by user are conserved and results are displayed by a *json* file containing the first 100 hits and a *text* file containing the first 10,000 hits.

##### Format of the *json* file :
```
{'sequence' : word, 'occurences' :
                                    {'org' : genome, 'all_ref' :
                                                                {'ref' : ref, 'coords' : [coordinates]
                                                               }
                                    }
}
```

##### Format of the *text* file :
```
#ALL GENOMES
#Genomes included : *list of genomes*  ; Genomes excluded : *list of genomes*
#Parameters: pam: *PAM* ; sgrna size: *size*
        genome_in_1     genome_in_2     genomes_in_3...
word    ref : coordinates   ref : coordinates   ref : coordinates...
.
.
.
```


## Dependencies
To execute this script, you need few dependencies :
* docopt&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;== 0.6.2
* ete3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;== 3.1.1
* biopython == 1.66
* tqdm (progress bar)

Then, you need the pycouchDB package which can not be installed by *pip* so the source code is in the *bin* folder.

## Example
Explain how to use the script
```
python post_processing.py -rfg ../reference_genomes -sl 20 -pam "NGG"\
-gi "genome1&genome2&genome3" -gni "genome4&genome5" -file coding.txt
```

# Insert new genomes

## A single genome
### Create pickle and index files

### Add sgrna sequences in the database


## From a node
Need a topology tree in *json* format and the name of the node. Then, the algorithm search after the subtree and from this subtree retrieve all leaves.
Then with the list, we retrieve GCF, ASM and TAXON id from the *genome_ref_taxid.json* file and create a configure file with it in a folder. Then, associated fasta file are copied in the folder and the process to add them in the database is launched.

```
def search_subtree(tree, value):
    # Node name is the searched node
    if tree["text"] == value: return tree
    # No children, so return None
    if "children" not in list(tree.keys()): return None
    # Traverse the tree via children
    for subref in tree["children"]:
        res = search_subtree(subref, value)
        if res != None:
            return res

def search_leaves(tree, list_child):
    # The node has not children, so it is a leave else traverse children
    if "children" not in list(tree.keys()):
        list_child.append(tree["text"])
        return list_child
    # Traverse the tree via children
    for subref in tree["children"]:
        list_child = search_leaves(subref, list_child)
    return list_child
```

## Date
March 15 2019
# Indexing ATCG words

### Usage

###### Specifying any output

```sh
python wordIntegerIndexing.py data/example.pickle --out /any/folder/toto.index
```

will produce `/any/folder/toto.index`

###### Otherwise output file name is guessed from input and wrote locally

```sh
python wordIntegerIndexing.py data/example.pickle
```

will produce `./example.index`

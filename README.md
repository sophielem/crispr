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

Then, you need the pycouchDB package which can not be installed by *pip* so the source code is in the *bin* folder.

## Example
Explain how to use the script
```
python post_processing.py -rfg ../reference_genomes -sl 20 -pam "NGG"\
-gi "genome1&genome2&genome3" -gni "genome4&genome5" -file coding.txt
```

## Date
March 15 2019

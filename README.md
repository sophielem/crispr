# CRISPR : search common sgrna sequences excluding genomes
This script is a part of a pipeline which allows to find common sgrna from a set of genomes of bacteria which cannot match with a set of other genomes.

With this, it is possible to select an heterogeneous population of bacteria and to kill a second population.

## Implemented functions
This python script decode the plain text given in argument with the alphabet ["A", "T", "C", "G"].<br>
Then, it does a research in the crispr database with the decoded sequences and retrieves all organism containing the sequence for each sequence with their coordinates in each organism. Only genomes selected by user are conserved and results are displayed by a *json* file containing the first 100 hits and a *text* file containing the first 10,000 hits.

##### Format of the *json* file :
```json
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

## Example normal mode
Explain how to use the script
```sh
usage: post_processing.py [-h] [-sl <int>] -pam <str> -gi <str> -gni <str> -r
                          <str> -c <int> -f <str> [--no-proxy] [-nb_top <int>]

Post-processing results

optional arguments:
  -h, --help     show this help message and exit
  -sl <int>      The length of the sgrna, excluding pam
  -pam <str>     The pam motif
  -gi <str>      The organisms to search inclusion in.
  -gni <str>     The organisms to search exclusion from
  -r <str>       The end point
  -c <int>       The length of the slice for the request
  -f <str>       The file index
  --no-proxy
  -nb_top <int>  The top hits to download
```

```sh
python bin/post_processing.py -sl 20 -pam "NGG"\
-gi "genome1&genome2&genome3" -gni "genome4&genome5"\
-r "http://localhost:1234" -c 1000 -f "data/example_outputC.txt"
```

## Example Specific gene
Retrieve sgRNA on a specific gene given by user. Retrieve result from SetCompare script and
create files for displaying on webservice.

```sh
usage: bin/specific_gene.py [-h] -gi <str> -gni <str> [-sl [<int>]]
                            [-pam [<str>]] -f <str> -blast <str>
                            [-nb_top <int>] [--no-proxy] -r <str> -c <int>

Specific gene program.

optional arguments:
  -h, --help     show this help message and exit
  -gi <str>      The organisms to search inclusion in.
  -gni <str>     The organisms to search exclusion from
  -sl [<int>]    The length of the sgrna, excluding pam
  -pam [<str>]   The pam motif
  -f <str>       The file index
  -blast <str>   The blast output file
  -nb_top <int>  The top hits to download
  --no-proxy
  -r <str>       The end point
  -c <int>       The length of the slice for the request
  ```

## Example Create Metafile (*pickle* and *index*)

```sh
usage: create_metafile.py [-h] -file <str> -out <str>

Create pickleand index metafile

optional arguments:
  -h, --help   show this help message and exit
  -file <str>  The fasta file to parse
  -out <str>   The output path without the extension
 ```
 The output path is the path to the output file without the extension.

## Example parse blastn output
Make db : makeblastdb -in file.fasta -title "Database with all genomes" -dbtype nucl

```sh
usage: parse_blast.py [-h] -blast <str> [-ip [IP]] -gi <str>

Parse the output of blastn

optional arguments:
  -h, --help    show this help message and exit
  -blast <str>  The path to the blastn output file
  -ip [IP]      identity percentage min for the research of homologous genes
                using blastn (default:70)
  -gi <str>     The organisms to search inclusion in
 ```

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

## Date
April 01 2019

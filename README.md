# CRISPR : search common sgrna sequences to target bacterial strains
This script is a part of a pipeline which allows to find common sgrna from a set of genomes of bacteria which cannot match with a set of other genomes.

With this, it is possible to select an heterogeneous population of bacteria and to kill a second population.
![alt text](https://github.com/sophielem/crispr/blob/dev_add_genome/schema_crispr_service.png)

## Dependencies
To execute this script, you need few dependencies :
* docopt&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;== 0.6.2
* ete3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;== 3.1.1
* biopython == 1.66

You also need blastn.

## Results files
Only genomes selected by user are conserved and results are displayed by a *json* file containing the first 100 hits and a *text* file containing the first 1,000 hits.
Another *json* file contains results sorted by organism name and not by sgRNAs.
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




## Usage
### Normal mode

```sh
usage: post_processing.py [-h] [-sl <int>] -pam <str> -gi <str> -gni <str> -r
                          <str> -c <int> -f <str> [--no-proxy] [-n
b_top <int>]

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

### Specific gene
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
### Check taxonomy
Check if a taxonomy ID given exists in the NCBI database. Then, check if this taxonomy ID is
present in a Taxonomy database than you can request with the package pyCouch and check
if the GCF given for this taxonomy ID already exists in the Taxonomy database.

### Create file for taxon database
Create pickle file to inset into Taxonomy database
It can create from a GCF and Taxonomy ID given or from the json file genome_ref_taxid.json file

OUTUPUT:
taxon_dt={"1234": {"GCF": ["GCF_111", "GCF_112", "GCF_113"], "date": "19-04-2019", "User": "Queen"},
          "2345": {"GCF": ["GCF_211", "GCF_212", "GCF_213"], "date": "19-04-2019", "User": "Queen"}}

### Update tree
Load the MaxiTree object from the taxon_tree_db and insert a new member.
Then, create a pickle file to insert it into the database.

### Maxi Tree object
Class MaxiTree object. This tree contains: a tree object from the package ete3 with the name
of the organism, the current GCF and its taxonomy Id if it's present into the CRISPR database. This taxon ID is
referenced in the name separate by a ':' and in a node's feature 'taxon'.
The name of a node is composed like this:
    name (without ' and /) GCF_ID : taxon_ID

Several constructor:
    1. from the database
    2. from the genomre_ref_taxid.json file
    3. from a MaxiTree pickle file

Several methods are implemented to insert a new node, to write a json file of the tree
full (with taxonID) or not, to get the json string full or not...

The main function allow to create this Tree from the genomre_ref_taxid.json file and to create a
pickle MaxiTree file to insert into the taxon_tree_db



### Create Metafile (*pickle* and *index*)
```sh
usage: create_metafile.py [-h] -file <str> [-out <str>] [-taxid <str>]
                          [-gcf <str>] [-rfg [<str>]] [-plasmid]

Create pickle and index metafile

optional arguments:
  -h, --help    show this help message and exit
  -file <str>   The fasta file to parse
  -out <str>    Path for the output file without extension
  -taxid <str>  Taxonomy ID
  -gcf <str>    GCF ID
  -rfg [<str>]  The path to the database for index and pickle file
  -plasmid      If present, indicates the genome is a plasmid
```

#### Word detection
Detect sgRNAs by regex and dump results in a pickle file.

#### Indexing ATCG words

```sh
python wordIntegerIndexing.py data/example.pickle --out /any/folder/toto.index
```

will produce `/any/folder/toto.index`

###### Otherwise output file name is guessed from input and wrote locally

```sh
python wordIntegerIndexing.py data/example.pickle
```

will produce `./example.index`



### Parse blastn output
If no database for blast software, use this command to create it:
```sh
makeblastdb -in file.fasta -title "Database with all genomes" -dbtype nucl
```

Parse the output of a Blastn on the reference database. Create a BlastReport object which
contains only hits for organism which are in a given list, here the included genomes list, and
hits which have a percentage identity superior to a given percentage identity, by default 70.
The BlastReport contains a dictionary of organism with references which contains BlastHit object
This object contains the coordinates of the hit and its length.
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

## Date
June 21 2019

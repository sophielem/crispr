# CRISPR workflow
![Worklofw script](https://github.com/sophielem/crispr/blob/dev_add_genome/doc/workflow_script.png)
## crispr workflow

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

## crispr workflow specific
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

## Add genome
### Check taxonomy
Check if a taxonomy ID given exists in the NCBI database. Then, check if this taxonomy ID is
present in a Taxonomy database than you can request with the package pyCouch and check
if the GCF given for this taxonomy ID already exists in the Taxonomy database.

### Create file for taxon database
Create pickle file to inset into Taxonomy database
It can create from a GCF and Taxonomy ID given or from the json file genome_ref_taxid.json file



### Update tree
Load the MaxiTree object from the taxon_tree_db and insert a new member.
Then, create a pickle file to insert it into the database.


## Date
June 24 2019

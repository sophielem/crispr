#!/bin/bash

error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
}


parse_logFile () {
    if grep "Program terminated" $1 > /dev/null;
    then
        perl -ne 'if ($_ =~ /Program terminated/){
            @error_split=split(/&/);
            $msg = $error_split[1];
            $msg =~ s/\n$//;
            print "{\"emptySearch\" :  \"$msg\" }";
        }' $1 > ./fail.log;
        cat ./fail.log;
        false
    else
         true
    fi
}


if [ "$URL_TAXON" != "" ]; then
    error_json
elif [ "$CRISPR_TOOL_SCRIPT_PATH" = "" ]; then
    error_json

elif  [ "$URL_CRISPR" = "" ]; then
    error_json
elif [ "$DB_TAXON_NAME" = "" ]; then
    error_json
else

    ### CREATE PICKLE AND INDEX METAFILE ###
    echo -u $CRISPR_TOOL_SCRIPT_PATH/python create_metafile.py -file $FASTA_FILE -taxid $TAXID -gcf $GCF -rfg $rfg 2> ./create_meta.err 1> ./create_meta.log
    echo cp $FASTA_FILE $rfg"/genome_fasta/"$TAXID"_"$GCF".fna"
    # Check if any sgRNA has been found in this genome
    parse_logFile ./create_meta.log
    PRG_TERMINATED=$?

    ### ADD INTO DATABASE AND UPDATE JSON_TREE ###
    if [ $PRG_TERMINATED = 0 ]; then
        NAME_FILE=`cat ./create_metafile.log`
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/couchBuild.py --url $URL_CRISPR --map $MAP_FILE --data $rfg"/genome_pickle/"$NAME_FILE".p"
        if [ -f error_add_db.err ]; then
            echo "{\"emptySearch\": \"There is a problem during the insertion into CRISPR database. Contact us \"}" > fail.log
            cat fail.log
        else
            echo python -u $CRISPR_TOOL_SCRIPT_PATH/create_file_taxondb.py single -user -gcf $GCF -taxid $TAXID -url $URL_TAXON"/"$DB_TAXON_NAME 2> ./taxondb_file.err 1> ./taxondb_file.log
            parse_logFile ./taxondb_file.log
            PRG_TERMINATED=$?
            if [ $PRG_TERMINATED = 0 ]; then
              ## Add into taxon database ##
              echo python -u $CRISPR_TOOL_SCRIPT_PATH/couchBuild $DB_TAXON_NAME --url $URL_TAXON --data "./taxonDB_data/taxon_dt.p"
              ## Update the Json_tree ##
              echo python -u $CRISPR_TOOL_SCRIPT_PATH/update_tree.py -url $URL_TREE"/"$DB_TREE_NAME -taxid $TAXID  2> ./update_tree.err 1> ./update_tree.log
              # Check if no problem to connect to the database (taxon_db and taxon_tree_db), to access to GCF
              parse_logFile ./update_tree.log
              PRG_TERMINATED=$?
            fi
            if [ $PRG_TERMINATED = 0 ]; then
              echo python -u $CRISPR_TOOL_SCRIPT_PATH/couchBuild $DB_TREE_NAME --url $URL_TREE --data "./treeDB_data/maxi_tree.p"
            fi
        fi
    fi
fi

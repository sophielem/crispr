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
elif [ "$DB_NAME" = "" ]; then
    error_json
else

    ### CHECK IF TAXON ID ALREADY PRESENT IN NCBI AND OUR TAXON DATABASE ###
    echo python check_taxonomy -taxid $TAXID -gcf $GCF -out $INFO_FILE -url $URL_TAXON -dbName $DB_NAME 2> ./check_taxon.err 1> ./check_taxon.log
    if grep "Be careful" ./check_taxon.log > /dev/null;
    then
        perl -ne 'if ($_ =~ /Be careful/){
            @error_split=split(/&/);
            $msg = $error_split[1];
            $msg =~ s/\n$//;
            print "{\"BeCareful\" :  \"$msg\" }";
        }' ./check_taxon.log > ./careful.log;
        cat ./careful.log;
    fi
    # Check if the taxonomy id given is present into NCBI database
    parse_logFile ./check_taxon.log
    PRG_TERMINATED=$?

    ### CREATE PICKLE AND INDEX METAFILE ###
    if [ $PRG_TERMINATED = 0 ]; then
        NAME_FILE=$(echo $ORG_NAME $GCF)

        echo python create_metafile.py -file $FASTA_FILE -out $NAME_FILE -rfg $rfg 2> ./create_meta.err 1> ./create_meta.log
        # Check if any sgRNA has been found in this genome
        parse_logFile ./create_meta.log
        PRG_TERMINATED=$?
    fi

    ### ADD INTO DATABASE AND UPDATE JSON_TREE ###
    if [ $PRG_TERMINATED = 0 ]; then
        echo python couchBuild.py --url $URL_CRISPR --map $MAP_FILE --data $NAME_FILE".p"
        if [ -f error_add_db.err ]; then
            echo "{\"emptySearch\": \"There is a problem during the insertion into CRISPR database. Contact us \"}" > fail.log
            cat fail.log
        else
            ## Add into taxon database ##
            echo python couchBuild $DB_NAME --url $URL_TAXON --data $INFO_FILE
            ## Update the Json_tree ##
            echo python update_tree.py 2> ./update_tree.err 1> ./update_tree.log

            if grep "Program terminated" ./update_tree.log > /dev/null;
            then
            perl -ne 'if ($_ =~ /Program terminated/){
                @error_split=split(/&/);
                $msg = $error_split[1];
                $msg =~ s/\n$//;
                print "{\"emptySearch\" :  \"$msg\" }";
            }' ./update_tree.log > ./fail.log;
            cat ./fail.log;
            else
                perl -ne 'if ($_ =~ /SUCCESS/){
                    @success_split=split(/&/);
                    $msg = $success_split[1];
                    $msg =~ s/\n$//;
                    print "{\"Success\" :  \"$msg\" }";
                };
                ' ./update_tree.log > ./success.log;
                cat ./success.log
            fi
        fi
    fi
fi

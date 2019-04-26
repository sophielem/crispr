#!/bin/bash

echo python check_taxonomy -taxid $TAXID -gcf $GCF -out $INFO_FILE -url $URL_TAXON -dbName $DB_TAXON_NAME 2> ./check_taxon.err 1> ./check_taxon.log
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

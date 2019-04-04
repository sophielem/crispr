#!/bin/bash

error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
}

if [ "$pam" != "NGG" ]; then
    error_json

#elif [ "$sl" != "20" ]; then
#    error_json

elif [ "$CRISPR_TOOL_SCRIPT_PATH" = "" ]; then
    error_json

elif  [ "$URL_CRISPR" = "" ]; then
    error_json

else

    slFlag=""
    #  shorter word
    if [ "$sl" != "20" ]; then
        ((_sl = $sl + 3))
        slFlag="-d 23 -c ${_sl}"
    fi
    printenv > env.log

    BASE_FOLDER=`pwd`

    #cp -r $BASE_FOLDER/* $WORKDIR
    #cd $WORKDIR
    pwd > pwd.log

    # Create Metafile
    echo python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $queryFasta -out $metafileQuery > sg.cmd
    python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $queryFasta -out $metafileQuery

    # Filter genomes
    gi=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gi")
    gni=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gni")
    echo $gi > f.gi

    # Set Compare
    fileSet="set_index.txt"
    setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet -s $metafileQuery".index" 2>> ./setCompare.err 1> ./setCompare.log

    # Blast N to find homologous genes
    echo blastn -outfmt 5 -query $queryFasta -db $dbBlast > $fileBlast >> sg.cmd
    blastn -outfmt 5 -query $queryFasta -db $dbBlast > $fileBlast

    # Parse the blast output
    echo python -u $CRISPR_TOOL_SCRIPT_PATH/parse_blast.py -blast $fileBlast -gi "$gi" -o $parseBlast >> sg.cmd
    python -u $CRISPR_TOOL_SCRIPT_PATH/parse_blast.py -blast $fileBlast -gi "$gi" -o $parseBlast 2>> ./parse_blast.err 1> ./parse_blast.log

    if grep "Program terminated" ./parse_blast.log > /dev/null;
    then
    perl -ne 'if ($_ =~ /Program terminated/){
        @error_split=split(/&/);
        $msg = $error_split[1];
        $msg =~ s/\n$//;
        print "{\"emptySearch\" :  \"$msg\" }";
    }' ./parse_blast.log > ./fail.log;
    cat ./fail.log;
    else

        # Post-processing with setCompare output and blast output
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/specific_gene.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -c 2000 --no-proxy -blast $parseBlast >> sg.cmd
        python -u $CRISPR_TOOL_SCRIPT_PATH/specific_gene.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -c 2000 --no-proxy -blast $parseBlast 2>> ./specific_gene.err 1> ./specific_gene.log


        if grep "Program terminated" ./specific_gene.log > /dev/null;
        then
        perl -ne 'if ($_ =~ /Program terminated/){
            @error_split=split(/&/);
            $msg = $error_split[1];
            $msg =~ s/\n$//;
            print "{\"emptySearch\" :  \"$msg\" }";
        }' ./specific_gene.log > ./fail.log;
        cat ./fail.log;
        else
            not_in=$(perl -ne 'chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;' ./specific_gene.log);
            number_hits=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 3){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);
            tag=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 2){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);
            echo "$not_in" > ./stuff.log;
            echo "$number_hits" >> ./stuff.log;
            echo "$tag" >> ./stuff.log;
            loc=$(pwd | perl -ne '@tmp = split(/\//, $_); print "$tmp[$#tmp - 1]/$tmp[$#tmp]";');
            echo "{\"out\" : {\"data\" : $(cat ./results.json),  \"not_in\" : \""$not_in"\",  \"number_hits\" : \""$number_hits"\", \"tag\" : \""$loc"\"}}"
        fi
    fi
fi

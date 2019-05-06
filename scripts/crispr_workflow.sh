#!/bin/bash
source activate crispr
source ../setenv
unset HTTPS_PROXY
unset HTTP_PROXY
CRISPR_TOOL_SCRIPT_PATH="/Users/slematre/Documents/crisprs/crispr/bin/"
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

    #echo "curl -X GET $URL_CRISPR/handshake" > handshake.cmd
    #curl -X GET $URL_CRISPR/handshake &> handshake.log
    gi=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gi")
    gni=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gni")

    echo $gi > f.gi
    fileSet="/Users/slematre/Documents/crisprs/set_index.txt"
    # setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet 2>> ./setCompare.err 1> ./setCompare.log

    echo python -u $CRISPR_TOOL_SCRIPT_PATH/post_processing.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -c 2000 --no-proxy > pp.cmd
    python -u $CRISPR_TOOL_SCRIPT_PATH/post_processing.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -c 2000 --no-proxy 2>> ./post_processing.err 1> ./post_processing.log
    #echo "post_processing.py -sl 20 -pam \"NGG\" -gi \"$gi\" -gni \"$gni\" -r \"$URL_CRISPR\"  -c 2000" >> ./cmd.txt
    #python $CRISPR_TOOL_SCRIPT_PATH/post_processing.py -sl 20 -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR" -c 2000 2>> ./post_processing.err 1> ./post_processing.log


    if grep "Program terminated" ./post_processing.log > /dev/null;
    then
    perl -ne 'if ($_ =~ /Program terminated/){
        @error_split=split(/&/);
        $msg = $error_split[1];
        $msg =~ s/\n$//;
        print "{\"emptySearch\" :  \"$msg\" }";
    }' ./post_processing.log > ./fail.log;
    cat ./fail.log;
    else
        not_in=$(perl -ne 'chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;' ./post_processing.log);
        number_hits=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 3){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./post_processing.log);
        tag=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 2){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./post_processing.log);
        echo "$not_in" > ./stuff.log;
        echo "$number_hits" >> ./stuff.log;
        echo "$tag" >> ./stuff.log;
        loc=$(pwd | perl -ne '@tmp = split(/\//, $_); print "$tmp[$#tmp - 1]/$tmp[$#tmp]";');
    	#number_hits=lines[2].strip()
        echo "{\"out\" : {\"data\" : $(cat ./results.json),  \"not_in\" : \""$not_in"\",  \"number_hits\" : \""$number_hits"\", \"tag\" : \""$loc"\"}}"
    fi


    #cp -r $WORKDIR/* $BASE_FOLDER
    #cd $BASE_FOLDER
fi

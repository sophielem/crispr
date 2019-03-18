#!/bin/bash

error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
    exit 1
}

if [ "$pam" != "NGG" ]; then
    error_json
fi

if [ "$sl" != "20" ]; then
    error_json
fi

if [ "$CRISPR_TOOL_SCRIPT_PATH" = "" ]; then
    error_json
fi

if [ "$URL_CRISPR" = "" ]; then
    error_json
fi

echo "setCompare -i \"$gi\" -o \"$gni\" -l $rfg/index/ -e index" > ./cmd.txt
setCompare -i "$gi" -o "$gni" -l $rfg/index/ -e index 2>> ./setCompare.err 1> ./setCompare.log
echo "post_processing.py -sl 20 -pam \"NGG\" -gi \"$gi\" -gni \"$gni\" -r \"$URL_CRISPR\"  -c 2000" >> ./cmd.txt
python $CRISPR_TOOL_SCRIPT_PATH/post_processing.py -sl 20 -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR" -c 2000 2>> ./post_processing.err 1> ./post_processing.log


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

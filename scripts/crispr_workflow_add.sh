#!/bin/bash
# nb_file=(`ls $FILE_PATH | wc -l`)
# echo "$nb_file"
#
# if [[ $nb_file = 1 ]]; then
#     for FASTA_FILE in $FILE_PATH; do
#         command="python $CRISPR_TOOL_SCRIPT_PATH/pre_treatment.py all -rfg $rfg -file $FASTA_FILE -gcf $gcf -asm $asm -taxid $taxid 2>> ./pre_treatment.err 1> ./pre_treatment.log"
#         echo "$command"
#     done
# else
#     for FASTA_FILE in $FILE_PATH; do
#         sbatch --export=rfg="$rfg",FASTA_FILE="$FASTA_FILE",gcf="$gcf",asm="$asm",taxid="$taxid" create_metafile.sbatch
#     done
#
#     for FASTA_FILE in $FILE_PATH; do
#         command="python $CRISPR_TOOL_SCRIPT_PATH/pre_treatment.py add -rfg $rfg -file $FASTA_FILE -gcf $gcf 2>> ./pre_treatment.err 1> ./pre_treatment.log"
#         echo "$command"
#     done
# fi

for FASTA_FILE in $FILE_PATH; do
    command="python $CRISPR_TOOL_SCRIPT_PATH/pre_treatment.py all -rfg $rfg -file $FASTA_FILE -gcf $gcf -asm $asm -taxid $taxid 2> ./pre_treatment.err 1> ./pre_treatment.log"
    echo "$command"
    if grep "Program terminated" ./pre_treatment.log > /dev/null;
    then
    perl -ne 'if ($_ =~ /Program terminated/){
        @error_split=split(/&/);
        $msg = $error_split[1];
        $msg =~ s/\n$//;
        print "{\"emptySearch\" :  \"$msg\" }";
    }' ./pre_treatment.log > ./fail.log;
    cat ./fail.log;
    else
        perl -ne 'if ($_ =~ /SUCCESS/){
            @success_split=split(/&/);
            $msg = $success_split[1];
            $msg =~ s/\n$//;
            print "{\"message\" :  \"$msg\" }";
        };
        ' ./pre_treatment.log > ./success.log;
        cat ./success.log
    fi

done

#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_transdecoderPredict.sh [input FASTA] [script]\n"
exit 0
fi

inputFASTA="$1"
script="$2"


qsub -v inputFASTA="$inputFASTA" "$script"

#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_transdecoderPredict.sh [input FASTA] [pfam Homology] [blastp Homology] [script]\n"
exit 0
fi

inputFASTA="$1"
pfamHomo="$2"
blastpHomo="$3"
script="$4"


qsub -v inputFASTA="$inputFASTA",pfamHomo="$pfamHomo",blastpHomo="$blastpHomo" "$script"

#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_blastp.sh [pep db] [refPtn] [output] [script]\n"
exit 0
fi

pep="$1"
refPtn="$2"
output="$3"
script="$4"


qsub -v pep="$pep",refPtn="$refPtn",output="$output" "$script"

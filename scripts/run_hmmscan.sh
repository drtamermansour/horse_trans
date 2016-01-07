#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_hmmscan.sh [pep db] [refPtn] [output] [script]\n"
exit 0
fi

pep="$1"
refPfam="$2"
output="$3"
script="$4"


qsub -v pep="$pep",refPfam="$refPfam",output="$output" "$script"

#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_blast.sh [input ptn/DNA sequence] [refPtn] [output] [blastp or blastx script]\n"
exit 0
fi

seq="$1"  ## ptn for blastp and DNA for blastx
refPtn="$2"
output="$3"
script="$4"


qsub -v seq="$seq",refPtn="$refPtn",output="$output" "$script"

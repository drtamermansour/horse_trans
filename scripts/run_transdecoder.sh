#!/bin/sh

if [ $# -lt 5 ]
then
  printf "\nUsage run_transdecoder.sh [samples list] [Reference genome] [Swissprot database] [Pfam database] [script]\n"
  exit 0
fi

sample_list="$1"
genome="$2"
refPtn="$3"
refPfam="$4"
script="$5"


while read assembly; do
  cd $assembly
  qsub -v genome="$genome",refPtn="$refPtn",refPfam="$refPfam" "$script"
done < $sample_list

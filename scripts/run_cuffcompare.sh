#!/bin/sh

if [ $# -lt 4 ]
then
  printf "\nUsage run_cuffcompare.sh [samples] [label] [Genes_GTF_file] [script]\n"
exit 0
fi

sample="$1"
label="$2"
Genes_GTF_file="$3"
script="$4"

qsub -v Genes_GTF_file="${Genes_GTF_file}",sample="${sample}",label="${label}" "${script}"


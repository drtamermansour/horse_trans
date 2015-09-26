#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_cufflinks2.sh [sample] [Genes_GTF_file] [sample label] [script]\n"
exit 0
fi

sample="$1"
Genes_GTF_file="$2"
label="$3"
script="$4"

qsub -v Genes_GTF_file=${Genes_GTF_file},label=${label},sample=${sample} "${script}"

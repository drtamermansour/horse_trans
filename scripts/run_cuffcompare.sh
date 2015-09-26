#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_cuffcompare.sh [samples list] [Genes_GTF_file] [script]\n"
exit 0
fi

sample_list="$1"
Genes_GTF_file="$2"
script="$3"

assemblies=()
while read assembly; do
  assemblies+=($assembly/*.gtf)
done < $sample_list

qsub -v Genes_GTF_file="${Genes_GTF_file}",samples="${assemblies[*]}" "${script}"


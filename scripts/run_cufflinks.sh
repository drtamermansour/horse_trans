#!/bin/sh

if [ $# -lt 4 ]
then
    printf "\nUsage run_cufflinks.sh [samples list] [Genes_GTF_file] [Genome] [script]\n"
    exit 0
fi

sample_list="$1"
Genes_GTF_file="$2"
script="$3"

while read f; do
     cd "$f"
     label=$"Sample_"${f#tophat_}
     echo $label
     qsub -v Genes_GTF_file=${Genes_GTF_file},label=${label} "${script}" 
     cd ..
done < $sample_list

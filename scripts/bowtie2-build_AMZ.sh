#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage bowtie2-build_AMZ.sh [Genome file in fasta format] [Genome index base]\n"
exit 0
fi


genome_file="$1"
genome_index_base="$2"

bowtie2-build $genome_file $genome_index_base


#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_bowtie2-build.sh [Genome file in fasta format] [Genome index base]\n"
exit 0
fi


genome_file="$1"
genome_index_base="$2"

module load bowtie2/2.1.0
bowtie2-build $genome_file $genome_index_base

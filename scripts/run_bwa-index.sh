#!/bin/sh

if [ $# -lt 1 ]
then
printf "\nUsage run_bwa-index.sh [Genome file in fasta format]\n"
exit 0
fi


genome_file="$1"

module load bwa/0.7.7.r441
bwa index -a bwtsw $genome_file


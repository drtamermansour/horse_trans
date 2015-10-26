#!/bin/sh

if [ $# -lt 1 ]
then
printf "\nUsage run_gatk-index.sh [Genome file in fasta format]\n"
exit 0
fi


module load picardTools/1.113
module load SAMTools/0.1.19

genome_file="$1"

java -Xmx4g -jar $PICARD/CreateSequenceDictionary.jar R= $genome_file O= ${genome_file%.fa}.dict
samtools faidx $genome_file
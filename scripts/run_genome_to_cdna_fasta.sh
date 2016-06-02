#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_genome_to_cdna_fasta.sh [input GTF] [Reference genome] [output fasta] [script]\n"
exit 0
fi

inputGTF="$1"
genome="$2"
outputFASTA="$3"
script="$4"
decoderUtil=$(dirname "${BASH_SOURCE[0]}")/decoderUtil 


qsub -v inputGTF="$inputGTF",genome="$genome",outputFASTA="$outputFASTA",decoderUtil="$decoderUtil" "$script"

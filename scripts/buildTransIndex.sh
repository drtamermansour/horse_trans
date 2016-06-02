#!/bin/sh
refGTF_file="$1"
transcriptome_index="$2"
Bowtie2_genome_index_base="$3"

module load TopHat2/2.0.14

tophat --GTF "$refGTF_file" --transcriptome-index "$transcriptome_index" "$Bowtie2_genome_index_base"
rm -R tophat_out

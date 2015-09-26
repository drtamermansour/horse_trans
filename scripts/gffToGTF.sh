#!/bin/sh

input_GFF="$1"
output_GTF="$2"


module load ucscUtils/262   # UCSC_kent_commands
gff3ToGenePred $input_GFF ${input_GFF%.*}.gpred
genePredToGtf file ${input_GFF%.*}.gpred $output_GTF


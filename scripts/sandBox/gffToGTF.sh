#!/bin/sh

input_GFF="$1"
output_GTF="$2"
script_path=$(dirname "${BASH_SOURCE[0]}")


#module load ucscUtils/262   # UCSC_kent_commands
$script_path/UCSC_kent_commands/gff3ToGenePred $input_GFF ${input_GFF%.*}.gpred
$script_path/UCSC_kent_commands/genePredToGtf file ${input_GFF%.*}.gpred $output_GTF


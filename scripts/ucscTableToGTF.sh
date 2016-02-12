#!/bin/sh

ucscTable="$1"
output_GTF="$2"


# module load ucscUtils/262
# UCSC_kent_commands
# zcat $ucscTable | cut -f2-16 | genePredToGtf file stdin $output_GTF

script_path=$(dirname "${BASH_SOURCE[0]}")   
zcat $ucscTable | cut -f2-16 | $script_path/UCSC_kent_commands/genePredToGtf file stdin $output_GTF

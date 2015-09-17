#!/bin/sh

ucscTable="$1"
output_GTF="$2"



module load ucscUtils/262   # UCSC_kent_commands
zcat $ucscTable | cut -f2-11 | genePredToGtf file stdin $output_GTF

#zcat $ucscTable | cut -f2-11 | $HOME/bin/UCSC_kent_commands/genePredToGtf file stdin $output_GTF

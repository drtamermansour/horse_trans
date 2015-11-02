#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################
## The map.chain file has the old genome as the target and the new genome as the query.
## The old genome (target) here should be the genome and the new (quary) is transcriptome

$HOME/bin/UCSC_kent_commands/genePredToFakePsl -chromSize=$genome_dir/$UCSCgenome.chrom.sizes file merged.gpred merged.psl merged.cds
## should I use "pslToChain" or "axtChain" ??
$HOME/bin/UCSC_kent_commands/pslToChain merged.psl genomeToTrans_map.chain
$HOME/bin/UCSC_kent_commands/chainSort genomeToTrans_map.chain genomeToTrans_map.sorted.chain



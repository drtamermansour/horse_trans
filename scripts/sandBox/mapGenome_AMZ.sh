#!/bin/sh
genome="$1"

script_path=$(dirname "${BASH_SOURCE[0]}") 
###########################################################################################
## The map.chain file has the old genome as the target and the new genome as the query.
## The old genome (target) here should be the NCBI and the new (quary) is UCSC

mkdir psl
# assembled chromosomes are too large and require very large memory. We can edit the psl files. All required data can be found in the fasta headers of ncbi files
wget --no-directories ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002305.2.assembly.txt
grep -v "^#" GCF_000002305.2.assembly.txt | awk '{ sub("\r$", ""); print }' > GCF_000002305.2.assembly_noCommonts.txt
cat GCF_000002305.2.assembly_noCommonts.txt | awk -F "\t" '$10 != "na"' | awk -F "\t" -v OFS='\t' '{ print $5,$7,$9,"-",$10,$9,"0",$9 }' > chromosome_map.txt ## genBankAcc refSeqAcc(i.e.ncbiName) ncbilength ucscScaf ucscName ucsclenght start end


# predict mapping of unassemblied scaffolds by size
cat GCF_000002305.2.assembly_noCommonts.txt | awk -F "\t" '$10 == "na"' | sort -t$'\t' -k9,9nr > NCBIunAssembled_map.txt
wget --no-directories http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/ctgPos2.txt.gz
gunzip ctgPos2.txt.gz
cat ctgPos2.txt | awk '$3 == "chrUn"' > UCSCunAssembled_map.txt
cat ctgPos2.txt | awk '$3 == "chrUn"' | awk -F "\t" -v OFS='\t' '{ print "chrUn",$4,$4+1400,$1 }' > ctgPos2_bait.txt
#module load BEDTools/2.24.0
bedtools getfasta -fi "$genome_dir/chrUn.fa" -bed ctgPos2_bait.txt -fo UCSC_baits.tab -name -tab
grep -A20 "^>" eca_ref_EquCab2.0_unplaced.fa | grep -v "^--" | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | awk 'BEGIN{RS=">"}{gsub("\n","\t",$0); print $0}' | sed '/^$/d' > ncbi_prey.tab
Rscript -e 'x=read.table("UCSC_baits.tab", sep="\t"); y=read.table("ncbi_prey.tab", sep="\t"); x$V3=tolower(x$V2); y$V3=tolower(y$V2); m=merge(x,y,by.x="V3",by.y="V3",all.x=TRUE,all.y=TRUE); m2=m[,c(2,4)]; z=data.frame(do.call('rbind', strsplit(as.character(m2$V1.y), "|", fixed=TRUE))); m2$ncbiName=z$X4; x2=read.table("NCBIunAssembled_map.txt", sep="\t"); y2=read.table("UCSCunAssembled_map.txt", sep="\t"); m3=merge(x2,m2,by.x="V7",by.y="ncbiName",all.x=TRUE,all.y=TRUE); m4=merge(m3,y2,by.x="V1.x",by.y="V1",all.x=TRUE,all.y=TRUE); m4$V6.y="117461955"; m5=m4[,c(7,2,13,1,14,17,15,16)]; write.table(m5, "NCBItoUCSCunAssembled_bait_map.txt",sep="\t", quote=F, col.names=F, row.names=F);'   ## genBankAcc refSeqAcc(i.e.ncbiName) ncbilength ucscScaf ucscName ucsclenght start end

cat chromosome_map.txt NCBItoUCSCunAssembled_bait_map.txt > NCBItoUCSC_map.txt
$script_path/NCBItoUCSCmapTopsl NCBItoUCSC_map.txt > NCBItoUCSC_map.psl

sed 's/gi|\(.*\)|ref|\(.*\)| \(.*$\)/\2/g' ncbi_genome.fa > ncbi_genome2.fa
# qsub -v genome="$genome" $script_path/run_axtChain.sh
bash $script_path/run_axtChain_AMZ.sh $genome
###$HOME/bin/UCSC_kent_commands/axtChain -linearGap=medium -psl NCBItoUCSC_map.psl -faT ncbi_genome2.fa -faQ $genome NCBItoUCSC_map.chain
$script_path/UCSC_kent_commands/chainSort NCBItoUCSC_map.chain NCBItoUCSC_map.sorted.chain
cat NCBItoUCSC_map.psl | awk -F "\t" -v OFS='\t' '{ print $14,$15 }' > ncbi.chromInfo
cat NCBItoUCSC_map.psl | awk -F "\t" -v OFS='\t' '{ print $10,$11 }' | uniq > ucsc.chromInfo
$script_path/UCSC_kent_commands/chainNet NCBItoUCSC_map.sorted.chain ncbi.chromInfo ucsc.chromInfo NCBItoUCSC_map.net UCSCtoNCBI_map.net
$script_path/UCSC_kent_commands/netChainSubset NCBItoUCSC_map.net NCBItoUCSC_map.sorted.chain NCBItoUCSC.liftOver





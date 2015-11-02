#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################

# assembled chromosomes are too large and require very large memory. We can edit the psl files. All required data can be found in the fasta headers of ncbi files
wget --no-directories ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002305.2.assembly.txt
grep -v "^#" GCF_000002305.2.assembly.txt | awk '{ sub("\r$", ""); print }' > GCF_000002305.2.assembly_noCommonts.txt
cat GCF_000002305.2.assembly_noCommonts.txt | awk -F "\t" '$10 != "na"' > chromosome_map.txt
$script_path/chromosomeMapToNCBIpsl chromosome_map.txt > psl/chromosome_map.psl

# predict mapping of unassemblied scaffolds by size
cat GCF_000002305.2.assembly_noCommonts.txt | awk -F "\t" '$10 == "na"' | sort -t$'\t' -k9,9nr > NCBIunAssembled_map.txt  ## size is col#9
wget --no-directories http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/ctgPos2.txt.gz
gunzip ctgPos2.txt.gz
cat ctgPos2.txt | awk '$3 == "chrUn"' > UCSCunAssembled_map.txt ## size is col#2
Rscript -e 'x=read.table("NCBIunAssembled_map.txt", sep="\t"); y=read.table("UCSCunAssembled_map.txt", sep="\t"); dup=unique(x[duplicated(x$V9),9]); xq=x[!(x$V9 %in%  dup),]; yq=y[!(y$V2 %in% dup),]; m=merge(xq,yq,by.x="V9",by.y="V2",all.x=TRUE,all.y=TRUE); m2=m[,c(6,8,1,11,12,13,14)]; write.table(m2, "NCBItoUCSCunAssembled_map.txt",sep="\t", quote=F, col.names=F, row.names=F); write.table(dup, "dupLengths.txt",sep="\t", quote=F, col.names=F, row.names=F)'   ## genBankAcc refSeqAcc(i.e.ncbiName) length ucscScaf ucscName start end
$script_path/unAssembledMapToNCBIpsl NCBItoUCSCunAssembled_map.txt > psl/unAssembled_map.psl

cat psl/unAssembled_map.psl | awk -F "\t" '{ print $10 }' > uniqLengths_ids.txt
grep -F -f uniqLengths_ids.txt eca_ref_EquCab2.0_unplaced.fa | sed 's/^>//g' > uniqLengths_Fullids.txt
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp eca_ref_EquCab2.0_unplaced.fa --output_fasta_fp eca_ref_EquCab2.0_unplaced_remain.fa --seq_id_fp uniqLengths_Fullids.txt --negate
#grep -F -w -f dupLengths.txt NCBIunAssembled_map.txt | awk -F "\t" '{ print $7 }' > dupLengths_ids.txt
#grep -F -f dupLengths_ids.txt eca_ref_EquCab2.0_unplaced.fa | sed 's/^>//g' > dupLengths_Fullids.txt
#module load QIIME/1.8.0
#filter_fasta.py --input_fasta_fp eca_ref_EquCab2.0_unplaced.fa --output_fasta_fp eca_ref_EquCab2.0_unplaced_remain.fa --seq_id_fp dupLengths_Fullids.txt

## make .ooc file
# http://redmine.soe.ucsc.edu/forum/index.php?t=msg&goto=1060&S=e9c689948dbd0ca84375e45882999599
module load ucscUtils/262
faSize $genome  ## 2484532062 bases (55741889 N's 2428790173 real 1433784199 upper 995005974 lower) in 34 sequences in 1 files
awk 'BEGIN{printf "%.6f\n", 2428790173 / 2897310462 * 1024}'  ## 858.410298
module load BLAT/36
blat $genome /dev/null /dev/null -tileSize=11 -makeOoc=equCab2.11.ooc -repMatch=850
## run blat
perl $script_path/splitFasta.pl eca_ref_EquCab2.0_unplaced.fa 500
for i in subset*_eca_ref_EquCab2.0_unplaced.fa eca_ref_EquCab2.0_chr*.fa; do
qsub -v genome=$genome,quary=$i,ooc=$"equCab2.11.ooc" $script_path/myblat.sh; done
# select hits where quary start (col#12)=0 and the quary end (col#13)=quary size (col# 11)
cd psl
mkdir finished
for f in subset*_eca_ref_EquCab2.0_unplaced.psl;do if [ $(wc -l $f|awk '{ print $1 }') != "0" ]; then
echo $f
cat $f | awk '($12 == 0 && $11 == $13)' >> finished/successful.psl
mv $f finished/.
fi; done
rm subset*_eca_ref_EquCab2.0_unplaced.psl
cd finished
grep -F -w -v -f ../../uniqLengths_ids.txt successful.psl > successful_forDup.psl
# make sure all target co-ordinates are matching expected co-ordinates in UCSCunAssembled_map.txt
cat successful_forDup.psl | awk -F "\t" '{ print $16 }' > successful_forDup_tStart.txt
grep -F -w -f successful_forDup_tStart.txt ../../UCSCunAssembled_map.txt > UCSCunAssembled_successDup_map.txt
wc -l successful_forDup.psl
wc -l UCSCunAssembled_successDup_map.txt
# correct for effect of repeats causing reciprocal inserations and deletions
$script_path/pslTopsl successful_forDup.psl > successful_forDup_cor.psl

# run blat
cd ../../
perl $script_path/splitFasta.pl eca_ref_EquCab2.0_unplaced_remain.fa 700
for i in subset*_eca_ref_EquCab2.0_unplaced_remain.fa; do
qsub -v genome="$genome_dir/chrUn.fa",quary=$i $script_path/myblat2.sh; done
# select hits where quary start (col#12)=0 and the quary end (col#13)=quary size (col# 11)
cd psl
mkdir finished2
for f in subset*_eca_ref_EquCab2.0_unplaced_remain.psl;do if [ $(wc -l $f|awk '{ print $1 }') != "0" ]; then
echo $f
cat $f | awk '($12 == 0 && $11 == $13)' >> finished2/successful.psl
mv $f finished2/.
fi; done
rm subset*_eca_ref_EquCab2.0_unplaced_remain.psl
cd finished2
grep -F -w -v -f ../../uniqLengths_ids.txt successful.psl > successful_forDup.psl
# make sure all target co-ordinates are matching expected co-ordinates in UCSCunAssembled_map.txt
cat successful_forDup.psl | awk -F "\t" '{ print $16 }' > successful_forDup_tStart.txt
grep -F -w -f successful_forDup_tStart.txt ../../UCSCunAssembled_map.txt > UCSCunAssembled_successDup_map.txt
wc -l successful_forDup.psl
wc -l UCSCunAssembled_successDup_map.txt
# correct for effect of repeats causing reciprocal inserations and deletions
$script_path/pslTopsl successful_forDup.psl > successful_forDup_cor.psl

# merge the psl files of the 2 blat runs then exclude replicates
cd ../
cat finished/successful_forDup_cor.psl finished2/successful_forDup_cor.psl > successful_forDup_cor.psl
sort successful_forDup_cor.psl | uniq -u > successful_forDup_cor2.psl

# exclude the recognized scaffolds and try to identify more by size matching
cd ../
grep -F -w -f dupLengths.txt NCBIunAssembled_map.txt > dupNCBIunAssembled_map.txt
grep -F -w -f dupLengths.txt UCSCunAssembled_map.txt > dupUCSCunAssembled_map.txt
Rscript -e 'x=read.table("dupNCBIunAssembled_map.txt", sep="\t"); y=read.table("dupUCSCunAssembled_map.txt", sep="\t"); fd=read.table("psl/successful_forDup_cor2.psl", sep="\t"); x2=x[!(x$V7 %in%  fd$V10),]; y2=y[!(y$V4 %in%  fd$V16),]; dup=unique(x2[duplicated(x2$V9),9]); xq=x2[!(x2$V9 %in%  dup),]; yq=y2[!(y2$V2 %in% dup),]; m=merge(xq,yq,by.x="V9",by.y="V2",all.x=TRUE,all.y=TRUE); m2=m[,c(6,8,1,11,12,13,14)]; write.table(m2, "NCBItoUCSCunAssembled_map2.txt",sep="\t", quote=F, col.names=F, row.names=F); write.table(dup, "dupLengths2.txt",sep="\t", quote=F, col.names=F, row.names=F)'   ## genBankAcc refSeqAcc(i.e.ncbiName) length ucscScaf ucscName start end
$script_path/unAssembledMapToNCBIpsl NCBItoUCSCunAssembled_map2.txt > psl/unAssembled_map2.psl

cat psl/unAssembled_map*.psl psl/successful_forDup_cor2.psl | awk -F "\t" '{ print $10 }' > uniqLengths_ids2.txt
grep -F -f uniqLengths_ids2.txt eca_ref_EquCab2.0_unplaced.fa | sed 's/^>//g' > uniqLengths_Fullids2.txt
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp eca_ref_EquCab2.0_unplaced.fa --output_fasta_fp eca_ref_EquCab2.0_unplaced_remain2.fa --seq_id_fp uniqLengths_Fullids2.txt --negate

# run blat
perl $script_path/splitFasta.pl eca_ref_EquCab2.0_unplaced_remain2.fa 500
for i in subset*_eca_ref_EquCab2.0_unplaced_remain2.fa; do
#blat "$genome" "$quary" -q=dna -extendThroughN -tileSize=18 -minMatch=4 -maxGap=0 -minIdentity=100 -noHead psl/`basename $quary .fa`.psl;done
qsub -v genome="$genome_dir/chrUn.fa",quary=$i $script_path/myblat2.sh; done


######################################################################
#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################
mkdir psl
# assembled chromosomes are too large and require very large memory. We can edit the psl files. All required data can be found in the fasta headers of ncbi files
wget --no-directories ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002305.2.assembly.txt
grep -v "^#" GCF_000002305.2.assembly.txt | awk '{ sub("\r$", ""); print }' > GCF_000002305.2.assembly_noCommonts.txt
cat GCF_000002305.2.assembly_noCommonts.txt | awk -F "\t" '$10 != "na"' > chromosome_map.txt
$script_path/chromosomeMapToNCBIpsl chromosome_map.txt > psl/chromosome_map.psl


# predict mapping of unassemblied scaffolds by size
cat GCF_000002305.2.assembly_noCommonts.txt | awk -F "\t" '$10 == "na"' | sort -t$'\t' -k9,9nr > NCBIunAssembled_map.txt
wget --no-directories http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/ctgPos2.txt.gz
gunzip ctgPos2.txt.gz
cat ctgPos2.txt | awk '$3 == "chrUn"' > UCSCunAssembled_map.txt
cat ctgPos2.txt | awk '$3 == "chrUn"' | awk -F "\t" -v OFS='\t' '{ print "chrUn",$4,$4+1400,$1 }' > ctgPos2_bait.txt
module load BEDTools/2.24.0
bedtools getfasta -fi "$genome_dir/chrUn.fa" -bed ctgPos2_bait.txt -fo UCSC_baits.tab -name -tab
grep -A20 "^>" eca_ref_EquCab2.0_unplaced.fa | grep -v "^--" | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | awk 'BEGIN{RS=">"}{gsub("\n","\t",$0); print $0}' | sed '/^$/d' > ncbi_prey.tab
Rscript -e 'x=read.table("UCSC_baits.tab", sep="\t"); y=read.table("ncbi_prey.tab", sep="\t"); x$V3=tolower(x$V2); y$V3=tolower(y$V2); m=merge(x,y,by.x="V3",by.y="V3",all.x=TRUE,all.y=TRUE); m2=m[,c(2,4)]; z=data.frame(do.call('rbind', strsplit(as.character(m2$V1.y), '|', fixed=TRUE))); m2$ncbiName=z$X4; x2=read.table("NCBIunAssembled_map.txt", sep="\t"); y2=read.table("UCSCunAssembled_map.txt", sep="\t"); m3=merge(x2,m2,by.x="V7",by.y="ncbiName",all.x=TRUE,all.y=TRUE); m4=merge(m3,y2,by.x="V1.x",by.y="V1",all.x=TRUE,all.y=TRUE); m5=m4[,c(7,2,13,1,14,15,16)]; write.table(m5, "NCBItoUCSCunAssembled_bait_map.txt",sep="\t", quote=F, col.names=F, row.names=F);'   ## genBankAcc refSeqAcc(i.e.ncbiName) ncbilength ucscScaf ucscName start end
$script_path/unAssembledMapToNCBIpsl NCBItoUCSCunAssembled_bait_map.txt > psl/unAssembled_bait_map.psl

cd psl
cat chromosome_map.psl unAssembled_bait_map.psl > NCBItoUCSC_map.psl

cd ../

sed 's/gi|\(.*\)|ref|\(.*\)| \(.*$\)/\2/g' ncbi_genome.fa > ncbi_genome2.fa
$HOME/bin/UCSC_kent_commands/axtChain -linearGap=medium -psl psl/NCBItoUCSC_map.psl -faT $genome -faQ ncbi_genome2.fa NCBItoUCSC_map.chain
$HOME/bin/UCSC_kent_commands/chainSort NCBItoUCSC_map.chain NCBItoUCSC_map.sorted.chain
cat psl/NCBItoUCSC_map.psl | awk -F "\t" -v OFS='\t' '{ print $10,$11 }' > ncbi.chromInfo
cat psl/NCBItoUCSC_map.psl | awk -F "\t" -v OFS='\t' '{ print $14,$15 }' | uniq > ucsc.chromInfo
$HOME/bin/UCSC_kent_commands/chainNet NCBItoUCSC_map.sorted.chain ucsc.chromInfo ncbi.chromInfo NCBItoUCSC_map.net UCSCtoNCBI_map.net
$HOME/bin/UCSC_kent_commands/netChainSubset NCBItoUCSC_map.net NCBItoUCSC_map.sorted.chain NCBItoUCSC.liftOver







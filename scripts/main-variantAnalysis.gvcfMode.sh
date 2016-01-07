#!/bin/sh
myRoot=$"/mnt/ls15/scratch/users/mansourt/Tamer"
source $myRoot/config.txt

## http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq
## https://www.broadinstitute.org/gatk/guide/article?id=3893
## What is a VCF and how should I interpret it?
## http://gatkforums.broadinstitute.org/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it
###########################################################################################
## Both Tophat and MergeSamFiles from picard tools sort the output BAM file by coordinates by default
## to double check you can run this samttol command for the target BAM file
#module load SAMTools/0.1.19
#samtools view -H accepted_hits_RG.bam | grep "SO:"

##########################
## mark duplicates
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/run_markDuplicates.sh "$sample_list" "$script_path/markDuplicates.sh";
done < $horse_trans/working_list_Retina.txt

## Check for successful mark duplicates
## To be added
## Tip: the line before last in .e file has "MarkDuplicates done"
## for f in $prepData/*/*/tophat_output/tophat_*/MarkDuplicates.e*; do grep "MarkDuplicates done" $f | wc -l; done
##########################
## reorder the BAM file to match the order of the GATK reference dictionary
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/run_reorderBAM.sh "$gatk_ref" "$sample_list" "$script_path/reorderBAM.sh";
done < $horse_trans/working_list_Cerebellum.txt

## Check for successful BAM reordering
## To be added
## Tip: the line before last in .e file has "ReorderSam done"
## for f in $prepData/*/*/tophat_output/tophat_*/reorderBAM.e*; do grep "ReorderSam done" $f | wc -l; done
##########################
## BuildBamIndex
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/run_BuildBamIndex.sh "$sample_list" "$script_path/buildBamIndex.sh";
done < $horse_trans/working_list_Cerebellum.txt

# Check for successful BAM indexing
## To be added
## Tip: the line before last in .e file has "BuildBamIndex done"
## for f in $prepData/*/*/tophat_output/tophat_*/buildBamIndex.e*; do grep "BuildBamIndex done" $f | wc -l; done
##########################
## Split'N'Trim and reassign mapping qualities
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/run_SplitNCigarReads.sh "$gatk_ref" "$sample_list" "$script_path/splitNCigarReads.sh";
done < $horse_trans/working_list_Cerebellum.txt

# Check for successful split
## To be added
## Tip: the line before last in .e file has "BuildBamIndex done"
## for f in $prepData/*/*/tophat_output/tophat_*/splitNCigarReads.e*; do echo $f; grep "Total runtime" $f | wc -l; done > temp
##########################
####  Realign Indels
#### This step used to be very important when the the variant callers were position-based (such as UnifiedGenotyper) but now that we have assembly-based variant callers (such as HaplotypeCaller) it is less important. We still perform indel realignment because we think it may improve the accuracy of the base recalibration model in the next step, but this step may be made obsolete in the near future.

## RealignerTargetCreator
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/tophat_output
#  sample_list=$work_dir/tophat_output/sample_list.txt
#  target_bam=$"split.bam"
#  bash ${script_path}/run_RealignerTargetCreator.sh "$knownIndels" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/realignerTargetCreator.sh";
#done < $horse_trans/working_list_Retina.txt

## indelRealigner
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/tophat_output
#  sample_list=$work_dir/tophat_output/sample_list.txt
#  target_bam=$"split.bam"
#  bash ${script_path}/run_indelRealigner.sh "$knownIndels" "gatk.intervals" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/indelRealigner.sh";
#done < $horse_trans/working_list_Retina.txt
##########################
## Variant calling per sample
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/tophat_output
#  sample_list=$work_dir/tophat_output/sample_list.txt
#  target_bam=$"split.bam"
#  bash ${script_path}/run_haplotypeCaller.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller.sh";
#done < $horse_trans/working_list_Retina.txt

## variant calling by sample
while read work_dir; do
  echo $work_dir
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_haplotypeCaller_GVCF.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_GVCF.sh"
done < $horse_trans/working_list_Retina.txt

# Check for successful variant calling
while read work_dir; do
  cd $work_dir/tophat_output
  for f in tophat_*/haplotypeCaller_multi.e*;do if [ -f $f ];then echo $f; grep "done" $f | wc -l;fi; done > haplotypeCaller_gvcf.temp
  grep -B1 "^0" haplotypeCaller_gvcf.temp | grep -v "^--" | grep -v "^0" > haplotypeCaller_gvcf.temp.redo
  while read f;do ls -tral $f;done < haplotypeCaller_gvcf.temp.redo > haplotypeCaller_gvcf.temp.redo.size
  for f in tophat_*;do if [ ! -f $f/haplotypeCaller_multi.e* ];then echo $work_dir/tophat_output/$f;fi; done > haplotypeCaller_gvcf.sample_list2.txt
done < $horse_trans/working_list_NoPBMCs.txt

while read work_dir; do
  cd $work_dir/tophat_output
  while read f;do
    output=$(echo $f | awk -F "/" '{ print $1 }');
    name=${output#tophat_};
    echo $work_dir/tophat_output/$output
    rm $output/{split.g.vcf,haplotypeCaller_multi.e*,haplotypeCaller_multi.o*}
  done < haplotypeCaller_gvcf.temp.redo > haplotypeCaller_gvcf.sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_haplotypeCaller_GVCF.sh "$knownSNPs" "$gatk_ref" "haplotypeCaller_gvcf.sample_list.txt" "$target_bam" "$script_path/haplotypeCaller_GVCF.sh"
  bash ${script_path}/run_haplotypeCaller_GVCF.sh "$knownSNPs" "$gatk_ref" "haplotypeCaller_gvcf.sample_list2.txt" "$target_bam" "$script_path/haplotypeCaller_GVCF.sh"
done < $horse_trans/working_list_NoPBMCs.txt

## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  rm -f $work_dir/tophat_output/gvcf_list.txt
  for f in $work_dir/tophat_output/tophat_*/*.g.vcf; do if [ -f $f ]; then
    echo $f >> $work_dir/tophat_output/gvcf_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

> $horse_trans/all_g.vcfs.txt
while read work_dir; do
  cat $work_dir/tophat_output/gvcf_list.txt >> $horse_trans/all_g.vcfs.txt
done < $horse_trans/working_list_NoPBMCs.txt

## joint genotyping
mkdir $horse_trans/Var_gvcf
cd $horse_trans/Var_gvcf
sample_list=$horse_trans/all_g.vcfs.txt
bash ${script_path}/run_genotypeGVCF.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$script_path/genotypeGVCF.sh"
##########################
bash $script_path/variantFiltration.sh "$gatk_ref" "GenotypeGVCFs_output.vcf"
grep "^#" GenotypeGVCFs_output_filtered.vcf > GenotypeGVCFs_PASS.vcf
grep "PASS" GenotypeGVCFs_output_filtered.vcf >> GenotypeGVCFs_PASS.vcf
grep -v "PASS" GenotypeGVCFs_output_filtered.vcf > GenotypeGVCFs_FAILED.vcf

awk '/#/{print;next}{if($5 !~ /,/){print}}' GenotypeGVCFs_PASS.vcf > GenotypeGVCFs_monoAllel.vcf
#awk '/#/{print;next}{if($5 ~ /,/){print}}' GenotypeGVCFs_PASS.vcf > GenotypeGVCFs_multiAllel.vcf

awk '/#/{print;next}{if($5 !~ /,/ && (length($5)>1 || length($4)>1)){print}}' GenotypeGVCFs_PASS.vcf > GenotypeGVCFs_monoAllel_indels.vcf
#grep "^#" GenotypeGVCFs_PASS.vcf > GenotypeGVCFs_PASS_indels.vcf
#grep -v "^#" GenotypeGVCFs_PASS.vcf | awk 'length($4)>1 || length($5)>1' >> GenotypeGVCFs_PASS_indels.vcf

awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' GenotypeGVCFs_PASS.vcf > GenotypeGVCFs_monoAllel_SNPs.vcf
#grep "^#" HC_output_ploidy1_haplo1_PASS.vcf > HC_output_ploidy1_haplo1_PASS_snps.vcf
#grep -v "^#" HC_output_ploidy1_haplo1_PASS.vcf | awk 'length($4)==1 && length($5)==1' >> HC_output_ploidy1_haplo1_PASS_snps.vcf

##########################
## explore the frequency of the ALT alleles
#INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
#INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
#INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
## AF=AC/AN
echo -e 'AC\tAF\tAN' > varFreq.txt
grep -v "^#" GenotypeGVCFs_monoAllel.vcf | awk '{ print $8 }' | awk -F ';' -v OFS='\t' '{ print $1,$2,$3 }' | sed 's/..=//g' >> varFreq.txt
tail -n+2 varFreq.txt | awk -F '\t' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF.freq
tail -n+2 varFreq.txt | awk -F '\t' '($1>2){A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF2.freq
tail -n+2 varFreq.txt | awk -F '\t' '($1>3){A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF3.freq
tail -n+2 varFreq.txt | awk -F '\t' '($1>4){A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF4.freq

echo -e 'AC\tAF\tAN' > varFreq_indels.txt
grep -v "^#" GenotypeGVCFs_monoAllel_indels.vcf | awk '{ print $8 }' | awk -F ';' -v OFS='\t' '{ print $1,$2,$3 }' | sed 's/..=//g' >> varFreq_indels.txt
tail -n+2 varFreq_indels.txt | awk -F '\t' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF_indels.freq
tail -n+2 varFreq_indels.txt | awk -F '\t' '($1>2){A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF_indels2.freq
tail -n+2 varFreq_indels.txt | awk -F '\t' '($1>3){A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF_indels3.freq
tail -n+2 varFreq_indels.txt | awk -F '\t' '($1>4){A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF_indels4.freq

##########################
## isolate common variants (AF >= 0.9)
grep "^#" GenotypeGVCFs_monoAllel.vcf > GenotypeGVCFs_monoAllel_sig.vcf
grep -v "^#" GenotypeGVCFs_monoAllel.vcf | awk -F '[\t=;]' '($9 >3 && $11 >= 0.8)'  >> GenotypeGVCFs_monoAllel_sig.vcf

grep -v "^#" GenotypeGVCFs_monoAllel_sig.vcf | wc -l  ##227717

grep "^#" GenotypeGVCFs_monoAllel_indels.vcf > GenotypeGVCFs_monoAllel_indels_sig.vcf
grep -v "^#" GenotypeGVCFs_monoAllel_indels.vcf | awk -F '[\t=;]' '($9 >3 && $11 >= 0.8)'  >> GenotypeGVCFs_monoAllel_indels_sig.vcf

grep -v "^#" GenotypeGVCFs_monoAllel_indels_sig.vcf | wc -l  ##38433

###########################
#### liftover the genome variance file to the transcriptome
assembly="$tissue_Cuffmerge/$cuffmerge_output/filtered"
mkdir $assembly/varFixed
cp $assembly/{merged.gtf,merged.gpred} $assembly/varFixed/.
cd $assembly/varFixed

## change merged.gtf to be all +ve starnded
cat merged.gtf | awk -F "\t" -v OFS='\t' '{ $7 = "+"; print }' > merged_allPositive.gtf

## create the chain file
$script_path/UCSC_kent_commands/gtfToGenePred merged_allPositive.gtf merged_allPositive.gpred
$script_path/UCSC_kent_commands/genePredToFakePsl -chromSize=$genome_dir/$UCSCgenome.chrom.sizes file merged_allPositive.gpred merged_allPositive.psl /dev/null
$script_path/UCSC_kent_commands/pslToChain merged_allPositive.psl genomeToTrans_map.chain ## should I use "pslToChain" or "axtChain" ??
$script_path/UCSC_kent_commands/chainSort genomeToTrans_map.chain genomeToTrans_map.sorted.chain
chain=$assembly/varFixed/genomeToTrans_map.sorted.chain

## Construct the transcript fasta file by decoderUtil
bash $script_path/run_genome_to_cdna_fasta.sh "merged_allPositive.gtf" "$genome" "transcripts_allPositive.fa" "$script_path/genome_to_cdna_fasta.sh"

## lift over the variants from the genome to the transcriptome
bash ${script_path}/run_gatk-index.sh transcripts_allPositive.fa
dict=$assembly/varFixed/transcripts_allPositive.dict
inputVCF="GenotypeGVCFs_monoAllel_sig.vcf"
outputVCF=$assembly/varFixed/${inputVCF%.vcf}_trans.vcf
bash $script_path/liftoverVariants.sh "$gatk_ref" "$horse_trans/Var_gvcf/$inputVCF" "$chain" "$dict" "$outputVCF" ## Converted 32249 records; failed to convert 195468 records.
index=$assembly/varFixed/transcripts_allPositive.fa.fai
outputVCF_refSorted=$assembly/varFixed/${inputVCF%.vcf}_transRefSorted.vcf
grep -v "^#" $outputVCF > temp
perl $script_path/sortByRef.pl temp $index > temp2
grep "^#" $outputVCF > $outputVCF_refSorted
cat temp2 >> $outputVCF_refSorted
rm temp*
finalVCF=$assembly/varFixed/${inputVCF%.vcf}_transFinal.vcf
bash $script_path/filterLiftedVariants.sh "transcripts_allPositive.fa" "$outputVCF_refSorted" "$finalVCF"  ## Filtered 1 records out of 32249 total records

## correct the transcripts with VCF
bash $script_path/fastaAlternateReferenceMaker.sh "transcripts_allPositive.fa" "$finalVCF" "corTranscripts_allPositive.fa"

## transfer the fasta header names to the new file
cat corTranscripts_allPositive.fa | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' > corTranscripts_allPositive_unwrap.fa
n=1
while true; do
  read -r f1 <&3 || break
  read -r f2 <&4 || break
  if [ $(($n % 2)) -eq 0 ];then echo "$f2"; else echo "$f1"; fi;
  let "n++"
done 3<transcripts_allPositive.fa 4<corTranscripts_allPositive_unwrap.fa > corTranscripts_allPositive_fixed.fa
#grep "^>" transcripts_allPositive.fa > originalHeaders
#grep "^>" corTranscripts_allPositive.fa > newHeaders
#paste originalHeaders newHeaders > header_map
#while IFS="$(printf '\t')" read OID NID;do sed -i "s/$NID$/$OID/" corTranscripts_allPositive.fa; done < header_map

## assess the rate of edits per transcript
grep "^>" transcripts_allPositive.fa | sed 's/>//' | awk '{ print $1 }' > ids
while read id;do grep $id <(grep -v "^#" $finalVCF) | wc -l; done < ids > counts
cat counts | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > freq

## determine the list of negativeStrandTransHeaders
cat merged.gpred | awk '{ if ( $3 == "-" ) { print $1 };}' > negativeStrandTransHeaders
grep -F -f negativeStrandTransHeaders corTranscripts_allPositive_fixed.fa | sed 's/>//' > negativeStrandTransCompleteHeaders

## split the corrected fasta into negativeStrandTrans & nonNegativeStrandTrans
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp corTranscripts_allPositive_fixed.fa --output_fasta_fp corNegativeStrandTrans.fa --seq_id_fp negativeStrandTransCompleteHeaders
filter_fasta.py --input_fasta_fp corTranscripts_allPositive_fixed.fa --output_fasta_fp corNonNegativeStrandTrans.fa --seq_id_fp negativeStrandTransCompleteHeaders --negate

## reverese complement the negativeStrandTrans then unwrap
module load EMBOSS/6.5.7
revseq fasta::corNegativeStrandTrans.fa corNegativeStrandTrans_rev.fa
cat corNegativeStrandTrans_rev.fa | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' > corNegativeStrandTrans_rev_unwrap.fa

## merge the files to restore one corrected fasta file
cat corNonNegativeStrandTrans.fa corNegativeStrandTrans_rev_unwrap.fa > transcripts.fasta
## may be we need to sort according to gff3 produced by transdecoder

cd $genome_dir
$script_path/UCSC_kent_commands/faToTwoBit $genome genome.2bit
## http://crc.ibest.uidaho.edu/help/Applications/BLAT.html
blat -makeOoc=11.ooc $genome_dir/genome.2bit /dev/null /dev/null

mkdir $assembly/varFixed/splitTrans
cd $assembly/varFixed/splitTrans
cp ../transcripts.fasta .
$script_path/splitFasta.pl transcripts.fasta 30
for f in subset*_transcripts.fasta;do
label=${f%_transcripts.fasta}
qsub -v genome="$genome_dir/genome.2bit",transcriptome=$f,oocFile="$genome_dir/11.ooc",output=$label."varFixed.psl" $script_path/blat.sh;
done
cat subset*.varFixed.psl >> ../varFixed.psl
cd ../
sort -k10,10 -k1,1rg varFixed.psl | sort -u -k10,10 --merge > varFixed_best.psl
$script_path/UCSC_kent_commands/pslToBed varFixed_best.psl varFixed_best.bed
$script_path/UCSC_kent_commands/bedToGenePred varFixed_best.bed varFixed_best.GenePred
$script_path/UCSC_kent_commands/genePredToGtf file varFixed_best.GenePred varFixed_best.gtf








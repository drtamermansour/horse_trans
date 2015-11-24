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
done < $horse_trans/working_list_Cerebellum.txt

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

##########################
## Variant calling per library
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/tophat_output
#  sample_list=$work_dir/tophat_output/sample_list.txt
#  target_bam=$"split.bam"
#  bash ${script_path}/run_haplotypeCaller_multi.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_multi.sh";
#done < $horse_trans/working_list_Retina.txt

##########################
## combined Variant calling
> $horse_trans/all_samples.txt
while read work_dir; do
  echo $work_dir
  cat $work_dir/tophat_output/sample_list.txt >> $horse_trans/all_samples.txt ## for the 1st run, I used only 49 samples (6 cerebellum samples were not included). This list of samples are saved as $horse_trans/all_samples_temp.txt
done < $horse_trans/working_list_NoPBMCs.txt
mkdir $horse_trans/Var_merge
cd $horse_trans/Var_merge
sample_list=$horse_trans/all_samples.txt
target_bam=$"split.bam"
bash ${script_path}/run_haplotypeCaller_multi.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_multi.sh"   ## output is HC_output_ploidy1_haplo1.vcf
bash $script_path/variantFiltration.sh "$gatk_ref" "HC_output_ploidy1_haplo1.vcf"
grep "^#" HC_output_ploidy1_haplo1_filtered.vcf > HC_output_ploidy1_haplo1_PASS.vcf
grep "PASS" HC_output_ploidy1_haplo1_filtered.vcf >> HC_output_ploidy1_haplo1_PASS.vcf
grep -v "PASS" HC_output_ploidy1_haplo1_filtered.vcf > HC_output_ploidy1_haplo1_FAILED.vcf

grep "^#" HC_output_ploidy1_haplo1_PASS.vcf > HC_output_ploidy1_haplo1_PASS_indels.vcf
grep -v "^#" HC_output_ploidy1_haplo1_PASS.vcf | awk 'length($4)>1 || length($5)>1' >> HC_output_ploidy1_haplo1_PASS_indels.vcf
grep "^#" HC_output_ploidy1_haplo1_PASS.vcf > HC_output_ploidy1_haplo1_PASS_snps.vcf
grep -v "^#" HC_output_ploidy1_haplo1_PASS.vcf | awk 'length($4)==1 && length($5)==1' >> HC_output_ploidy1_haplo1_PASS_snps.vcf


##########################
## explore the frequency of the ALT alleles
cd $horse_trans/Var_merge
echo -e 'AC\tAF\tAN' > varFreq.txt
grep -v "^#" HC_output_ploidy1_haplo1_PASS.vcf | awk '{ print $8 }' | awk -F ';' -v OFS='\t' '{ print $1,$2,$3 }' | sed 's/..=//g' >> varFreq.txt
tail -n+2 varFreq.txt | awk -F '\t' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF.freq

echo -e 'AC\tAF\tAN' > varFreq_indels.txt
grep -v "^#" HC_output_ploidy1_haplo1_PASS_indels.vcf | awk '{ print $8 }' | awk -F ';' -v OFS='\t' '{ print $1,$2,$3 }' | sed 's/..=//g' >> varFreq_indels.txt
tail -n+2 varFreq_indels.txt | awk -F '\t' '{A[$2]++}END{for(i in A)print i,A[i]}' | sort -k1,1nr > AF_indels.freq
##########################
## isolate common indels (AF >= 0.9)
cd $horse_trans/Var_merge
grep "^#" HC_output_ploidy1_haplo1_PASS_indels.vcf > HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf
grep -v "^#" HC_output_ploidy1_haplo1_PASS_indels.vcf | awk -F '[\t=;]' '$11 >= 0.9'  >> HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf

grep -v "^#" HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf | wc -l  ##4725

## isolate common SNPs (AF >= 0.9)
grep "^#" HC_output_ploidy1_haplo1_PASS_snps.vcf > HC_output_ploidy1_haplo1_PASS_snps_AF0.9.vcf
grep -v "^#" HC_output_ploidy1_haplo1_PASS_snps.vcf | awk -F '[\t=;]' '$11 >= 0.9'  >> HC_output_ploidy1_haplo1_PASS_snps_AF0.9.vcf

grep -v "^#" HC_output_ploidy1_haplo1_PASS_snps_AF0.9.vcf | wc -l  ##4725

##########################
## RealignerTargetCreator
cd $horse_trans/Var_merge
bash ${script_path}/run_RealignerTargetCreator_forKnowns.sh "HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf" "$gatk_ref" "$script_path/realignerTargetCreator_forKnowns.sh";

## indelRealigner
indels=$horse_trans/Var_merge/HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf
intervals=$horse_trans/Var_merge/gatk.intervals
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_indelRealigner.sh "$indels" "$intervals" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/indelRealigner_forKnowns.sh";
done < $horse_trans/working_list_NoPBMCs.txt

# Check for successful BAM indexing
## To be added
## for f in $prepData/*/*/tophat_output/tophat_*/indelRealigner.e*; do grep "Total runtime" $f | wc -l; done

##########################
#### Base Recalibration
#### We do recommend running base recalibration (BQSR). Even though the effect is also marginal when applied to good quality data, it can absolutely save your butt in cases where the qualities have systematic error modes.

## create list of known variants
known_var=$horse_trans/Var_merge/known_var.txt
> $known_var
echo "$knownSNPs" >> $known_var
echo "$horse_trans/Var_merge/HC_output_ploidy1_haplo1_PASS_snps_AF0.9.vcf" >> $known_var
echo "$horse_trans/Var_merge/HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf" >> $known_var

## baseRecalibrator- 1st round
## it might be better to pass all the library samples into one recalibaryion job
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.realigned.bam"
  bash ${script_path}/run_baseRecalibrator.sh "$known_var" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/baseRecalibrator_1st.sh"
done < $horse_trans/working_list_NoPBMCs.txt

# Check for successful baseRecalibrator
## To be added
## for f in $prepData/*/*/tophat_output/baseRecalibrator-1stR.e*; do echo $f; grep "Total runtime" $f | wc -l; done
## for f in $prepData/*/*/tophat_output/baseRecalibrator-1stR.e*; do echo $f; grep "reads were filtered out during the traversal" $f; done

## baseRecalibrator- 2nd round
## it might be better to pass all the library samples into one recalibaryion job
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.realigned.bam"
  bash ${script_path}/run_baseRecalibrator.sh "$known_var" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/baseRecalibrator_2nd.sh";
done < $horse_trans/working_list_NoPBMCs.txt

## Generate before/after plots
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  bash ${script_path}/run_createRecalPlots.sh "$gatk_ref" "$script_path/createRecalPlots.sh";
done < $horse_trans/working_list_NoPBMCs.txt

## Apply the recalibration to your sequence data
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.realigned.bam"
  bash ${script_path}/run_applyRecalib.sh "$gatk_ref" "$sample_list" "$target_bam" "$script_path/applyRecalib.sh"
done < $horse_trans/working_list_NoPBMCs.txt
###########################
## variant calling by library
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"recal_reads.bam"
  bash ${script_path}/run_haplotypeCaller_multi.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_multi.sh"   ## output is HC_output_ploidy1_haplo1.vcf
done < $horse_trans/working_list_NoPBMCs.txt

# Check for successful baseRecalibrator
## To be added
## for f in $prepData/*/*/tophat_output/haplotypeCaller_multi.e*; do echo $f; grep "Total runtime" $f | wc -l; done
## for f in $prepData/*/*/tophat_output/haplotypeCaller_multi.e*; do echo $f; grep "reads were filtered out during the traversal" $f; done

## Combine Variants
> $horse_trans/all_variants.txt
while read work_dir; do
  echo $work_dir/tophat_output/HC_output_ploidy1_haplo1.vcf >> $horse_trans/all_variants.txt
done < $horse_trans/working_list_NoPBMCs.txt
mkdir -p $horse_trans/Var_merge2
cd $horse_trans/Var_merge2
sample_list=$horse_trans/all_variants.txt
bash ${script_path}/run_combineVariants.sh "$gatk_ref" "$sample_list" "$script_path/combineVariants.sh"   ## output is merge_HC_output_ploidy1_haplo1.vcf

## fix annoatation
target_bam="recal_reads.bam"
sample_list="$horse_trans/all_samples.txt"
bash ${script_path}/run_variantAnnotator.sh "$gatk_ref" "$sample_list" "$target_bam" "$script_path/variantAnnotator.sh"   ## output is HC_output_ploidy1_haplo1_Ann.vcf
###########################
## combined Variant calling
#mkdir $horse_trans/Var_merge3
#cd $horse_trans/Var_merge3
#sample_list=$horse_trans/all_samples.txt
#target_bam=$"recal_reads.bam"
#bash ${script_path}/run_haplotypeCaller_multi.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_multi.sh"   ## output is HC_output_ploidy1_haplo1.vcf
###########################

bash $script_path/variantFiltration.sh "$gatk_ref" "HC_output_ploidy1_haplo1.vcf"
grep "^#" HC_output_ploidy1_haplo1_filtered.vcf > HC_output_ploidy1_haplo1_PASS.vcf
grep "PASS" HC_output_ploidy1_haplo1_filtered.vcf >> HC_output_ploidy1_haplo1_PASS.vcf
grep -v "PASS" HC_output_ploidy1_haplo1_filtered.vcf > HC_output_ploidy1_haplo1_FAILED.vcf

grep "^#" HC_output_ploidy1_haplo1_PASS.vcf > HC_output_ploidy1_haplo1_PASS_indels.vcf
grep -v "^#" HC_output_ploidy1_haplo1_PASS.vcf | awk 'length($4)>1 || length($5)>1' >> HC_output_ploidy1_haplo1_PASS_indels.vcf
grep "^#" HC_output_ploidy1_haplo1_PASS.vcf > HC_output_ploidy1_haplo1_PASS_snps.vcf
grep -v "^#" HC_output_ploidy1_haplo1_PASS.vcf | awk 'length($4)==1 && length($5)==1' >> HC_output_ploidy1_haplo1_PASS_snps.vcf


###########################
#### liftover the genome variance file to the transcriptome
assembly="$tissue_Cuffmerge/all_tissues/nonGuided_Cufflinks/nonGuided_Cuffmerge"
assembly="$tissue_Cuffmerge/Skin/refGeneGuided_Cufflinks/nonGuided_Cuffmerge"
cd $assembly

## change merged.gtf to be all +ve starnded
cat merged.gtf | awk -F "\t" -v OFS='\t' '{ $7 = "+"; print }' > merged_allPositive.gtf

## create the chain file
$script_path/UCSC_kent_commands/gtfToGenePred merged_allPositive.gtf merged_allPositive.gpred
$script_path/UCSC_kent_commands/genePredToFakePsl -chromSize=$genome_dir/$UCSCgenome.chrom.sizes file merged_allPositive.gpred merged_allPositive.psl /dev/null
$script_path/UCSC_kent_commands/pslToChain merged_allPositive.psl genomeToTrans_map.chain ## should I use "pslToChain" or "axtChain" ??
$script_path/UCSC_kent_commands/chainSort genomeToTrans_map.chain genomeToTrans_map.sorted.chain
chain=$assembly/genomeToTrans_map.sorted.chain

## Construct the transcript fasta file by decoderUtil
bash $script_path/run_genome_to_cdna_fasta.sh "merged_allPositive.gtf" "$genome" "transcripts_allPositive.fa" "$script_path/genome_to_cdna_fasta.sh"

## correct the transcripts with VCF
bash ${script_path}/run_gatk-index.sh transcripts_allPositive.fa
dict=$assembly/transcripts_allPositive.dict
inputVCF="HC_output_ploidy1_haplo1_PASS_indels_AF0.9.vcf"
outputVCF=$assembly/${inputVCF%.vcf}_trans.vcf
bash $script_path/liftoverVariants.sh "$gatk_ref" "$horse_trans/Var_merge/$inputVCF" "$chain" "$dict" "$outputVCF"
index=$assembly/transcripts_allPositive.fa.fai
outputVCF_refSorted=$assembly/${inputVCF%.vcf}_transRefSorted.vcf
grep -v "^#" $outputVCF > temp
perl $script_path/sortByRef.pl temp $index > temp2
grep "^#" $outputVCF > $outputVCF_refSorted
cat temp2 >> $outputVCF_refSorted
rm temp*
finalVCF=$assembly/${inputVCF%.vcf}_transFinal.vcf
bash $script_path/filterLiftedVariants.sh "transcripts_allPositive.fa" "$outputVCF_refSorted" "$finalVCF"
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
seqret fastq-illumina::s1_pe.fq fastq::${out_path}/quake_data/s1_pe33.fq
revseq fasta::corNegativeStrandTrans.fa corNegativeStrandTrans_rev.fa
cat corNegativeStrandTrans_rev.fa | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' > corNegativeStrandTrans_rev_unwrap.fa

## merge the files to restore one corrected fasta file
mkdir $assembly/corTransdecoder
cat corNonNegativeStrandTrans.fa corNegativeStrandTrans_rev_unwrap.fa > $assembly/corTransdecoder/transcripts.fasta
## may be we need to sort according to gff3 produced by transdecoder









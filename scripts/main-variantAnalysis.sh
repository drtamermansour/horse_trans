#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
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
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_RealignerTargetCreator.sh "$knownIndels" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/realignerTargetCreator.sh";
done < $horse_trans/working_list_Retina.txt


## indelRealigner
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_indelRealigner.sh "$knownIndels" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/indelRealigner.sh";
done < $horse_trans/working_list_Retina.txt

##########################
#### Base Recalibration
#### We do recommend running base recalibration (BQSR). Even though the effect is also marginal when applied to good quality data, it can absolutely save your butt in cases where the qualities have systematic error modes.


## baseRecalibrator- 1st round
## it might be better to pass all the library samples into one recalibaryion job
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.realigned.bam"
  bash ${script_path}/run_baseRecalibrator.sh "$knownSNPs" "$knownIndels" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/baseRecalibrator_1st.sh";
done < $horse_trans/working_list_Retina.txt

## baseRecalibrator- 2nd round
## it might be better to pass all the library samples into one recalibaryion job
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.realigned.bam"
  bash ${script_path}/run_baseRecalibrator.sh "$knownSNPs" "$knownIndels" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/baseRecalibrator_2nd.sh";
done < $horse_trans/working_list_Retina.txt


## Generate before/after plots

## Apply the recalibration to your sequence data


##########################
## Variant calling per sample
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_haplotypeCaller.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller.sh";
done < $horse_trans/working_list_Retina.txt

##########################
## Variant calling per library
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  target_bam=$"split.bam"
  bash ${script_path}/run_haplotypeCaller_multi.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_multi.sh";
done < $horse_trans/working_list_Retina.txt

##########################
## combined Variant calling
> $horse_trans/all_samples.txt
while read work_dir; do
  echo $work_dir
  cat $work_dir/tophat_output/sample_list.txt >> $horse_trans/all_samples.txt
done < $horse_trans/working_list_NoPBMCs.txt
mkdir $horse_trans/Var_merge
cd $horse_trans/Var_merge
sample_list=$horse_trans/all_samples.txt
target_bam=$"split.bam"
bash ${script_path}/run_haplotypeCaller_multi.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_multi.sh";













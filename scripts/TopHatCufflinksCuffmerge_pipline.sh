#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################
## pipeline_OneSampleAtaTime_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/trimmed_RNA_reads ]; then
  rm -f $work_dir/trimmed_RNA_reads/sample_list.txt
  for f in $work_dir/trimmed_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
    echo $f >> $work_dir/trimmed_RNA_reads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## run Tophat on each sample
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/tophat_output
  cd $work_dir/tophat_output
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
  bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful tophat runs and trouble shooting the failed tophat jobs (require tophat-[SP]E.e & .o)
while read work_dir; do
  cd $work_dir/tophat_output
  sample_list=$work_dir/trimmed_RNA_reads/tophat_failedSamples.txt           ## define the path of empty file
  bash $script_path/check_tophat.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')
    echo "Failed tophat jobs in: "$work_dir
    bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
  fi
done < $horse_trans/working_list_NoPBMCs.txt
##################
## create summary for tophat run
while read work_dir; do
  for f in $work_dir/tophat_output/tophat_*; do
    echo ${f} >> $work_dir/tophat_output/allsample_summary.txt
    cd ${f}
    grep "overall read mapping rate" align_summary.txt >> ../allsample_summary.txt
    grep "concordant pair alignment rate" align_summary.txt >> ../allsample_summary.txt
  done
done < $horse_trans/working_list_NoPBMCs.txt
##################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  rm -f $work_dir/tophat_output/sample_list.txt
  for f in $work_dir/tophat_output/tophat_*; do if [ -d $f ]; then
  echo $f >> $work_dir/tophat_output/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/run_cufflinks.sh "$sample_list" "$refGTF_file" "$script_path/cufflinks.sh";
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful Cufflinks runs and trouble shooting the failed Cufflinks jobs (requires cufflinks.e)
while read work_dir; do
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/failedSamples.txt           ## define the path of empty file
  bash $script_path/check_cufflinks.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    echo "Failed Cufflinks jobs in: "$work_dir
    cat $sample_list
    bash ${script_path}/run_cufflinks.sh "$sample_list" "$refGTF_file" "$script_path/cufflinks.sh";
fi; done < $horse_trans/working_list_NoPBMCs.txt
############
## Assess computational utilization of cufflinks
cufflinks_utlization=$prepData/cufflinks_utlization.txt
> $cufflinks_utlization
while read work_dir; do
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/assess_cufflinks_utlization.sh "$sample_list" "$cufflinks_utlization"
done < $horse_trans/working_list_NoPBMCs.txt
##################
### Run cuffmerge: merge the sample assemblies and output merged.gtf in tophat_output/cuffmerge_output
output_noGTF=$"cuffmerge_output/withoutGTF"
#output_withRefGene=$"cuffmerge_output/withRefGene"
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cd $work_dir/tophat_output
  for dir in $work_dir/tophat_output/tophat_*; do if [ -d "$dir" ]; then
    echo "$dir"/transcripts.gtf; fi; done > assemblies.txt
  rm -fR $output_noGTF
  #rm -fR $output_withRefGene
  bash ${script_path}/cuffmerge_noGTF.sh "$genome" "$output_noGTF" $"assemblies.txt"
  #bash ${script_path}/cuffmerge_withRefGene.sh "$genome" "$output_withRefGene" $"assemblies.txt" "$refGTF_file"
fi; done < $horse_trans/working_list_NoPBMCs_NoCereb.txt
## http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#transfrag-class-codes
##################
### cuffmerge the lab-specific tissue assemblies into tissue specific assemblies
# prepare list of target tissues
rm -f $horse_trans/multi_lib_tissues.txt
for tissue_dir in $prepData/*; do if [[ -d $tissue_dir && $(ls $tissue_dir | wc -l) -gt 1 ]]; then
  echo $tissue_dir >> $horse_trans/multi_lib_tissues.txt; fi; done;

# For every target tissue, cuffmerge the lab-spefific assemblies into one assembly in $tissue_merge/cuffmerge/TISSUE.NAME/withORwithoutRefGuidence.gtf
tissue_Cuffmerge=$tissue_merge/cuffmerge
mkdir -p $tissue_Cuffmerge
cd $tissue_Cuffmerge
while read tissue_dir; do
  tissue=$(basename $tissue_dir)
  mkdir -p $tissue
  rm -f $tissue_dir/All_libraries_assemblies.txt
  for dir in $tissue_dir/*; do if [ -d "$dir" ]; then
    cat "$dir"/tophat_output/assemblies.txt >> $tissue_dir/All_libraries_assemblies.txt; fi; done
  output_noGTF=$tissue/withoutGTF
  rm -fR $output_noGTF
  bash ${script_path}/cuffmerge_noGTF.sh "$genome" "$output_noGTF" "$tissue_dir/All_libraries_assemblies.txt"
  #output_withRefGene=$tissue/withRefGene
  #rm -fR $output_withRefGene
  #bash ${script_path}/cuffmerge_withRefGene.sh "$genome" "$output_withRefGene" "$tissue_dir/All_libraries_assemblies.txt" "$refGTF_file"
done < $horse_trans/multi_lib_tissues.txt
##################
### cuffmerge all assemblies into one total assembly
# save the assembly in $horse_trans/total_merge/cuffmerge/withORwithoutRefGuidence.gtf
cd $tissue_Cuffmerge
rm -f $prepData/All_tissues_assemblies.txt
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cat "$work_dir"/tophat_output/assemblies.txt >> $prepData/All_tissues_assemblies.txt
done < $horse_trans/working_list.txt

mkdir -p all_tissues
output_noGTF=$"all_tissues/withoutGTF"
bash ${script_path}/cuffmerge_noGTF.sh "$genome" "$output_noGTF" "$prepData/All_tissues_assemblies.txt"
#output_withRefGene=$"all_tissues/withRefGene"
#bash ${script_path}/cuffmerge_withRefGene.sh "$genome" "$output_withRefGene" "$prepData/All_tissues_assemblies.txt" "$refGTF_file"
###################
## create list of assemblies from each library
## This is where you can edit the list to restrict the processing for certain target(s)
rm -f $prepData/merged_assemblies.txt
while read work_dir; do
  for dir in $work_dir/tophat_output/cuffmerge_output/withoutGTF; do if [ -d $dir ]; then
    echo ${dir#$prepData/} >> $prepData/merged_assemblies.txt;
  fi; done;
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## create list of assemblies for tissues of multiple libraries
rm -f $tissue_Cuffmerge/tissue_assemblies.txt
for tissue in $tissue_Cuffmerge/*; do
  echo ${tissue#$tissue_Cuffmerge/}/withoutGTF >> $tissue_Cuffmerge/tissue_assemblies.txt;
done
####################
## convert the gtf files into BigBed files & copy the BigBed files to the track hub directory
rm -f $horse_trans/merged_and_tissue_assemblies.txt
while read assembly; do
  echo $assembly
  cd $prepData/$assembly
  targetAss=$"merged.gtf"
  bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  if [ -f $"merged.BigBed" ];then
    identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    cp merged.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    echo $prepData/$assembly >> $horse_trans/merged_and_tissue_assemblies.txt;
fi; done < $prepData/merged_assemblies.txt

while read assembly; do
  echo $assembly
  cd $tissue_Cuffmerge/$assembly
  targetAss=$"merged.gtf"
  bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  if [ -f $"merged.BigBed" ];then
    identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    cp merged.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    echo $tissue_Cuffmerge/$assembly >> $horse_trans/merged_and_tissue_assemblies.txt;
fi; done < $tissue_Cuffmerge/tissue_assemblies.txt

## run icommand to push the file to iplant
## https://pods.iplantcollaborative.org/wiki/display/DS/Using+iCommands
## http://bioinformatics.plantbiology.msu.edu/display/IP/Moving+Data+from+HPCC+to+iPlant
#icd /iplant/home/drtamermansour/horseTrans
#while read assembly; do
#  echo $assembly
#  iput $assembly/*.BigBed
#done < $prepData/merged_assemblies.txt

## initiate a given track hub
hub_name=$"HorseTrans1"
shortlabel=$"TopCuff_Cuffmerge"
longlabel=$"Single samlpe reference guided Tophat/Cufflinks followed by Cuffmerge"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$prepData/merged_assemblies.txt
tiss_assemblies=$tissue_Cuffmerge/tissue_assemblies.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

## create the HTML file page for every track

## add metadata like closest Ref gene
grep "exon_number \"1\"" merged.gtf > merged_ex1.gtf
grep "class_code \"u\"" merged_ex1.gtf > merged_ex1_u.gtf
grep -v "class_code \"u\"" merged_ex1.gtf > merged_ex1_nu.gtf

#######################
## run Transdecoder to predict UTRs with homology options
sample_list="$horse_trans/merged_and_tissue_assemblies.txt"
bash $script_path/run_transdecoder.sh $sample_list $genome $refPtn $refPfam $script_path/transdecoder.sh

## calculate the phase of Transdecoder GFF3 files
#while read assembly; do if [ -f $assembly/transdecoder/transcripts.fasta.transdecoder.genome.gff3 ];then
#  echo $assembly
#  cd $assembly/transdecoder
#  bash $script_path/cdsphase.sh transcripts.fasta.transdecoder.genome.gff3
#fi; done < $horse_trans/merged_and_tissue_assemblies.txt
#######################
## create list of assemblies from each library
## This is where you can edit the list to restrict the processing for certain target(s)
rm -f $prepData/merged_decoder_assemblies.txt
while read work_dir; do
  for dir in $work_dir/tophat_output/cuffmerge_output/withoutGTF/transdecoder; do if [ -d $dir ]; then
    echo ${dir#$prepData/} >> $prepData/merged_decoder_assemblies.txt;
  fi; done;
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## create list of assemblies for tissues of multiple libraries
rm -f $tissue_Cuffmerge/tissue_decoder_assemblies.txt
for tissue in $tissue_Cuffmerge/*; do
  echo ${tissue#$tissue_Cuffmerge/}/withoutGTF/transdecoder >> $tissue_Cuffmerge/tissue_decoder_assemblies.txt;
done

## convert the bed files into BigBed files & copy to the track hub directory & create assembly list
rm -f $horse_trans/merged_and_tissue_decoder_assemblies.txt
while read assembly; do
  echo $assembly
  cd $prepData/$assembly
  targetAss=$"transcripts.fasta.transdecoder.genome.bed"
  bash $script_path/bedToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes"
  if [ -f $"transcripts.fasta.transdecoder.genome.BigBed" ];then
    identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    cp transcripts.fasta.transdecoder.genome.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    echo $prepData/$assembly >> $horse_trans/merged_and_tissue_decoder_assemblies.txt;
fi; done < $prepData/merged_decoder_assemblies.txt

while read assembly; do
  echo $assembly
  cd $tissue_Cuffmerge/$assembly
  targetAss=$"transcripts.fasta.transdecoder.genome.bed"
  bash $script_path/bedToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes"
  if [ -f $"transcripts.fasta.transdecoder.genome.BigBed" ];then
  identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
  cp transcripts.fasta.transdecoder.genome.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
  echo $tissue_Cuffmerge/$assembly >> $horse_trans/merged_and_tissue_decoder_assemblies.txt;
fi; done < $tissue_Cuffmerge/tissue_decoder_assemblies.txt

#########################
## initiate a given track hub
hub_name=$"HorseTrans2"
shortlabel=$"TopCuffmerge_TD"
longlabel=$"Single samlpe ref guided Tophat/Cufflinks followed by Cuffmerge/transdecoder"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$prepData/merged_decoder_assemblies.txt
tiss_assemblies=$tissue_Cuffmerge/tissue_decoder_assemblies.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

###########################################################################################

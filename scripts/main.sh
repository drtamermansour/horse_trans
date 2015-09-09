myRoot=$"/mnt/ls12/Tamer"
horse_trans=$myRoot/horse_trans
script_path=$myRoot/horse_trans/scripts
resources=$myRoot/horse_trans/resources
rawData=$myRoot/horse_trans/rawdata
prepData=$myRoot/horse_trans/prepdata
genome_dir=$myRoot/horse_trans/refGenome
track_hub=$myRoot/horse_trans/track_hub
###########################################################################################
## assess_public_annotations
bash $script_path/assess_public_annotations.sh $resources

###########################################################################################
#### prepare the raw fastq files:
## Every tissue has a separate folder carrying its name (maximum 14 letter) in $prepData.
## Every libarary should have a separate folder into the corresponding tissue folder.
## The libarary folder name should start with PE_ or SE_ according to the sequencoing type
## Then it should should have the read length followed by underscore
## Then it should have fr.unstranded_ , fr.firststrand_ , fr.secondstrand_ according to lib type
## Then the owner name (or names separted by dots) followed by underscore
## Then the date of sequencing as MMDDYYYY
## The raw data files should be prepared so that they have enconding "Sanger / illumina 1.9"
## all sample replicates should be merged into one sample
## The file names should fit the format *_R1_*.fastq.gz & *_R2_*.fastq.gz for PE reads or *_SR_*.fastq.gz for SE
## The first syllbus (the part of the file name before _R1_ , _R2_ or _SR_) should be unique
## the final form of the data files should be kept in a folder named fastq_data in the tissue folder
bash $script_path/prep_in_home_seq_files.sh $rawData $prepData $script_path

proj_rawData=$rawData/PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014
SRA_URL=$"ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653"
bash $script_path/prep_proj.sh $proj_rawData $SRA_URL $prepData $script_path
###########################################################################################
## define the list of working directory and the list samples in each. This is where you can edit the output list file(s) to restrict the processing for certain target(s)
rm -f $horse_trans/working_list.txt
for work_dir in $prepData/*/{PE_*,SE_*}; do if [ -d $work_dir/fastq_data ]; then
  echo $work_dir >> $horse_trans/working_list.txt
  rm -f $work_dir/fastq_data/sample_list.txt
  for f in $work_dir/fastq_data/{*_R1_*.fastq.gz,*_SR_*.fastq.gz}; do if [ -f $f ]; then
    echo $f >> $work_dir/fastq_data/sample_list.txt; fi; done;
fi; done;

#### trimming with sliding window
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/trimmed_RNA_reads
  cd $work_dir/trimmed_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                   ## PE or SE
  sample_list=$work_dir/fastq_data/sample_list.txt
  bash ${script_path}/run_adapter_trimmer.sh $sample_list $lib $script_path
done < $horse_trans/working_list.txt

## Check for successful trimming and trouble shooting the failed trimming jobs (requires T_Trim.e)
while read work_dir; do
  cd $work_dir/trimmed_RNA_reads
  sample_list=$work_dir/fastq_data/failedSamples.txt           ## define the path of empty file
  bash $script_path/check_ad_trim.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    echo "Failed jobs in: "$work_dir
    bash ${script_path}/run_adapter_trimmer.sh $sample_list $lib $script_path
  fi
done < $horse_trans/working_list.txt

#### merge singletones
while read work_dir; do if [ -d $work_dir/trimmed_RNA_reads ]; then
  cd $work_dir/trimmed_RNA_reads
  ## change /2 to /1 in s2_se then combine single reads
  for f in $work_dir/trimmed_RNA_reads/*_R1_*.se.fq; do if [ -f $f ]; then echo $f; f2=$(basename $f | sed 's/_R1_/_R2_/'); newf=$(basename $f | sed 's/_R1_/_R_/'); sed 's/\/2$/\/1/g' $f2 > $f2.temp; cat $f $f2.temp > $newf; fi; done;
  rm -f *.se.fq.temp
  ## merge the single reads to the end of s1_pe file to make s1_pe_se
  for f in $work_dir/trimmed_RNA_reads/*_R1_*.pe.fq; do if [ -f $f ]; then echo $f; fr=$(basename $f | sed 's/_R1_/_R_/'); fr2=$(echo $fr | sed 's/.pe.fq/.se.fq/'); newf=$(basename $f | sed 's/.pe.fq/.pe.se.fq/'); cat $f $fr2 > $newf; fi; done;
  #rm {*_R1_*.pe.fq,*.se.fq}   ## all what you need *_R1_*.pe.se.fq and *_R2_*.pe.fq
fi
done < $horse_trans/working_list.txt
###########################################################################################
## get the referenece genome and prepare Bowtie2Index
cd $genome_dir
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab2/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
tar xvzf chromFa.tar.gz
mkdir -p $genome_dir/Bowtie2Index && cd $genome_dir/Bowtie2Index
cat ../*.fa > genome.fa
module load bowtie2/2.1.0
bowtie2-build genome.fa genome
Bowtie2_genome_index_base=$genome_dir/Bowtie2Index/genome

## define the GTF/GFF files
## generation of GTF from UCSC tables using the guidelines of genomewiki.ucsc
## http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
## The example on the wiki is based ok knowngene table which does not exist for horses. Instead there are refGene and ensGene tables
## The commands of wikigenome are modified to match the table schemes
#$HOME/bin/UCSC_kent_commands/genePredToGtf equCab2 refGene refGene.gtf    ## Couldn't connect to database equCab2 on genome-mysql.cse.ucsc.edu as genomep
wget http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/refGene.txt.gz
zcat refGene.txt.gz | cut -f2-11 | $HOME/bin/UCSC_kent_commands/genePredToGtf file stdin refGene.gtf
Genes_GTF_file=$genome_dir/refGene.gtf
#Genes_GFF_file=$resources/NCBI/GFF/ref_EquCab2.0_top_level.gff3

###########################################################################################
## pipeline_OneSampleAtaTime_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
## build the transcriptome-index
transcriptome_index=$genome_dir/trans_index/equ
module load TopHat2/2.0.14              ## bowtie2/2.1.0 => bowtie2/2.2.3
tophat --GTF "${Genes_GTF_file}" --transcriptome-index "${transcriptome_index}" "${Bowtie2_genome_index_base}"
rm -R tophat_out

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
  bash ${script_path}/run_cufflinks.sh "$sample_list" "$Genes_GTF_file" "$script_path/cufflinks.sh";
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
    bash ${script_path}/run_cufflinks.sh "$sample_list" "$Genes_GTF_file" "$script_path/cufflinks.sh";
fi; done < $horse_trans/working_list_NoPBMCs.txt
############
## Assess computational utilization of cufflinks
cufflinks_utlization=$prepData/cufflinks_utlization.txt
> $cufflinks_utlization
while read work_dir; do
  cd $work_dir/tophat_output
  echo $(pwd) >> $cufflinks_utlization
  for dir in tophat_*; do
    f=$dir/cufflinks.e*
    if [ -f $f ]; then
      f2=$(echo $f | sed 's/cufflinks.e/cufflinks.o/')
      echo $dir >> $cufflinks_utlization
      grep "resources_used.vmem =" $f2 >> $cufflinks_utlization
      grep "resources_used.walltime =" $f2 >> $cufflinks_utlization
      grep "Resource_List.nodes =" $f2 >> $cufflinks_utlization
    else echo "no cufflinks report in $dir" >> $cufflinks_utlization
fi; done; done < $horse_trans/working_list_NoPBMCs.txt
#############
waiting for
27040191.mgr-04.i       mansourt    main     cufflinks         69153     1     16   48gb 168:00:00 R  38:24:45
27040192.mgr-04.i       mansourt    main     cufflinks         74838     1     16   48gb 168:00:00 R  37:53:49
tophat_6_MO_C_AGTCAA
tophat_Nancy_TGACCA
#############


##################
### Run cuffmerge: merge the sample assemblies and output merged.gtf in tophat_output/cuffmerge_output
module load cufflinks/2.2.1
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cd $work_dir/tophat_output
  for dir in $work_dir/tophat_output/tophat_*; do if [ -d "$dir" ]; then
    echo "$dir"/transcripts.gtf; fi; done > assemblies.txt
  ## with the --ref-gtf option, merge input assemblies together with the reference GTF
  #cuffmerge -s $genome_dir/Bowtie2Index/genome.fa -o cuffmerge_output/withRefGene --num-threads 4 --ref-gtf "${Genes_GTF_file}" assemblies.txt > cuffmerge_withRef.txt 2>&1
  ## without the --ref-gtf option
  cuffmerge -s $genome_dir/Bowtie2Index/genome.fa -o cuffmerge_output/withoutGTF --num-threads 4 assemblies.txt > cuffmerge_withoutGTF.txt 2>&1
fi; done < $horse_trans/working_list_NoPBMCs_NoCereb.txt
## http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#transfrag-class-codes
##################
### cuffmerge the lab-specific tissue assemblies into tissue specific assemblies
# prepare list of target tissues
rm -f $horse_trans/multi_lib_tissues.txt
for tissue_dir in $prepData/*; do if [[ -d $tissue_dir && $(ls $tissue_dir | wc -l) -gt 1 ]]; then
  echo $tissue_dir >> $horse_trans/multi_lib_tissues.txt; fi; done;

# For every target tissue, cuffmerge the lab-spefific assemblies into one assembly in $horse_trans/tissue_merge/cuffmerge/TISSUE.NAME/withORwithoutRefGuidence.gtf
tissue_Cuffmerge=$horse_trans/tissue_merge/cuffmerge
while read tissue_dir; do
  tissue=$(basename $tissue_dir)
  mkdir -p $tissue_Cuffmerge/$tissue
  rm -f $tissue_dir/All_libraries_assemblies.txt
  for dir in $tissue_dir/*; do if [ -d "$dir" ]; then
    cat "$dir"/tophat_output/assemblies.txt >> $tissue_dir/All_libraries_assemblies.txt; fi; done

  ## with the --ref-gtf option, merge input assemblies together with the reference GTF
  #cuffmerge -s $genome_dir/Bowtie2Index/genome.fa -o $tissue_Cuffmerge/$tissue/withRefGene --num-threads 4 --ref-gtf "${Genes_GTF_file}" $tissue_dir/All_libraries_assemblies.txt > $tissue_merge/cuffmerge_withRef.txt 2>&1
  ## without the --ref-gtf option
  cuffmerge -s $genome_dir/Bowtie2Index/genome.fa -o $tissue_Cuffmerge/$tissue/withoutGTF --num-threads 4 $tissue_dir/All_libraries_assemblies.txt > $tissue_merge/cuffmerge_withoutGTF.txt 2>&1
done < $horse_trans/multi_lib_tissues.txt
##################
### cuffmerge all assemblies into one total assembly
# save the assembly in $horse_trans/total_merge/cuffmerge/withORwithoutRefGuidence.gtf
rm -f $prepData/All_tissues_assemblies.txt
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cat "$work_dir"/tophat_output/assemblies.txt >> $prepData/All_tissues_assemblies.txt
done < $horse_trans/working_list.txt

mkdir -p $tissue_Cuffmerge/all_tissues
## with the --ref-gtf option, merge input assemblies together with the reference GTF
#cuffmerge -s $genome_dir/Bowtie2Index/genome.fa -o $tissue_Cuffmerge/all_tissues/withRefGene --num-threads 4 --ref-gtf "${Genes_GTF_file}" $prepData/All_tissues_assemblies.txt > $Alltissue_merge/cuffmerge_withRef.txt 2>&1
## without the --ref-gtf option
cuffmerge -s $genome_dir/Bowtie2Index/genome.fa -o $tissue_Cuffmerge/all_tissues/withoutGTF --num-threads 4 $prepData/All_tissues_assemblies.txt > $Alltissue_merge/cuffmerge_withoutGTF.txt 2>&1
####################
## predict UTR
## Add some code to do the UTR prediction
####################
#### Format the data (create the BigBed files)

## fetch the UCSC database to get the chromosome sizes
UCSCgenome=$"equCab2"
module load ucscUtils/262
fetchChromSizes $UCSCgenome > $genome_dir/$UCSCgenome.chrom.sizes
## Create the basic directory structure of the track hubs
mkdir -p $track_hub/$UCSCgenome/BigBed

## initiate a given track hub
hub_name=$"HorseTrans1"
shortlabel=$"TopCuff_Cuffmerge"
longlabel=$"Single samlpe reference guided Tophat/Cufflinks followed by Cuffmerge"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh $UCSCgenome $hub_name $shortlabel $longlabel $email

## create list of assemblies from each library
## This is where you can edit the list to restrict the processing for certain target(s)
rm -f $prepData/merged_assemblies.txt
while read work_dir; do if [ -d $work_dir/tophat_output/cuffmerge_output ]; then
  for dir in $work_dir/tophat_output/cuffmerge_output/withoutGTF; do if [ -d $dir ]; then
    echo ${dir#$prepData/} >> $prepData/merged_assemblies.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs_NoCereb.txt
## convert the gtf files into BigBed files
while read assembly; do
  echo $assembly
  cd $prepData/$assembly && rm -f *.BigBed
  bash $script_path/gtfToBigBed.sh "$genome_dir/$UCSCgenome.chrom.sizes" "$assembly" "$script_path"
done < $prepData/merged_assemblies.txt
## copy the BigBed files to the track hub directory
while read assembly; do
  echo $assembly
  cp $prepData/$assembly/*.BigBed $track_hub/$UCSCgenome/BigBed/.
done < $prepData/merged_assemblies.txt

## create list of assemblies for tissues of multiple libraries
rm -f $tissue_Cuffmerge/tissue_assemblies.txt
for tissue in $tissue_Cuffmerge/*; do
  echo ${tissue#$tissue_Cuffmerge/}/withoutGTF >> $tissue_Cuffmerge/tissue_assemblies.txt; done
## convert the gtf files into BigBed files
while read assembly; do
  echo $assembly
  cd $tissue_Cuffmerge/$assembly && rm -f *.BigBed
  bash $script_path/gtfToBigBed.sh "$genome_dir/$UCSCgenome.chrom.sizes" "$assembly" "$script_path"
done < $tissue_Cuffmerge/tissue_assemblies.txt
## copy the BigBed files to the track hub directory
while read assembly; do
  echo $assembly
  cp $tissue_Cuffmerge/$assembly/*.BigBed $track_hub/$UCSCgenome/BigBed/.
done < $tissue_Cuffmerge/tissue_assemblies.txt

## run icommand to push the file to iplant
## https://pods.iplantcollaborative.org/wiki/display/DS/Using+iCommands
## http://bioinformatics.plantbiology.msu.edu/display/IP/Moving+Data+from+HPCC+to+iPlant
#icd /iplant/home/drtamermansour/horseTrans
#while read assembly; do
#  echo $assembly
#  iput $assembly/*.BigBed
#done < $prepData/merged_assemblies.txt

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$prepData/merged_assemblies.txt
tiss_assemblies=$tissue_Cuffmerge/tissue_assemblies.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

## create the HTML file page for every track

## add metadata like closest Ref gene
grep "exon_number \"1\"" merged.gtf > merged_ex1.gtf
grep "class_code \"u\"" merged_ex1.gtf > merged_ex1_u.gtf
grep -v "class_code \"u\"" merged_ex1.gtf > merged_ex1_nu.gtf

###########################################################################################
###########################################################################################
###########################################################################################

## pipeline_diginormAllsamples_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
cd ~/khmer
git checkout horseTrans         ## the branch was created on 08/24/2015 for reproducibility
## interleave PE files
while read work_dir; do
  echo $work_dir
  cd $work_dir/trimmed_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    bash ${script_path}/run_interleave.sh "$sample_list" $script_path/interleave.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful interleave runs and trouble shooting the failed jobs (requires interleave.e)
while read work_dir; do
  cd $work_dir/trimmed_RNA_reads
  sample_list=$work_dir/trimmed_RNA_reads/interleave_failedSamples.txt           ## define the path of empty file
  bash $script_path/check_interleave.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    echo "Failed interleave jobs in: "$work_dir
    cat $sample_list
    bash ${script_path}/run_interleave.sh "$sample_list" $script_path/interleave.sh;
fi
done < $horse_trans/working_list_NoPBMCs.txt

## run digital normalization of lab specific tissues (need to be updated to use sample list and check for success)
kmer=20
cutoff=50
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/normalizied_RNA_reads
  cd $work_dir/normalizied_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/trimmed_RNA_reads" "$kmer" "$cutoff" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful diginorm and trouble shooting the failed jobs (requires T_Trim.e)
sample_list=$prepData/failed_diginorm.txt                   ## define the path of empty file
> $sample_list
while read work_dir; do
  cd $work_dir/normalizied_RNA_reads
  bash $script_path/check_diginorm.sh "$sample_list"
done < $horse_trans/working_list_NoPBMCs.txt
x=$(cat $sample_list | wc -l)
if [ $x -ne 0 ]; then
  echo "Failed jobs in: "
  cat $sample_list
  while read work_dir; do
    cd $work_dir/normalizied_RNA_reads
    lib=$(basename $work_dir | cut -d"_" -f 1)
    bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/trimmed_RNA_reads" "$kmer" "$cutoff" "$script_path"
  done < $sample_list
fi

## split the interleaved reads
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    bash ${script_path}/run_split_reads.sh "$sample_list" $script_path/split_reads.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge singletons and change the file names to fit the tophat script
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    singletones=1
    for f in *_R1_001.pe.fq; do
      base=${f%_R1_001.pe.fq}
      if [ $singletones -eq 1 ]; then cat $f allsingletons.fq.keep > "$base"_R1_001.pe.se.fq; singletones=0;
      else mv $f "$base"_R1_001.pe.se.fq; fi; done
  elif [ "$lib" = $"SE" ]; then
    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge the files in tophat compatible format
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/normalizied_RNA_reads
#  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
#  if [ "$lib" = $"PE" ]; then
#    cat *_R1_001.pe.fq allsingletons.fq.keep > allsamples_R1_002.pe.se.fq
#    cat *_R2_001.pe.fq > allsamples_R2_002.pe.fq
#  elif [ "$lib" = $"SE" ]; then
#    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
#fi; done < $horse_trans/working_list_NoPBMCs.txt

## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/normalizied_RNA_reads ]; then
  rm -f $work_dir/normalizied_RNA_reads/sample_list.txt
  for f in $work_dir/normalizied_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
    echo $f >> $work_dir/normalizied_RNA_reads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## run Tophat on each sample independantly
## we meight need to run tophat on one sample after another feeding each run by the junctions of the previous run
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/digi_tophat_output
  cd $work_dir/digi_tophat_output
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  sample_list=$work_dir/normalizied_RNA_reads/sample_list.txt
  bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful tophat runs and trouble shooting the failed tophat jobs (require tophat-[SP]E.e & .o)
while read work_dir; do
  cd $work_dir/digi_tophat_output
  sample_list=$work_dir/normalizied_RNA_reads/tophat_failedSamples.txt           ## define the path of empty file
  bash $script_path/check_tophat.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')
    echo "Failed tophat jobs in: "$work_dir
    bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
fi; done < $horse_trans/working_list_NoPBMCs.txt

##################
## create summary for tophat run
while read work_dir; do
  for f in $work_dir/digi_tophat_output/tophat_*; do
    echo ${f} >> $work_dir/digi_tophat_output/allsample_summary.txt
    cd ${f}
    grep "overall read mapping rate" align_summary.txt >> ../allsample_summary.txt
    grep "concordant pair alignment rate" align_summary.txt >> ../allsample_summary.txt
  done
done < $horse_trans/working_list_NoPBMCs.txt
##################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/digi_tophat_output ]; then
  rm -f $work_dir/digi_tophat_output/sample_list.txt
  for f in $work_dir/digi_tophat_output/tophat_*; do if [ -d $f ]; then
    echo $f >> $work_dir/digi_tophat_output/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## Merge BAM files
module load SAMTools/0.1.19
while read work_dir; do
  echo $work_dir
  cd $work_dir/digi_tophat_output
  samples=()
  while read sample; do samples+=($sample/accepted_hits.bam); done < sample_list.txt
  len=${#samples[@]}
  if [ $len -gt 1 ]; then samtools merge merged.bam ${samples[*]};
  elif [ $len -eq 1 ]; then cp ${samples[0]} merged.bam;
  else echo "can find bam files"; fi
done < $horse_trans/working_list_NoPBMCs.txt

### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
while read work_dir; do
  echo $work_dir
  cd $work_dir/digi_tophat_output
  sample=$work_dir/digi_tophat_output/merged.bam
  label=$(basename $work_dir)
  bash ${script_path}/run_cufflinks2.sh "$sample" "$Genes_GTF_file" "$label" "$script_path/cufflinks2.sh";
done < $horse_trans/working_list_NoPBMCs.txt

### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/digi_tophat_output
#  sample_list=$work_dir/digi_tophat_output/sample_list.txt
#  bash ${script_path}/run_cufflinks.sh "$sample_list" "$Genes_GTF_file" "$script_path/cufflinks2.sh";
#done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful Cufflinks runs and trouble shooting the failed Cufflinks jobs (requires cufflinks.e)
#while read work_dir; do
#  cd $work_dir/digi_tophat_output
#  sample_list=$work_dir/digi_tophat_output/failedSamples.txt           ## define the path of empty file
#  bash $script_path/check_cufflinks.sh "$sample_list"
#  x=$(cat $sample_list | wc -l)
#  if [ $x -ne 0 ]; then
#    echo "Failed Cufflinks jobs in: "$work_dir
#    cat $sample_list
#     bash ${script_path}/run_cufflinks.sh "$sample_list" "$Genes_GTF_file" "$script_path/cufflinks.sh";
#fi; done < $horse_trans/working_list_NoPBMCs.txt
###################
## Assess computational utilization of cufflinks
#cufflinks_utlization=$prepData/digi_cufflinks_utlization.txt
#> $cufflinks_utlization
#while read work_dir; do
#  cd $work_dir/digi_tophat_output
#  echo $(pwd) >> $cufflinks_utlization
#  for dir in tophat_*; do
#    f=$dir/cufflinks.e*
#    if [ -f $f ]; then
#      f2=$(echo $f | sed 's/cufflinks.e/cufflinks.o/')
#      echo $dir >> $cufflinks_utlization
#      grep "resources_used.vmem =" $f2 >> $cufflinks_utlization
#      grep "resources_used.walltime =" $f2 >> $cufflinks_utlization
#      grep "Resource_List.nodes =" $f2 >> $cufflinks_utlization
#    else echo "no cufflinks report in $dir" >> $cufflinks_utlization
#fi; done; done < $horse_trans/working_list_NoPBMCs.txt

## multi-tissue diginorm
## start with tissues with longer reads (& higher mapping effeciency)

## initiate a given track hub
hub_name=$"HorseTrans2"
shortlabel=$"rawdigi_TopCuff"
longlabel=$"diginorm of raw data followed by reference guided Tophat/Cufflinks"
cd $track_hub
bash $script_path/create_trackHub.sh $UCSCgenome $hub_name $shortlabel $longlabel $email


###########################################################################################
##### Suggested piplines
### With Tophat
## Compare Ref-guided versus reference free
## Try very senstive Bowtie search
## Try --coverage-search
## For diginorm, try to start with mapped reads only (to minimize error accumulation ==> decrease RAM required to achive minimal false postive rate and allow either final mapping). You can diginorm the non-mapped reads separetly

###############################################################################
## convert the tophat output BAM files into Fastq files
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/mapped_trimmed
  cd $work_dir/mapped_trimmed
  sample_list=$work_dir/tophat_output/sample_list.txt
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_BamToFastq.sh "$sample_list" "$lib" "$script_path/restore_mapped_trimmed.sh"
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## Check for successful jobs (requires restore_mapped.e & restore_mapped.o)
while read work_dir; do
  cd $work_dir/mapped_trimmed
  sample_list=$work_dir/mapped_trimmed/failedSamples.txt           ## define the path of empty file
  bash $script_path/check_BamToFastq.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    echo "Failed jobs in: "$work_dir
    bash ${script_path}/run_BamToFastq.sh "$sample_list" "$lib" "$script_path/restore_mapped_trimmed.sh"
fi; done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## run digital normalization of lab specific tissues (need to be updated to use sample list and check for success)
kmer=20
cutoff=50
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/normalizied_mappedRNA_reads
  cd $work_dir/normalizied_mappedRNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/mapped_trimmed" "$kmer" "$cutoff" "$script_path"
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## Check for successful diginorm and trouble shooting the failed jobs (requires T_Trim.e)
sample_list=$prepData/failed_diginorm.txt                   ## define the path of empty file
> $sample_list
while read work_dir; do
  cd $work_dir/normalizied_mappedRNA_reads
  bash $script_path/check_diginorm.sh "$sample_list"
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt
x=$(cat $sample_list | wc -l)
if [ $x -ne 0 ]; then
    echo "Failed jobs in: "
    cat $sample_list
    while read work_dir; do
      cd $work_dir/normalizied_mappedRNA_reads
      lib=$(basename $work_dir | cut -d"_" -f 1)
      bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/mapped_trimmed" "$kmer" "$cutoff" "$script_path"
    done < $sample_list
fi

## split the interleaved reads
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_mappedRNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    bash ${script_path}/run_split_reads.sh "$sample_list" $script_path/split_reads.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge singletons and change the file names to fir the tophat script
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_mappedRNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    singletones=1
    for f in *_R1_001.pe.fq; do
      base=${f%_R1_001.pe.fq}
      if [ $singletones -eq 1 ]; then cat $f allsingletons.fq.keep > "$base"_R1_001.pe.se.fq; singletones=0;
      else mv $f "$base"_R1_001.pe.se.fq; fi; done
  elif [ "$lib" = $"SE" ]; then
    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
fi; done < $horse_trans/working_list_NoPBMCs.txt

## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/normalizied_RNA_reads ]; then
rm -f $work_dir/normalizied_RNA_reads/sample_list.txt
for f in $work_dir/normalizied_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
echo $f >> $work_dir/normalizied_RNA_reads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## run Tophat on each sample
while read work_dir; do
echo $work_dir
mkdir -p $work_dir/digi_tophat_output
cd $work_dir/digi_tophat_output
lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
sample_list=$work_dir/normalizied_RNA_reads/sample_list.txt
bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt














while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/unmapped_trimmed
  cd $work_dir/unmapped_trimmed
  sample_list=$work_dir/tophat_output/sample_list.txt
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_BamToFastq.sh "$sample_list" "$lib" "$script_path/restore_unmapped_trimmed.sh"
done < $horse_trans/working_list_NoPBMCs_NoCereb_NoSkin.txt



#################################################################################
## Split GTF and BAM files by chromosomes
cat ${Genes_GTF_file} | awk '{ print $1}' | sort | uniq > $genome_dir/refChrom_list.txt
mkdir $genome_dir/GTFbyChrom
while read refChrom; do
  cat ${Genes_GTF_file} | awk -v ref="$refChrom" '($1 == ref)' >> $genome_dir/GTFbyChrom/"$refChrom"_refGene.gtf
done < $genome_dir/refChrom_list.txt


module load BAMTools/2.2.3
bamtools split -in accepted_hits.bam -reference


module load cufflinks/2.2.1
cufflinks --GTF-guide $genome_dir/GTFbyChrom/chr10_refGene.gtf --num-threads 4 --verbose accepted_hits.REF_chr10.bam






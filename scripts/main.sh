#!/bin/sh

## construction of the basic diretory structure
## you need to change myRoot to match you actual working directory
myRoot=$"/mnt/ls12/Tamer"
mkdir -p $myRoot/horse_trans/{scripts,resources,rawdata,prepdata,tissue_merge,refGenome,public_assemblies,track_hub}

## create a config file to contain all the pathes to be used by all pipelines
echo "horse_trans=$myRoot/horse_trans" >> $myRoot/config.txt
echo "script_path=$myRoot/horse_trans/scripts" >> $myRoot/config.txt
echo "resources=$myRoot/horse_trans/resources" >> $myRoot/config.txt
echo "rawData=$myRoot/horse_trans/rawdata" >> $myRoot/config.txt
echo "prepData=$myRoot/horse_trans/prepdata" >> $myRoot/config.txt
echo "tissue_merge=$myRoot/horse_trans/tissue_merge" >> $myRoot/config.txt
echo "genome_dir=$myRoot/horse_trans/refGenome" >> $myRoot/config.txt
echo "pubAssemblies=$myRoot/horse_trans/public_assemblies" >> $myRoot/config.txt
echo "track_hub=$myRoot/horse_trans/track_hub" >> $myRoot/config.txt
source $myRoot/config.txt
###########################################################################################
#### prepare the raw fastq files:
## Every tissue has a separate folder carrying its name (maximum 14 letter) in $prepData.
## Every RNAseq libarary should have a separate folder into the corresponding tissue folder.
## The libarary folder name should start with PE_ or SE_ according to the sequencoing type
## Then it should should have the read length followed by underscore
## Then it should have fr.unstranded_ , fr.firststrand_ , fr.secondstrand_ according to lib type
## Then the owner name (or names separted by dots) followed by underscore
## Then the date of sequencing as MMDDYYYY
## The raw data files should be kept in a folder named fastq_data in the libarary folder
## The raw data files should be prepared so that they have enconding "Sanger / illumina 1.9"
## The file names should fit the format *_R1_*.fastq.gz & *_R2_*.fastq.gz for PE reads or *_SR_*.fastq.gz for SE
## All reads in every given file should belong to ONE sequencing lane.
## All sample replicates from the same lane should be merged into one sample
## If there are sample replicates from different lanes, you can add a text file called "replicates.txt"
## A line in this file should have the names of one sample replicates with space separation. (only the *_R1_*.fastq.gz for PE reads or *_SR_*.fastq.gz for SE)
## The first syllbus (the part of the file name before _R1_ , _R2_ or _SR_) should be unique

## a representive work out to one SRA reposatories
newLib=$"PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014"
SRA_URL=$"ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653"
mkdir -p $rawData/$newLib
cd $rawData/PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014
bash $script_path/prep_proj.sh $SRA_URL $prepData $script_path ## the script assume no sample replicates
###########################################################################################
## define the list of working directory and the list samples in each. This is where you can edit the output list file(s) to restrict the processing for certain target(s)
rm -f $horse_trans/working_list.txt
for work_dir in $prepData/*/{PE_*,SE_*}; do if [ -d $work_dir/fastq_data ]; then
  echo $work_dir >> $horse_trans/working_list.txt
  rm -f $work_dir/fastq_data/sample_list.txt
  for f in $work_dir/fastq_data/{*_R1_*.fastq.gz,*_SR_*.fastq.gz}; do if [ -f $f ]; then
    echo $f >> $work_dir/fastq_data/sample_list.txt; fi; done;
fi; done;

## get read length statistics
headers=$(Rscript -e 'cat("Tissue", "Library", "No_of_samples", "No_of_reads(M)", "Total_bp_count(Gb)", "Min_length", "Max_length", "Median_length", "Mean_length", sep="\t");')
echo "$headers" > $horse_trans/raw_statistics.txt
while read work_dir; do
  echo $work_dir
  cd $work_dir/fastq_data
  sample_list=$work_dir/fastq_data/sample_list.txt
  cat sample_list.txt | xargs zcat | awk '{if(NR%4==2) print length($1)}' > input.readslength.txt
  stat=$(Rscript -e 'd<-scan("input.readslength.txt", quiet=TRUE); cat(round(length(d)/10^6,2), round(sum(d)/10^9,2), min(d), max(d), median(d), mean(d), sep="\t");')
  lib=$(basename $work_dir)
  tissue=$(dirname $work_dir | xargs basename)
  no=$(cat $sample_list | wc -l)
  echo "$tissue"$'\t'"$lib"$'\t'"$no"$'\t'"$stat" >> $horse_trans/raw_statistics.txt
done < $horse_trans/working_list_NoPBMCs.txt

###########################################################################################
## prepare sorted working list accoding to read length
tail -n+2 $horse_trans/raw_statistics.txt | sort -k8,8nr | awk -v myRoot=$prepData '{ print myRoot"/"$1"/"$2 }' > $horse_trans/working_list_sorted.txt
###########################################################################################
#### Read trimming
## This step requires input working_list & sample_list and output folder
## the module should produce trimmed reads in the form of (*_R1_*.pe.se.fq and *_R2_*.pe.fq) for PE lib and (*_SR_*.se.fq) for SE
## The output is an sample list with *_R1_*.pe.se.fq and *_SR_*.se.fq files

## Mild trimming with Trimmomatic using sliding window
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

## define the list samples for subsequent analysis
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/trimmed_RNA_reads ]; then
  rm -f $work_dir/trimmed_RNA_reads/sample_list.txt
  for f in $work_dir/trimmed_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
    echo $f >> $work_dir/trimmed_RNA_reads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

###########################################################################################
## Assess read length statistics after trimming
headers=$(Rscript -e 'cat("Tissue", "Library", "No_of_samples", "No_of_reads(M)", "Total_bp_count(Gb)", "Min_length", "Max_length", "Median_length", "Mean_length", sep="\t");')
echo "$headers" > $horse_trans/trimmed_statistics.txt
while read work_dir; do
  echo $work_dir
  cd $work_dir/trimmed_RNA_reads
  sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
  lib=$(basename $work_dir)
  tissue=$(dirname $work_dir | xargs basename)
  no=$(cat $sample_list | wc -l)
  libType=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$libType" = $"PE" ]; then
    ## stats of the R1 files
    cat sample_list.txt | sed 's/.pe.se.fq/.pe.fq/'| xargs cat | awk '{if(NR%4==2) print length($1)}' > trimmed_R1.readslength.txt
    stat=$(Rscript -e 'd<-scan("trimmed_R1.readslength.txt", quiet=TRUE); cat(round(length(d)/10^6,2), round(sum(d)/10^9,2), min(d), max(d), median(d), mean(d), sep="\t");')
    echo "$tissue"$'\t'"$lib"$'_R1\t'"$no"$'\t'"$stat" >> $horse_trans/trimmed_statistics.txt

    ## stats of the R2 files
    cat sample_list.txt | sed 's/\(.*\)_R1_\(.*\).pe.se.fq/\1_R2_\2.pe.fq/'| xargs cat | awk '{if(NR%4==2) print length($1)}' > trimmed_R2.readslength.txt
    stat=$(Rscript -e 'd<-scan("trimmed_R2.readslength.txt", quiet=TRUE); cat(round(length(d)/10^6,2), round(sum(d)/10^9,2), min(d), max(d), median(d), mean(d), sep="\t");')
    echo "$tissue"$'\t'"$lib"$'_R2\t'"$no"$'\t'"$stat" >> $horse_trans/trimmed_statistics.txt

   ## stats of the singletones files
    cat sample_list.txt | sed 's/\(.*\)_R1_\(.*\).pe.se.fq/\1_R_\2.se.fq/'| xargs cat | awk '{if(NR%4==2) print length($1)}' > trimmed_SE.readslength.txt
    stat=$(Rscript -e 'd<-scan("trimmed_SE.readslength.txt", quiet=TRUE); cat(round(length(d)/10^6,2), round(sum(d)/10^9,2), min(d), max(d), median(d), mean(d), sep="\t");')
    echo "$tissue"$'\t'"$lib"$'_SE\t'"$no"$'\t'"$stat" >> $horse_trans/trimmed_statistics.txt
  elif [ "$libType" = $"SE" ]; then
    cat sample_list.txt | xargs cat | awk '{if(NR%4==2) print length($1)}' > trimmed.readslength.txt
    stat=$(Rscript -e 'd<-scan("trimmed.readslength.txt", quiet=TRUE); cat(round(length(d)/10^6,2), round(sum(d)/10^9,2), min(d), max(d), median(d), mean(d), sep="\t");')
    echo "$tissue"$'\t'"$lib"$'\t'"$no"$'\t'"$stat" >> $horse_trans/trimmed_statistics.txt
fi; done < $horse_trans/working_list_NoPBMCs.txt
###########################################################################################
## get the referenece genome
cd $genome_dir
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab2/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
tar xvzf chromFa.tar.gz

## prepare Bowtie2Index (for Tophat mapping)
mkdir -p $genome_dir/Bowtie2Index && cd $genome_dir/Bowtie2Index
cat ../*.fa > genome.fa
bash ${script_path}/run_bowtie2-build.sh genome.fa genome
echo "genome=$genome_dir/Bowtie2Index/genome.fa" >> $myRoot/config.txt
echo "Bowtie2_genome_index_base=$genome_dir/Bowtie2Index/genome" >> $myRoot/config.txt
source $myRoot/config.txt

## prepare BWA index (for GATK variant analysis)
mkdir -p $genome_dir/BwaIndex && cd $genome_dir/BwaIndex
cat ../*.fa > genome.fa
bash ${script_path}/run_bwa-index.sh genome.fa
echo "Bwa_ref=$genome_dir/BwaIndex/genome.fa" >> $myRoot/config.txt
source $myRoot/config.txt

## prepare GATK dictionary and index (for GATK variant analysis)
mkdir -p $genome_dir/gatkIndex && cd $genome_dir/gatkIndex
cat ../*.fa > genome.fa
bash ${script_path}/run_gatk-index.sh genome.fa
echo "gatk_ref=$genome_dir/gatkIndex/genome.fa" >> $myRoot/config.txt
echo "gatk_ref_index=$genome_dir/gatkIndex/genome.fa.fai" >> $myRoot/config.txt
source $myRoot/config.txt
###########################################################################################
## create liftover files
## http://genomewiki.ucsc.edu/index.php/LiftOver_Howto

# Download the genome files
mkdir $genome_dir/ncbi && cd $genome_dir/ncbi
wget -r --no-directories ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/Assembled_chromosomes/seq/eca_ref_EquCab2.0_*.fa.gz
gunzip eca_ref_EquCab2.0_*.fa.gz
cat eca_ref_EquCab2.0_*.fa > ncbi_genome.fa
## map the genomes
bash $script_path/mapGenome.sh              ## ends by creating ncbi/NCBItoUCSC_map.sorted.chain

###########################################################################################
## Create GTF file based of refGenes
## generation of GTF from UCSC tables using the guidelines of genomewiki.ucsc
## http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
## The example on the wiki is based ok knowngene table which does not exist for horses. Instead there are refGene and ensGene tables
## The commands of wikigenome are modified to match the table schemes
## Note: using genepredToGTF can resolve duplicates of transcript ids but not for gene ids so the 10 fields genepred format which uses transcript names as gene names produce no duplicates but the extended genepred uses separte gene names from column 12 and susbtable for gene name duplication
cd $genome_dir
wget http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/refGene.txt.gz
ucscTable=$"refGene.txt.gz"
output_GTF=$"refGene.gtf"
bash ${script_path}/ucscTableToGTF.sh $ucscTable $output_GTF
echo "refGTF_file=$genome_dir/refGene.gtf" >> $myRoot/config.txt
output_GTF=$"refGene_transcripts.gtf"
bash ${script_path}/ucscTableToGTF2.sh $ucscTable $output_GTF
echo "refTransGTF_file=$genome_dir/refGene_transcripts.gtf" >> $myRoot/config.txt
zcat $ucscTable | cut -f2-16 | $script_path/genePredToBed > refGene.bed
echo "refBED_file=$genome_dir/refGene.bed" >> $myRoot/config.txt
source $myRoot/config.txt

## Get the NCBI annotation files
wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/GFF/ref_EquCab2.0_top_level.gff3.gz
gunzip ref_EquCab2.0_top_level.gff3.gz
$HOME/bin/UCSC_kent_commands/gff3ToGenePred -useName ref_EquCab2.0_top_level.gff3 ref_EquCab2.0_top_level.gpred
## exclude non RNA entries e.g. CDs with no parant transcripts, gene_segments, ..
egrep "^rna|^NM|^NR|^XM|^XR" ref_EquCab2.0_top_level.gpred > ref_EquCab2.0_top_level_rna.gpred
$HOME/bin/UCSC_kent_commands/liftOver ref_EquCab2.0_top_level_rna.gpred ncbi/NCBItoUCSC_map.sorted.chain ref_EquCab2.0_top_level_mapped_rna.gpred unMapped -genePred
$HOME/bin/UCSC_kent_commands/genePredToGtf file ref_EquCab2.0_top_level_mapped_rna.gpred ref_EquCab2.0_top_level_rna.gtf
echo "ncbiGTF_file=$genome_dir/ref_EquCab2.0_top_level_rna.gtf" >> $myRoot/config.txt
cat ref_EquCab2.0_top_level_mapped_rna.gpred | $script_path/genePredToBed > ref_EquCab2.0_top_level_mapped_rna.bed
echo "ncbiBED_file=$genome_dir/ref_EquCab2.0_top_level_mapped_rna.bed" >> $myRoot/config.txt

$HOME/bin/UCSC_kent_commands/gff3ToGenePred ref_EquCab2.0_top_level.gff3 ref_EquCab2.0_top_level_noName.gpred
## exclude non RNA entries e.g. CDs with no parant transcripts, gene_segments, ..
grep "^rna" ref_EquCab2.0_top_level_noName.gpred > ref_EquCab2.0_top_level_noName_rna.gpred
$HOME/bin/UCSC_kent_commands/liftOver ref_EquCab2.0_top_level_noName_rna.gpred ncbi/NCBItoUCSC_map.sorted.chain ref_EquCab2.0_top_level_mapped_noName_rna.gpred unMapped_noName -genePred
$HOME/bin/UCSC_kent_commands/genePredToGtf file ref_EquCab2.0_top_level_mapped_noName_rna.gpred ref_EquCab2.0_top_level_noName_rna.gtf
echo "ncbiNoNameGTF_file=$genome_dir/ref_EquCab2.0_top_level_noName_rna.gtf" >> $myRoot/config.txt
cat ref_EquCab2.0_top_level_mapped_noName_rna.gpred | $script_path/genePredToBed > ref_EquCab2.0_top_level_mapped_noName_rna.bed
echo "ncbiNoNameBED_file=$genome_dir/ref_EquCab2.0_top_level_mapped_noName_rna.bed" >> $myRoot/config.txt
source $myRoot/config.txt

## Get the ensemble GTF files
wget http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/ensGene.txt.gz
ucscTable=$"ensGene.txt.gz"
output_GTF=$"ensGene.gtf"
bash ${script_path}/ucscTableToGTF.sh $ucscTable $output_GTF
echo "ensGTF_file=$genome_dir/ensGene.gtf" >> $myRoot/config.txt
zcat $ucscTable | cut -f2-16 | $script_path/genePredToBed > ensGene.bed
echo "ensBED_file=$genome_dir/ensGene.bed" >> $myRoot/config.txt
source $myRoot/config.txt
#wget ftp://ftp.ensembl.org/pub/release-80/gtf/equus_caballus/Equus_caballus.EquCab2.80.gtf.gz
#gunzip Equus_caballus.EquCab2.80.gtf.gz
#echo "ensGTF_file=$genome_dir/Equus_caballus.EquCab2.80.gtf" >> $myRoot/config.txt
#source $myRoot/config.txt
###########################################################################################
## get known variants
mkdir $genome_dir/knowVar
cd $genome_dir/knowVar
wget --timestamping 'ftp://ftp.ensembl.org/pub/release-82/variation/vcf/equus_caballus/Equus_caballus.vcf.gz' -O Equus_caballus.vcf.gz
wget --timestamping 'ftp://ftp.ensembl.org/pub/release-82/variation/vcf/equus_caballus/Equus_caballus_structural_variations.vcf.gz' -O Equus_caballus_structural_variations.vcf.gz
wget --timestamping 'ftp://ftp.ensembl.org/pub/release-82/variation/vcf/equus_caballus/Equus_caballus_incl_consequences.vcf.gz' -O Equus_caballus_incl_consequences.vcf.gz
gunzip *.gz

## change the name of the chromosomes to match the UCSC genome (bet keep the file co-ordinates 1-based)
grep -v "^#" Equus_caballus_structural_variations.vcf | awk -F "\t" -v OFS='\t' '{ print "chr"$1,$2,$3,$4,$5,$6,$7,$8 }' > Equus_caballus_structural_variations_fixedChrNames.vcf
sed -i 's/chrMT/chrM/g' Equus_caballus_structural_variations_fixedChrNames.vcf
perl $script_path/sortByRef.pl Equus_caballus_structural_variations_fixedChrNames.vcf $gatk_ref_index > Equus_caballus_structural_variations_fixedChrNames_sorted.vcf
grep "^#" Equus_caballus_structural_variations.vcf > Equus_caballus_structural_variations_final.vcf
cat Equus_caballus_structural_variations_fixedChrNames_sorted.vcf >> Equus_caballus_structural_variations_final.vcf
echo "knownIndels=$genome_dir/knowVar/Equus_caballus_structural_variations_final.vcf" >> $myRoot/config.txt
grep -v "^#" Equus_caballus.vcf | awk -F "\t" -v OFS='\t' '{ print "chr"$1,$2,$3,$4,$5,$6,$7,$8 }' > Equus_caballus_fixedChrNames.vcf
sed -i 's/chrMT/chrM/g' Equus_caballus_fixedChrNames.vcf
perl $script_path/sortByRef.pl Equus_caballus_fixedChrNames.vcf $gatk_ref_index > Equus_caballus_fixedChrNames_sorted.vcf
grep "^#" Equus_caballus.vcf > Equus_caballus_final.vcf
cat Equus_caballus_fixedChrNames_sorted.vcf >> Equus_caballus_final.vcf
echo "knownSNPs=$genome_dir/knowVar/Equus_caballus_final.vcf" >> $myRoot/config.txt
source $myRoot/config.txt
###########################################################################################
### Initiate the basic structure for horse track hubs
echo "UCSCgenome=equCab2" >> $myRoot/config.txt
source $myRoot/config.txt
## fetch the UCSC database to get the chromosome sizes
chromSizes=$genome_dir/$UCSCgenome.chrom.sizes
bash ${script_path}/calcChromSizes.sh $UCSCgenome $chromSizes
## Create the basic directory structure of the track hubs
mkdir -p $track_hub/$UCSCgenome/BigBed
###########################################################################################
### Install homology search databases
## download protein database such as Swissprot (fast) or Uniref90 (slow but more comprehensive)
mkdir -p $genome_dir/ptnDB
cd $genome_dir/ptnDB
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
#gunzip uniref90.fasta.gz
echo "refPtn=$genome_dir/ptnDB/uniprot_sprot.fasta" >> $myRoot/config.txt
#echo "refPtn=$genome_dir/ptnDB/uniref90.fasta" >> $myRoot/config.txt
source $myRoot/config.txt
bash $script_path/make_ptnDB.sh $refPtn
#bash $script_path/make_ptnDB.sh $refPtn

## download pfam database
wget ftp://ftp.broadinstitute.org/pub/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
echo "refPfam=$genome_dir/ptnDB/Pfam-A.hmm" >> $myRoot/config.txt
source $myRoot/config.txt
bash $script_path/make_PfamDB.sh $refPfam
###########################################################################################
## build Tophat transcriptome-index
echo "transcriptome_index=$genome_dir/trans_index/equ" >> $myRoot/config.txt
source $myRoot/config.txt
bash ${script_path}/buildTransIndex.sh "$refGTF_file" "$transcriptome_index" "$Bowtie2_genome_index_base"
###########################################################################################
#### Mapping
## This step requires input  working_list & sample_list and output folder
## the module should produce mapped reads in the form of tophat_<sampleID>/accepted.bam
## The output is an sample list with accepted.bam files

## refGene guided Tophat mapping per sample
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/tophat_output
  cd $work_dir/tophat_output
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
  bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"   ## Tophat sort by coordinate
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful tophat runs and trouble shooting the failed tophat jobs
while read work_dir; do
  cd $work_dir/tophat_output
  sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
  failedSample_list=$work_dir/trimmed_RNA_reads/tophat_failedSamples.txt
  > $failedSample_list                                                ## erase previouslly failed samples if any
  bash $script_path/check_tophat.sh "$failedSample_list"        ## require tophat-[SP]E.e & .o
  bash $script_path/check2_tophat.sh "$sample_list" "$failedSample_list"   ## check output log files
  x=$(cat $failedSample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')
    echo "Failed tophat jobs in: "$work_dir
    bash ${script_path}/run_tophat.sh "$failedSample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
  fi
done < $horse_trans/working_list_Cerebellum.txt
##################
## create summary for tophat run
headers=$(Rscript -e 'cat("Tissue", "Library", "min_mapping", "max_mapping", "min_concordance", "max_concordance", sep="\t");')
echo "$headers" > $horse_trans/tophat_summary.txt
while read work_dir; do
  > $work_dir/tophat_output/allsample_summary.txt
  for f in $work_dir/tophat_output/tophat_*; do
    echo ${f} >> $work_dir/tophat_output/allsample_summary.txt
    cd ${f}
    grep "overall read mapping rate" align_summary.txt >> ../allsample_summary.txt
    grep "concordant pair alignment rate" align_summary.txt >> ../allsample_summary.txt
  done
  mapping=$(grep "overall read mapping rate" $work_dir/tophat_output/allsample_summary.txt | awk '{ print $1 }' | sort -n | sed -e 1b -e '$!d' | tr "\n" "\t")
  conc=$(grep "concordant pair alignment rate" $work_dir/tophat_output/allsample_summary.txt | awk '{ print $1 }' | sort -n | sed -e 1b -e '$!d' | tr "\n" "\t")
  lib=$(basename $work_dir)
  tissue=$(dirname $work_dir | xargs basename)
  echo "$tissue"$'\t'"$lib"$'\t'"$mapping""$conc" >> $horse_trans/tophat_summary.txt
done < $horse_trans/working_list_NoPBMCs.txt
##################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  rm -f $work_dir/tophat_output/preMerge_sample_list.txt
  for f in $work_dir/tophat_output/tophat_*/accepted_hits.bam; do if [ -f $f ]; then
    echo $f >> $work_dir/tophat_output/preMerge_sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_Cerebellum.txt
###########################################################################################
## Add Read group headers
while read work_dir; do
  echo $work_dir
  lib=$(basename $work_dir)
  sample_list=$work_dir/tophat_output/preMerge_sample_list.txt
  bash ${script_path}/run_readGroupInfo_illumina.sh "$sample_list" "$lib" "${script_path}/readGroupInfo_illumina.sh"
done < $horse_trans/working_list_NoPBMCs_NoCereb_NoRetina.txt  ## done with cerebellum

## Check for successful adding of groupread
## To be added
## Tip: the line before last in .e file has "AddOrReplaceReadGroups done"
## for f in $prepData/*/*/tophat_output/tophat_*/AddOrReplaceReadGroups.e*; do grep "AddOrReplaceReadGroups done" $f | wc -l; done
###########################################################################################
## merge replicates
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  replicates_list=$work_dir/fastq_data/replicates.txt
  if [ -f $replicates_list ]; then
bash ${script_path}/run_mergeBAM.sh "$replicates_list" "${script_path}/mergeBAM.sh" ## picard tools sort by coordinate
  else echo "No replicates file in" $work_dir;
fi; done < $horse_trans/working_list_Cerebellum.txt

## check for successful merging
## To be added
##################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  rm -f $work_dir/tophat_output/sample_list.txt
  for f in $work_dir/tophat_output/tophat_*; do if [ -d $f ]; then
    echo $f >> $work_dir/tophat_output/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_Cerebellum.txt
###########################################################################################
#### pipeline_OneSampleAtaTime_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
#cufflinks_run="refGeneGuided_Cufflinks"
cufflinks_run="nonGuided_Cufflinks"
while read work_dir; do
  echo $work_dir
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  #bash ${script_path}/run_cufflinks.sh "$sample_list" "$refGTF_file" "$script_path/cufflinks.sh";
  bash ${script_path}/run_cufflinks_noRef.sh "$sample_list" "$script_path/cufflinks_noRef.sh";
done < $horse_trans/working_list_Cerebellum.txt    #working_list_NoPBMCs.txt

## Check for successful Cufflinks runs and trouble shooting the failed Cufflinks jobs (requires cufflinks.e)
while read work_dir; do
  cd $work_dir/tophat_output
  failedSample_list=$work_dir/tophat_output/failedSamples.txt           ## define the path of empty file
  bash $script_path/check_cufflinks.sh "$failedSample_list"
  x=$(cat $failedSample_list | wc -l)
  if [ $x -ne 0 ]; then
    echo "Failed Cufflinks jobs in: "$work_dir
    cat $failedSample_list
    #bash ${script_path}/run_cufflinks.sh "$failedSample_list" "$refGTF_file" "$script_path/cufflinks.sh";
#    bash ${script_path}/run_cufflinks_noRef.sh "$failedSample_list" "$refGTF_file" "$script_path/cufflinks_noRef2.sh";
fi; done < $horse_trans/working_list_NoPBMCs.txt
############
## Assess computational utilization of cufflinks
cufflinks_utlization=$prepData/${cufflinks_run}_cufflinks_utlization.txt
> $cufflinks_utlization
while read work_dir; do
  cd $work_dir/tophat_output
  sample_list=$work_dir/tophat_output/sample_list.txt
  bash ${script_path}/assess_cufflinks_utlization.sh "$sample_list" "$cufflinks_utlization"
done < $horse_trans/working_list_NoPBMCs.txt
############
## relocate the cufflinks analysis results
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cd $work_dir/tophat_output
  for dir in $work_dir/tophat_output/tophat_Barich_Assar_*; do if [ -f "$dir"/transcripts.gtf ]; then
    cd $dir
    mkdir $cufflinks_run && \
    mv "$dir"/{transcripts.gtf,skipped.gtf,*.fpkm_tracking,cufflinks.[oe]*} $cufflinks_run/.; fi; done
fi; done < $horse_trans/working_list_Cerebellum.txt
########################
### Run cuffmerge: merge the sample assemblies and output merged.gtf in tophat_output/$cufflinks_run/$cuffmerge_run
#cuffmerge_run="refGeneGuided_Cuffmerge"
cuffmerge_run="nonGuided_Cuffmerge"
cuffmerge_output=$cufflinks_run/$cuffmerge_run
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cd $work_dir/tophat_output
  for dir in $work_dir/tophat_output/tophat_*; do if [ -f "$dir"/$cufflinks_run/transcripts.gtf ]; then
    echo "$dir"/$cufflinks_run/transcripts.gtf; fi; done > ${cufflinks_run}_assemblies.txt
  rm -fR $cuffmerge_output
  #bash ${script_path}/cuffmerge_withRefGene.sh "$genome" "$cuffmerge_output" "${cufflinks_run}_assemblies.txt" "$refGTF_file"
  bash ${script_path}/cuffmerge_noGTF.sh "$genome" "$cuffmerge_output" "${cufflinks_run}_assemblies.txt"
fi; done < $horse_trans/working_list_Cerebellum.txt
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
  rm -f $tissue_dir/${cufflinks_run}_assemblies.txt
  for dir in $tissue_dir/*; do if [ -d "$dir" ]; then
    cat "$dir"/tophat_output/${cufflinks_run}_assemblies.txt >> $tissue_dir/${cufflinks_run}_assemblies.txt; fi; done
  cuffmerge_output=$tissue/$cufflinks_run/$cuffmerge_run
  rm -fR $cuffmerge_output
  #bash ${script_path}/cuffmerge_withRefGene.sh "$genome" "$cuffmerge_output" "$tissue_dir/${cufflinks_run}_assemblies.txt" "$refGTF_file"
  bash ${script_path}/cuffmerge_noGTF.sh "$genome" "$cuffmerge_output" "$tissue_dir/${cufflinks_run}_assemblies.txt"
done < $horse_trans/multi_lib_tissues.txt
##################
### cuffmerge all assemblies into one total assembly
# save the assembly in $horse_trans/total_merge/cuffmerge/withORwithoutRefGuidence.gtf
cd $tissue_Cuffmerge
rm -f $prepData/${cufflinks_run}_assemblies.txt
while read work_dir; do if [ -d $work_dir/tophat_output ]; then
  cat "$work_dir"/tophat_output/${cufflinks_run}_assemblies.txt >> $prepData/${cufflinks_run}_assemblies.txt
fi; done < $horse_trans/working_list_NoPBMCs.txt

mkdir -p all_tissues
cuffmerge_output=$"all_tissues"/$cufflinks_run/$cuffmerge_run
#bash ${script_path}/cuffmerge_withRefGene.sh "$genome" "$cuffmerge_output" "$prepData/${cufflinks_run}_assemblies.txt" "$refGTF_file"
bash ${script_path}/cuffmerge_noGTF.sh "$genome" "$cuffmerge_output" "$prepData/${cufflinks_run}_assemblies.txt"
###################
## create list of assemblies from each library
## This is where you can edit the list to restrict the processing for certain target(s)
rm -f $prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt
while read work_dir; do
  dir=$work_dir/tophat_output/$cufflinks_run/$cuffmerge_run
  if [ -d $dir ]; then
    echo ${dir#$prepData/} >> $prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## create list of assemblies for tissues of multiple libraries
rm -f $tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt
for tissue in $tissue_Cuffmerge/*/$cufflinks_run/$cuffmerge_run; do if [ -d $tissue ]; then
  echo ${tissue#$tissue_Cuffmerge/} >> $tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt;
fi; done
####################
## convert the gtf files into BigBed files & copy the BigBed files to the track hub directory
rm -f $horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt
while read assembly; do
  echo $assembly
  cd $prepData/$assembly
  targetAss=$"merged.gtf"
  bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  if [ -f $"merged.BigBed" ];then
    identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    cp merged.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    echo $prepData/$assembly >> $horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt;
fi; done < $prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt

while read assembly; do
  echo $assembly
  cd $tissue_Cuffmerge/$assembly
  targetAss=$"merged.gtf"
  bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  if [ -f $"merged.BigBed" ];then
    identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    cp merged.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    echo $tissue_Cuffmerge/$assembly >> $horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt;
fi; done < $tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt

## run icommand to push the file to iplant
## https://pods.iplantcollaborative.org/wiki/display/DS/Using+iCommands
## http://bioinformatics.plantbiology.msu.edu/display/IP/Moving+Data+from+HPCC+to+iPlant
#icd /iplant/home/drtamermansour/horseTrans
#while read assembly; do
#  echo $assembly
#  iput $assembly/*.BigBed
#done < $prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt

## initiate a given track hub for cufflinks_run="refGeneGuided_Cufflinks"
hub_name=$"HorseTrans_TopGuidedCuff"
shortlabel=$"TopGuidedCuff"
longlabel=$"Single samlpe refGuided Tophat/Cufflinks - nonGuided Cuffmerge"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt
tiss_assemblies=$tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies
##########################
## initiate a given track hub for cufflinks_run="nonGuided_Cufflinks"
hub_name=$"HorseTrans_TopNonGuidCuff"
shortlabel=$"TopNonGuidCuff"
longlabel=$"Single samlpe refGuided Tophat - Guided Cufflinks/Cuffmerge"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt
tiss_assemblies=$tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

## add metadata like closest Ref gene
#grep "exon_number \"1\"" merged.gtf > merged_ex1.gtf
#grep "class_code \"u\"" merged_ex1.gtf > merged_ex1_u.gtf
#grep -v "class_code \"u\"" merged_ex1.gtf > merged_ex1_nu.gtf
#######################
## Run Cuffcompare
mkdir $horse_trans/cuffcompare_$shortlabel
cd $horse_trans/cuffcompare_$shortlabel
#sample_list=$horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt
#bash ${script_path}/run_cuffcompare.sh "$sample_list" "$refGTF_file" "$script_path/cuffcompare.sh"
while read assembly; do
  echo $assembly
  identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
  bash ${script_path}/run_cuffcompare2.sh "$prepData/$assembly/merged.gtf" "$identifier" "$refGTF_file" "$script_path/cuffcompare2.sh"
done < $prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt
while read assembly; do
  echo $assembly
  identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
  bash ${script_path}/run_cuffcompare2.sh "$tissue_Cuffmerge/$assembly/merged.gtf" "$identifier" "$refGTF_file" "$script_path/cuffcompare2.sh"
done < $tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt

mkdir $horse_trans/cuffcompare_R_$shortlabel
cd $horse_trans/cuffcompare_R_$shortlabel
while read assembly; do
  echo $assembly
  identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
  bash ${script_path}/run_cuffcompare2.sh "$prepData/$assembly/merged.gtf" "$identifier" "$refGTF_file" "$script_path/cuffcompare3.sh"
done < $prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt
while read assembly; do
  echo $assembly
  identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
  bash ${script_path}/run_cuffcompare2.sh "$tissue_Cuffmerge/$assembly/merged.gtf" "$identifier" "$refGTF_file" "$script_path/cuffcompare3.sh"
done < $tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt
#######################
## compare bed files
mkdir $horse_trans/compareBed
cd $horse_trans/compareBed

## define the list of home made annoatations
alltissueGTF_file=$tissue_Cuffmerge/$"all_tissues"/$cuffmerge_run/merged.gtf

## get a copy of the annotation you want to compare
for f in refGTF_file ncbiNoNameGTF_file ensGTF_file alltissueGTF_file; do
  cp ${!f} ${f}.gtf
  cat ${!f} | gzip > ${f}.gtf.gz
done

## create ExonMerge trackhup
## merge transcripts per loci
#for f in *GTF_file.gtf; do
#  filename=$(basename ${f%.gtf})
#  cat $f \
#    | cgat gtf2gtf --method=sort --sort-order=gene \
#    | cgat gtf2gtf --method=merge-exons --with-utr \
#    | cgat gtf2gtf --method=set-transcript-to-gene \
#    | cgat gtf2gtf --method=sort --sort-order=position \
#    > ${filename}_mergeExons_withUTR.gtf
#done
for f in ncbiNoNameGTF_file ensGTF_file alltissueGTF_file; do
  filename=$(basename ${f%.gtf})
  cat $f \
    | cgat gtf2gtf --method=sort --sort-order=gene \
    | cgat gtf2gtf --method=merge-exons \
    | cgat gtf2gtf --method=set-transcript-to-gene \
    | cgat gtf2gtf --method=sort --sort-order=position \
    > ${filename}_mergeExons.gtf
done
## create list of ExonMerge assemblies
## & convert reference GTF to bed files
rm -f exonMerge_assemblies.txt
for f in *_mergeExons.gtf; do
  bash $script_path/gtfToBigBed.sh "$f" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  identifier=${f%.gtf}
  cp ${identifier}.BigBed $track_hub/$UCSCgenome/BigBed/.
  echo ${identifier} >> exonMerge_assemblies.txt;
done
## initiate a given track hub
hub_name=$"HorseTrans_exonMerge_assemblies"
shortlabel=$"exonMerge_assemblies"
longlabel=$"Assemblies presented with exon merge"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"
## edit the trackDb
> $horse_trans/emptyTemp.txt
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$horse_trans/compareBed/exonMerge_assemblies.txt
tiss_assemblies=$horse_trans/emptyTemp.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

## compare gene sets
python /CAGT_devel/cgat-code/scripts/diff_gtf.py alltissueGTF_file.gtf ncbiNoNameGTF_file.gtf ensGTF_file.gtf > threeway.tsv
python /CAGT_devel/cgat-code/scripts/diff_gtf.py --update=threeway.tsv  alltissueGTF_file.gtf ncbiNoNameGTF_file.gtf ensGTF_file.gtf refGTF_file.gtf > fourway.tsv

mkdir allVSncbi && cd allVSncbi
python /CAGT_devel/cgat-code/scripts/gtfs2tsv.py ../alltissueGTF_file.gtf ../ncbiNoNameGTF_file.gtf > allVSncbi.tsv
mkdir ../allVSens && cd ../allVSens
python /CAGT_devel/cgat-code/scripts/gtfs2tsv.py ../alltissueGTF_file.gtf ../ensGTF_file.gtf > allVSens.tsv
mkdir ../ensVSncbi && cd ../ensVSncbi
python /CAGT_devel/cgat-code/scripts/gtfs2tsv.py ../ensGTF_file.gtf ../ncbiNoNameGTF_file.gtf > ensVSncbi.tsv

## do the bed comparison
module load ucscUtils/262
filename=${alltissueGTF_file%.gtf}
gtfToGenePred $alltissueGTF_file ${filename}.gpred
cat ${filename}.gpred | $script_path/genePredToBed > ${filename}.bed
alltissueBED_file=${filename}.bed

for f in alltissue ncbiNoName ens ref;do
echo $f >> assAnn.txt
bed="$f"BED_file
echo "No_of_transcripts=" $(wc -l ${!bed}) >> assAnn.txt       ## No of transcripts 157567
cat ${!bed} | awk -F $'\t' '{A["Transcripts with "$10" exons= "]++}END{for(i in A)print i,A[i]}' | sort -n > $f.exonsPerTranscript.count
echo "Trans_W_1_exon=" $(cat ${!bed} | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## Trans w 1 exon=  58571
echo "Trans_W_2_exon=" $(cat ${!bed} | awk -F $'\t' '$10 == 2' | wc -l) >> assAnn.txt ## Trans w 2 exons=  8624
echo "Trans_W_>2_exon=" $(cat ${!bed} | awk -F $'\t' '$10 > 1' | wc -l) >> assAnn.txt ## Trans w >2 exons=90372

mergeExon="$f"GTF_file_mergeExons.bed
echo "No_of_genes=" $(wc -l $mergeExon) >> assAnn.txt     ## No of genes 73808
cat $mergeExon | awk -F $'\t' '{A["Gene with "$10" exons= "]++}END{for(i in A)print i,A[i]}' | sort -n > $f.exonsPerGene.count
echo "Genes_W_1_exon=" $(cat $mergeExon | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## genes w 1 exon=54822
echo "Genes_W_2_exon=" $(cat $mergeExon | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## genes w 2 exons=2707
echo "Genes_W_>2_exon=" $(cat $mergeExon | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## genes w >2 exons=90372
done


#cat $alltissueBED_file | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print i,A[i]}' | sort > id_type.count
module load cufflinks/2.2.1
gffread -E $alltissueGTF_file -o ${filename}.gff
alltissueGFF_file=${filename}.gff

module load BEDTools/2.24.0


#######################
## run Transdecoder to predict UTRs with homology options
sample_list="$horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt"
bash $script_path/run_transdecoder.sh $sample_list $genome $refPtn $refPfam $script_path/transdecoder.sh

## calculate the phase of Transdecoder GFF3 files
#while read assembly; do if [ -f $assembly/transdecoder/transcripts.fasta.transdecoder.genome.gff3 ];then
#  echo $assembly
#  cd $assembly/transdecoder
#  bash $script_path/cdsphase.sh transcripts.fasta.transdecoder.genome.gff3
#fi; done < $horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt
#######################
## create list of assemblies from each library
## This is where you can edit the list to restrict the processing for certain target(s)
rm -f $prepData/merged_decoder_assemblies.txt
while read work_dir; do
  for dir in $work_dir/tophat_output/$cufflinks_run/$cuffmerge_run/transdecoder; do if [ -d $dir ]; then
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
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

## pipeline_OneSampleAtaTime_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
bash ${script_path}/TopHatCufflinksCuffmerge_pipline.sh

###########################################################################################
## pipeline_diginormAllsamples_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
bash ${script_path}/DigiTopHatCufflinks_pipline.sh

###########################################################################################
## pipeline_mapped_diginormAllsamples_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
bash ${script_path}/MapDigiTopHatCufflinks_pipline.sh

###########################################################################################
## Other assemblies
mkdir -p $pubAssemblies/Hestand_2014 && cd $pubAssemblies/Hestand_2014
wget http://server1.intrepidbio.com/FeatureBrowser/gtffilereader/record/-4027466309079733025/7666675698.gtf

## create list of public assemblies
rm -f $pubAssemblies/public_assemblies.txt
for tissue in $pubAssemblies/*; do
  echo ${tissue#$pubAssemblies/} >> $pubAssemblies/public_assemblies.txt;
done

####################
## convert the gtf files into BigBed files & copy the BigBed files to the track hub directory
#rm -f $horse_trans/public_assemblies.txt
while read assembly; do
  echo $assembly
  cd $prepData/$assembly
  targetAss=*.gtf
  bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  if [ -f *.BigBed ];then
    identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
    cp *.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    #echo $prepData/$assembly >> $horse_trans/public_assemblies.txt;
fi; done < $pubAssemblies/public_assemblies.txt

## initiate a given track hub
hub_name=$"HorseTrans_public_assemblies"
shortlabel=$"public_assemblies"
longlabel=$"Publically available assemblies"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$pubAssemblies/public_assemblies.txt
#tiss_assemblies=$(> $horse_trans/temp.txt)
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

## create the HTML file page for every track

###########################################################################################


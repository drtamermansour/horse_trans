#!/bin/sh

## construction of the basic diretory structure
## you need to change myRoot to match you actual working directory
myRoot=$"/mnt/ls12/Tamer"
mkdir -p $myRoot/horse_trans/{scripts,resources,rawdata,prepdata,tissue_merge,refGenome,track_hub}

## create a config file to contain all the pathes to be used by all pipelines
echo "horse_trans=$myRoot/horse_trans" >> $myRoot/config.txt
echo "script_path=$myRoot/horse_trans/scripts" >> $myRoot/config.txt
echo "resources=$myRoot/horse_trans/resources" >> $myRoot/config.txt
echo "rawData=$myRoot/horse_trans/rawdata" >> $myRoot/config.txt
echo "prepData=$myRoot/horse_trans/prepdata" >> $myRoot/config.txt
echo "tissue_merge=$myRoot/horse_trans/tissue_merge" >> $myRoot/config.txt
echo "genome_dir=$myRoot/horse_trans/refGenome" >> $myRoot/config.txt
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
## The raw data files should be prepared so that they have enconding "Sanger / illumina 1.9"
## The file names should fit the format *_R1_*.fastq.gz & *_R2_*.fastq.gz for PE reads or *_SR_*.fastq.gz for SE
## All sample replicates should be merged into one sample
## The first syllbus (the part of the file name before _R1_ , _R2_ or _SR_) should be unique
## the final form of the data files should be kept in a folder named fastq_data in the tissue folder

## a representive work out to one SRA reposatories
newLib=$"PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014"
SRA_URL=$"ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653"
mkdir -p $rawData/$newLib
cd $rawData/PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014
bash $script_path/prep_proj.sh $SRA_URL $prepData $script_path  ## The script is missing a step to ensure merging of sample duplicates if they do exist
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

## Add the paths to the common config file
echo "Bowtie2_genome_index_base=$genome_dir/Bowtie2Index/genome" >> $myRoot/config.txt
echo "genome=$genome_dir/Bowtie2Index/genome.fa" >> $myRoot/config.txt
source $myRoot/config.txt
###########################################################################################
## Create GTF file based of refGenes
## generation of GTF from UCSC tables using the guidelines of genomewiki.ucsc
## http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
## The example on the wiki is based ok knowngene table which does not exist for horses. Instead there are refGene and ensGene tables
## The commands of wikigenome are modified to match the table schemes
cd $genome_dir
wget http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/refGene.txt.gz
ucscTable=$"refGene.txt.gz"
output_GTF=$"refGene.gtf"
bash ${script_path}/ucscTableToGTF.sh $ucscTable $output_GTF

## Get the NCBI GFF files
#NCBI_GFF_file=$resources/NCBI/GFF/ref_EquCab2.0_top_level.gff3

## Add the paths to the common config file
echo "refGTF_file=$genome_dir/refGene.gtf" >> $myRoot/config.txt
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

## download pfam database
wget ftp://ftp.broadinstitute.org/pub/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

## Add the paths to the common config file
echo "refPtn=$genome_dir/ptnDB/uniprot_sprot.fasta" >> $myRoot/config.txt
#echo "refPtn=$genome_dir/ptnDB/uniref90.fasta" >> $myRoot/config.txt
echo "refPfam=$genome_dir/ptnDB/Pfam-A.hmm" >> $myRoot/config.txt
source $myRoot/config.txt

## make the databases
bash $script_path/make_ptnDB.sh $refPtn
#bash $script_path/make_ptnDB.sh $refPtn
bash $script_path/make_PfamDB.sh $refPfam
###########################################################################################
## build Tophat transcriptome-index
echo "transcriptome_index=$genome_dir/trans_index/equ" >> $myRoot/config.txt
source $myRoot/config.txt
bash ${script_path}/buildTransIndex.sh "$refGTF_file" "$transcriptome_index" "$Bowtie2_genome_index_base"
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
## Split GTF and BAM files by chromosomes
cat ${refGTF_file} | awk '{ print $1}' | sort | uniq > $genome_dir/refChrom_list.txt
mkdir $genome_dir/GTFbyChrom
while read refChrom; do
  cat ${refGTF_file} | awk -v ref="$refChrom" '($1 == ref)' >> $genome_dir/GTFbyChrom/"$refChrom"_refGene.gtf
done < $genome_dir/refChrom_list.txt

module load BAMTools/2.2.3
bamtools split -in accepted_hits.bam -reference

module load cufflinks/2.2.1
cufflinks --GTF-guide $genome_dir/GTFbyChrom/chr10_refGene.gtf --num-threads 4 --verbose accepted_hits.REF_chr10.bam
###########################################################################################
##### Suggested piplines
### With Tophat
## Compare Ref-guided versus reference free
## Try very senstive Bowtie search
## Try --coverage-search
## For diginorm, try to start with mapped reads only (to minimize error accumulation ==> decrease RAM required to achive minimal false postive rate and allow either final mapping). You can diginorm the non-mapped reads separetly






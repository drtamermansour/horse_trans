SRA_URL="$1"
extData="$2"
script_path="$3"
proj_rawData=$(pwd)

wget -r $SRA_URL/*
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR6534[2-9][0-9]
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR6535[0-9][0-9]
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR6536[0-9][0-9]
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR6537[0-9][0-9]
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR6538[0-9][0-9]
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR6539[0-7][0-9]
#wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR653/ERR653980

#module load SRAToolkit/2.3.4.2
mkdir -p $proj_rawData/temp1
cd $proj_rawData/temp1
for dir in $proj_rawData/$SRA_URL/*; do if [ -d $dir ]; then
        echo $dir
        qsub -v inputdir=${dir} ${script_path}/fastq-dump.sh
        #fastq-dump --split-files --gzip $dir/*.sra
fi done;

## check for failed jobs and re-submet them
mkdir -p dump_reports/failed_reports
mv fastq-dump.* dump_reports/.
cd dump_reports
> failedJobs.txt
> failedSamples.txt
for f in fastq-dump.e*; do
  x=$(du -b $f | cut -f1)
  if [ $x -ne 0 ]; then
    mv $f failed_reports/.
    echo $f >> failedJobs.txt; fi; done;
x=$(cat failedJobs.txt | wc -l)
if [ $x -ne 0 ]; then
  while read e; do
    f=$(echo $e | sed 's/fastq-dump.e/fastq-dump.o/')
    start=$(grep -n "submit_args" $f | cut -d":" -f1)
    end=$(($(grep -n "start_time" $f| cut -d":" -f1)-1))
    # grep the paragraph of submit_args
    sed ''"$start"','"$end"'!d' $f > temp
    # remove the spaces at the start of the lines
    cat temp | sed 's/^[ \t]*//' > temp2
    # remove the new line characters
    sed ':a;N;$!ba;s/\n//g' temp2 > temp3
    # grap the file name
    cat temp3 | cut -d" " -f4 | cut -d"=" -f2
    #sed ''"$start"','"$end"'!d' $f | sed 's/^[ \t]*//' | sed ':a;N;$!ba;s/\n//g' | cut -d" " -f4 | cut -d"=" -f2
  done < failedJobs.txt > failedSamples.txt
  cd $proj_rawData/temp1
  mkdir -p failedSamples
  while read dir; do
    if [ -d $dir ]; then
      echo $dir
      mv "$(basename $dir)"_* failedSamples/.
      qsub -v inputdir=${dir} ${script_path}/fastq-dump.sh
      #fastq-dump --split-files --gzip $dir/*.sra
    fi; done; < dump_reports/failedSamples.txt
fi
## check again until you make sure that there are no failed samples

## check encoding
#cd fastq_data/ERR653420
#gunzip -c *_1.fastq.gz | awk 'NR % 4 == 0' | head -n 1000000 | python $script_path/guess-encoding.py

## change the name of the data files to match the standard format
cd $proj_rawData/temp1
for R1 in *_1.fastq.gz; do
echo $R1; R2=$(echo $R1 | sed s/_1.fastq.gz/_2.fastq.gz/);
newR1=$(echo $R1 | sed s/_1.fastq.gz/_R1_001.fastq.gz/);
newR2=$(echo $R1 | sed s/_1.fastq.gz/_R2_001.fastq.gz/);
mv $R1 $newR1; mv $R2 $newR2; done;

bash ${script_path}/run_fastqc.sh "$proj_rawData/temp1"

mkdir -p $extData/PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014/fastq_data
mv $rawData/PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014/temp1/* $extData/PBMCs/PE_49_fr.unstranded_bioproj.265983_10302014/fastq_data/.




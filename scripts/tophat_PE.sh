#!/bin/bash -login
#PBS -l walltime=30:00:00,nodes=1:ppn=4,mem=16Gb
#mdiag -A ged
#PBS -m abe			
#PBS -N tophat		

module load TopHat2/2.0.14

#output_dir="$1"
#Bowtie2_genome_index_base="$2"
#transcriptome_index="$3"
#lib_type="$4"
#input_file1="$5"
#input_file2="$6"

cd $PBS_O_WORKDIR
tophat \
--transcriptome-index "${transcriptome_index}" \
--num-threads 4 \
--output-dir "${output_dir}" \
--microexon-search \
--library-type "${lib_type}" \
"${Bowtie2_genome_index_base}" $f1 $f2

qstat -f ${PBS_JOBID}


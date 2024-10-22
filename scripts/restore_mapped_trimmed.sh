#!/bin/bash -login
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=16Gb
#mdiag -A ged
#PBS -m abe
#PBS -N restore_mapped

module load SAMTools/1.0

#sample="$1"
#label="$2"
#lib_type="$3"
#output_dir="$4"

cd $output_dir

output_pe=$label"_R_001.pe.fq"
output_se=$label"_R_001.se.fq"
# 1. Shuffle the reads in the bam file
samtools bamshuf -uOn 128 $sample tmp_$label > "$label"_shuf_accepted_hits.bam
# 2. Revert the BAM file to FastQ
samtools bam2fq "$label"_shuf_accepted_hits.bam -s $output_se > $output_pe
# 3. delete the empty pe output in SE libraries
if [ "$lib_type" = $"SE" ]; then rm $output_pe; fi


qstat -f ${PBS_JOBID}


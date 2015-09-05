#!/bin/bash -login
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=16Gb
#mdiag -A ged
#PBS -m abe
#PBS -N restore_mapped

module load SAMTools/1.0

#input_dir="$1"
#sample="$2"
#output_pe="$3"
#output_se="$4"
#lib_type="$5"


cd $PBS_O_WORKDIR

# 1. Shuffle the reads in the bam file
samtools bamshuf -uOn 128 $input_dir/accepted_hits.bam tmp_$sample > "$sample"_shuf_accepted_hits.bam
# 2. Revert the BAM file to FastQ
samtools bam2fq "$sample"_shuf_accepted_hits.bam -s $output_se > $output_pe
# 3. delete the empty pe output in SE libraries
if [ "$lib_type" = $"SE" ]; then rm $output_pe; fi


qstat -f ${PBS_JOBID}


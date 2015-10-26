#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N MarkDuplicates

module load picardTools/1.113

cd $PBS_O_WORKDIR

java -Xmx10g -jar $PICARD/MarkDuplicates.jar INPUT=$sample OUTPUT="dedup_reads.bam" METRICS_FILE=metrics.txt

qstat -f ${PBS_JOBID}



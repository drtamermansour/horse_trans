#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=2,mem=10Gb
#mdiag -A ged
#PBS -m abe
#PBS -N applyRecalib

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar \
-T PrintReads \
-R $gatk_ref \
-I $sample \
-BQSR ../recal_data.table \
-o recal_reads.bam

qstat -f ${PBS_JOBID}
#!/bin/bash -login
#PBS -l walltime=48:00:00,nodes=1:ppn=2,mem=80Gb
#mdiag -A ged
#PBS -m abe
#PBS -N indelRealigner

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx78g -jar $GATK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $gatk_ref \
-I $sample \
-targetIntervals $intervals \
-nWayOut '.realigned.bam' \
-known $indels \
-model USE_SW \
-LOD 0.4 \
-dcov 500

qstat -f ${PBS_JOBID}
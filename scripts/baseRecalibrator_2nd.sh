#!/bin/bash -login
#PBS -l walltime=12:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N baseRecalibrator-2nd

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $gatk_ref \
-I $sample \
-L 20 \
-known $snps \
-known $indels \
-BQSR recal_data.table \
-o post_recal_data.table

qstat -f ${PBS_JOBID}




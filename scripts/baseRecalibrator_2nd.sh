#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=2,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N baseRecalibrator-2nd

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $gatk_ref \
$(echo $samples) \
$(echo $variants) \
-BQSR recal_data.table \
-o post_recal_data.table

qstat -f ${PBS_JOBID}




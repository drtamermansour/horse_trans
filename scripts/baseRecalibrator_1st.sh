#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=2,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N baseRecalibrator-1stR

module load GATK/3.4.46

cd $PBS_O_WORKDIR

#java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
#-T BaseRecalibrator \
#-R $gatk_ref \
#-I $sample \
#$(echo $variants) \
#-o recal_data.table

java -Xmx42g -jar $GATK/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $gatk_ref \
$(echo $samples) \
$(echo $variants) \
-o recal_data.table


qstat -f ${PBS_JOBID}




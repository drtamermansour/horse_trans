#!/bin/bash -login
#PBS -l walltime=1:00:00,nodes=1:ppn=2,mem=4Gb
#mdiag -A ged
#PBS -m abe
#PBS -N createRecalPlots

module load GATK/3.4.46
module load R/3.0.1

cd $PBS_O_WORKDIR

java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R $gatk_ref \
-before recal_data.table \
-after post_recal_data.table \
-plots recalibration_plots.pdf

qstat -f ${PBS_JOBID}
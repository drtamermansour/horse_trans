#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N combineVariants

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $gatk_ref \
$(echo $samples) \
-o HC_output_ploidy1_haplo1.vcf \
-genotypeMergeOptions UNIQUIFY

qstat -f ${PBS_JOBID}











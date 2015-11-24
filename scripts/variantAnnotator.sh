#!/bin/bash -login
#PBS -l walltime=10:00:00,nodes=1:ppn=8,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N haplotypeCaller_multi

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx45g -jar $GATK/GenomeAnalysisTK.jar \
-R $gatk_ref \
-T VariantAnnotator \
$(echo $samples) \
-V HC_output_ploidy1_haplo1.vcf \
-A Coverage \
-A FisherStrand \
-A QualByDepth \
-nt 7 \
-o HC_output_ploidy1_haplo1_ann.vcf


qstat -f ${PBS_JOBID}





















#!/bin/bash -login
#PBS -l walltime=12:00:00,nodes=1:ppn=4,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N haplotypeCaller_multi

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $gatk_ref \
-I $sample \
--emitRefConfidence GVCF \
--dbsnp $snps \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-nct 3 \
-o $output ##samplename.vcf ##HC_output_ploidy1_haplo1.vcf

qstat -f ${PBS_JOBID}

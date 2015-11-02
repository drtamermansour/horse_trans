#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N haplotypeCaller

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $gatk_ref \
-I $sample \
--dbsnp $snps \
-dontUseSoftClippedBases \
-stand_call_conf 20.0 \
-stand_emit_conf 20.0 \
--minPruning 10 \
-ploidy 1 \
-maxNumHaplotypesInPopulation 1 \
-o HC_output.vcf

qstat -f ${PBS_JOBID}

## best practice for RNAseq
#java -jar GenomeAnalysisTK.jar \
#-T HaplotypeCaller \
#-R ref.fasta \
#-I input.bam \
#-dontUseSoftClippedBases \
#-stand_call_conf 20.0 \
#-stand_emit_conf 20.0 \
#-o output.vcf

## my natal teeth project
#java -Xmx10g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar \
#    -T HaplotypeCaller \
#    -R $ref \
#    -I recal_G13.bam \
#    --emitRefConfidence GVCF \
#    --variant_index_type LINEAR \
#    --variant_index_parameter 128000 \
#    --dbsnp $knownIndels_dir/dbsnp_138.b37.vcf \
#    -o Gr13.raw.snps.indels.g.vcf





#-minPruning (Default=2): The basic idea is that sections of the graph that are supported by very few reads are most probably the result of stochastic errors, so we are going to remove any sections that are supported by fewer than a certain threshold number of reads. increase to 10

#-maxNumHaplotypesInPopulation (Default=128): limit the number of haplotypes that will be considered for each value of k choosing the most likely haplotypes

#-ploidy (Default=2): Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).

















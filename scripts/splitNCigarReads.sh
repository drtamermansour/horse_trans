#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N splitNCigarReads

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R $gatk_ref -I $sample -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

qstat -f ${PBS_JOBID}

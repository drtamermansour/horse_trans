#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N realignerTargetCreator

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R "$gatk_ref" \
$(echo $samples) \
-o gatk.intervals \
-known "$indels"

qstat -f ${PBS_JOBID}
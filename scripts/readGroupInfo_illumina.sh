#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N AddOrReplaceReadGroups

module load picardTools/1.113

cd $PBS_O_WORKDIR

java -Xmx10g -jar $PICARD/AddOrReplaceReadGroups.jar INPUT=$bam OUTPUT=$output RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM

qstat -f ${PBS_JOBID}



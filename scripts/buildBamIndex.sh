#!/bin/bash -login
#PBS -l walltime=00:30:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N buildBamIndex

module load picardTools/1.113

cd $PBS_O_WORKDIR

java -Xmx10g -jar $PICARD/BuildBamIndex.jar INPUT=$sample

qstat -f ${PBS_JOBID}

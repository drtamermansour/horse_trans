#!/bin/bash -login
#PBS -l walltime=01:00:00,nodes=1:ppn=2,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N reorderBAM

module load picardTools/1.113

cd $PBS_O_WORKDIR

java -Xmx10g -jar $PICARD/ReorderSam.jar I=$sample R=$gatk_ref O=reordered.bam

qstat -f ${PBS_JOBID}

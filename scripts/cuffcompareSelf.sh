#!/bin/bash -login
#PBS -l walltime=00:30:00,nodes=1:ppn=1,mem=4Gb
#mdiag -A ged
#PBS -m abe
#PBS -N cuffcompare

module load cufflinks/2.2.1


#sample="$1"
#label="$2"
#log="$3"

#cd $PBS_O_WORKDIR
cd $(dirname $sample)

cuffcompare -T -V -o ${label} ${sample} &> $log

qstat -f ${PBS_JOBID}



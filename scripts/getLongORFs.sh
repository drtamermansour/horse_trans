#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=1,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N getLongOrfs

module load TransDecoder/2.0.1

cd $PBS_O_WORKDIR

# extract the long open reading frames
TransDecoder.LongOrfs -t $inputFASTA


qstat -f ${PBS_JOBID}

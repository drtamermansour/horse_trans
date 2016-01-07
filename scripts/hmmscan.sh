#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=4,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N hmmscan

module load HMMER/3.0

cd $PBS_O_WORKDIR

hmmscan --cpu 4 --domtblout $output $refPfam $pep


qstat -f ${PBS_JOBID}

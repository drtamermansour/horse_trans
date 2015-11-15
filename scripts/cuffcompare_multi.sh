#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N cuffcompare

module load cufflinks/2.2.1


#Genes_GTF_file="$1"
#samples="$2"

cd $PBS_O_WORKDIR

cuffcompare -V -r ${Genes_GTF_file} ${samples[@]}

qstat -f ${PBS_JOBID}



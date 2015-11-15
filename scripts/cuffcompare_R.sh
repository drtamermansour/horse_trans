#!/bin/bash -login
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=8Gb
#mdiag -A ged
#PBS -m abe
#PBS -N cuffcompare

module load cufflinks/2.2.1


#Genes_GTF_file="$1"
#sample="$2"
#label="$3"

cd $PBS_O_WORKDIR

cuffcompare -V -R -r ${Genes_GTF_file} -o ${label} ${sample}

qstat -f ${PBS_JOBID}



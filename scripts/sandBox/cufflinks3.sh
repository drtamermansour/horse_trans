#!/bin/bash -login
#PBS -l walltime=7:00:00:00,nodes=1:ppn=32,mem=180Gb
#mdiag -A ged
#PBS -m abe
#PBS -N cufflinks

module load cufflinks/2.2.1


#Genes_GTF_file="$1"
#label="$2"
#sample="$3"

cd $PBS_O_WORKDIR

cufflinks --GTF-guide ${Genes_GTF_file} --label ${label} --num-threads 32 --verbose ${sample}

qstat -f ${PBS_JOBID}

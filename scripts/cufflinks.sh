#!/bin/bash -login
#PBS -l walltime=72:00:00,nodes=1:ppn=8,mem=64Gb
#mdiag -A ged	
#PBS -m abe			
#PBS -N cufflinks		

module load cufflinks/2.2.1


#Genes_GTF_file="$1"
#label="$2" 

cd $PBS_O_WORKDIR

cufflinks --GTF-guide ${Genes_GTF_file} --label ${label} --num-threads 8 --verbose accepted_hits.bam

qstat -f ${PBS_JOBID}



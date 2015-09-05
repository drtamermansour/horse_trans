#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N fastq-dump		#give name to the job


module load SRAToolkit/2.3.4.2

cd $PBS_O_WORKDIR

inputdir=$inputdir

fastq-dump --split-files --gzip ${inputdir}/*.sra

qstat -f ${PBS_JOBID}




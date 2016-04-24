#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N tabix		#give name to the job

module load tabix/0.2.6

cd $PBS_O_WORKDIR

tabix -p vcf $f

qstat -f ${PBS_JOBID}



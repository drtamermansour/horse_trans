#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N BLAT		#give name to the job


module load BLAT/35

cd $PBS_O_WORKDIR

blat -t=dna $genome -q=rna $transcriptome -ooc=$oocFile -fine -noHead $output

qstat -f ${PBS_JOBID}


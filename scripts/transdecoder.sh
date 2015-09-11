#!/bin/bash -login
#PBS -l walltime=2:00:00:00,nodes=1:ppn=8,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N TransDecoder

module load TransDecoder/2.0.1
module load BLAST+/2.2.30
module load HMMER/3.0
decoderUtil=$"/opt/software/TransDecoder/2.0.1--GCC-4.4.5/util"

cd $PBS_O_WORKDIR

$decoderUtil/cufflinks_gtf_genome_to_cdna_fasta.pl merged.gtf $genome > transcripts.fasta
$decoderUtil/cufflinks_gtf_to_alignment_gff3.pl merged.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep  -db $refPtn  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6
hmmscan --cpu 8 --domtblout pfam.domtblout $refPfam transcripts.fasta.transdecoder_dir/longest_orfs.pep

qstat -f ${PBS_JOBID}

ssembly="$tissue_Cuffmerge/$cuffmerge_output/filtered"
cd $assembly/varFixed
qsub -v f="GenotypeGVCFs_monoAllel_sig_transFinal.vcf" ${script_path}/sandBox/bgzipFiles.sh
#tabix -p vcf GenotypeGVCFs_output_max50.raw_SNPs.vcf.gz
qsub -v f="GenotypeGVCFs_monoAllel_sig_transFinal.vcf.gz" ${script_path}/sandBox/tabixFiles.sh


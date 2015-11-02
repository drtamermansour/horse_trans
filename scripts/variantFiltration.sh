gatk_ref="$1"
input="$2"

module load GATK/3.4.46


output=${input%.vcf}_filtered.vcf

java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $gatk_ref \
-V $input \
-window 35 -cluster 3 \
-filterName FS -filter "FS > 30.0" \
-filterName QD -filter "QD < 2.0" \
-o $output




## -window 35 -cluster 3: filter clusters of at least 3 SNPs that are within a window of 35 bases between them
## FS: filtering based on Fisher Strand values (FS > 30.0)
## QD: Qual By Depth values (QD < 2.0).











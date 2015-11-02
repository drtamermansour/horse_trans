gatk_ref="$1"
input="$2"
output="$3"

module load GATK/3.4.46

java -Xmx1g -jar $GATK/GenomeAnalysisTK.jar \
-T FastaAlternateReferenceMaker \
-R $gatk_ref \
-V $input \
-o $output


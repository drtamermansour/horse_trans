gatk_ref="$1"
input="$2"
chain="$3"
dict="$4"
output="$5"

module load GATK/3.4.46

java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar \
-T LiftoverVariants \
-R "$gatk_ref" \
-V "$input" \
-chain "$chain" \
-dict "$dict" \
-o "$output"

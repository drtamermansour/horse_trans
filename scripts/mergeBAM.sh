outputName="$1"
bamList=( "$@" )
unset bamList[0]

echo "${bamList[*]}"

#module load SAMTools/0.1.19
#samtools merge "$outputName" "${bamList[*]}"
#Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users must provide the correct header with -h, or uses Picard which properly maintains the header dictionary in merging.


module load picardTools/1.113

inputs=()
for bam in ${bamList[*]}; do
  input=$"INPUT="$bam
  inputs+=($input)
done

java -Xmx10g -jar $PICARD/MergeSamFiles.jar ${inputs[*]} OUTPUT=$outputName


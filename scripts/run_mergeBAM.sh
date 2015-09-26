replicates_list="$1"
script="$2"

mkdir -p singleSamples
pathToSingleSamples=$(pwd)/singleSamples
while read sample; do
  replicates=($(echo $sample))
  echo ${replicates[@]}
  len=${#replicates[@]}
  if [ $len -gt 1 ]; then
    bamList=()
    newName=""
    for rep in "${replicates[@]}"; do
      base=${rep%_R1_*.fastq.gz}
      mv tophat_$base singleSamples/.
      bam=$pathToSingleSamples/tophat_$base/accepted_hits_RG.bam
      if [ ! -f "$bam" ];then
        bam=$pathToSingleSamples/tophat_$base/accepted_hits.bam; fi
      bamList+=($bam)
      newName=$newName$base"_"
    done
    mkdir tophat_${newName%_}
    cd tophat_${newName%_}
    bash $script "$(basename $bam)" "${bamList[*]}"
    cd ../
fi; done < $replicates_list
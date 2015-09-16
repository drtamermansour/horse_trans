sample_list="$1"

## Check for successful trimming
> $sample_list
success=$"Processed"
mkdir -p failed_reports
f=cufflinks.e*
if [ -f $f ]; then
  line=$(tail -1 $f | cut -d" " -f1)
  if [ "$line" != "$success" ]; then
    echo $(pwd)/merged.bam >> $sample_list
    mv $f failed_reports/.
  else echo "no cufflinks reports in $dir"; fi
fi
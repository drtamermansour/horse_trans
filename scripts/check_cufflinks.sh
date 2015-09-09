sample_list="$1"

## Check for successful trimming
> $sample_list
success=$"Processed"
for dir in tophat_*; do if [ -d $dir ]; then
  mkdir -p $dir/failed_reports
  f=$dir/cufflinks.e*
  if [ -f $f ]; then
    line=$(tail -1 $f | cut -d" " -f1)
    if [ "$line" != "$success" ]; then
      echo $(pwd)/$dir >> $sample_list
      mv $f $dir/failed_reports/.; fi
  else echo "no cufflinks reports in $dir"; fi
fi; done;
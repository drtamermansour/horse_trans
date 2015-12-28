sample_list="$1"
failedSample_list="$2"

for dir in tophat_*; do if [ -d $dir ]; then
  echo $dir
  f=$dir/logs/tophat.log
  if [ -f $f ]; then
    x=$(grep -c "Run complete" $f)
    if [ $x -eq 0 ]; then
      key=${dir#tophat_}
      grep $key $sample_list >> $failedSample_list
      echo "Sample failed"
    else echo "Run complete"; fi
  else echo "can't find tophat.log file"; fi
fi; done
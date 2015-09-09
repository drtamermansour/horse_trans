 sample_list="$1"

## Check for successful trimming
mkdir -p failed_reports
success=$"fp rate"
f=$"digi_norm.e*"
if [ -f $f ]; then
  line=$(tail -1 $f | cut -d" " -f1,2)
  if [ "$line" != "$success" ]; then
    dirname $(pwd) >> $sample_list
    mv $f failed_reports/.; fi
else echo "no diginorm reports in $(pwd)"; fi

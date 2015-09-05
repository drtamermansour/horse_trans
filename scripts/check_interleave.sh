sample_list="$1"

## Check for successful interleave
mkdir -p job_reports/failed_reports
mv interleave.* job_reports/. 2>/dev/null
cd job_reports

> $sample_list
success=$"output written to block device"
for f in interleave.e*; do if [ -f $f ]; then
  line=$(tail -1 $f)
  if [ "$line" != "$success" ]; then
    grep -A1 "Interleaving:" $f | tail -1 | sed 's/^[ \t]*//' | sed 's/.pe.fq/.pe.se.fq/' >> $sample_list
    mv $f failed_reports/.
fi; fi; done;
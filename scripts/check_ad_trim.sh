sample_list="$1"

## Check for successful trimming
mkdir -p job_reports/failed_reports
mv T_Trim.* job_reports/. 2>/dev/null
cd job_reports
> $sample_list
for f in T_Trim.e*; do if [ -f $f ]; then
  x=$(grep -c "Completed successfully" $f) #line=$(tail -1 $f)
  if [ $x -eq 0 ]; then
    grep "Started with arguments" $f | cut -d" " -f 8 >> $sample_list
    mv $f failed_reports/.
  fi; fi; done;


## Some alteranative code to Check for successful trimming
#cd $work_dir/trimmed_RNA_reads
#for f in *_R1_*.pe.fq; do
#f2=$(echo $f | sed 's/_R1_/_R2_/')
#x1=$(wc -l $f | cut -d" " -f1)
#x2=$(wc -l $f2 | cut -d" " -f1)
#if [ $x1 != $x2 ]; then
#echo $f >> failedPairing; fi
#z1=$(tail -n 4 $f | grep -n "^@" | cut -d":" -f1)
#z2=$(tail -n 4 $f2 | grep -n "^@" | cut -d":" -f1)
#if [ z1 -ne 0 || z2 -ne 0 ]; then
#echo $f >> failedPairing2; fi
#g1=$(tail -n 4 $f | grep -v "^[@+]" | head -1 | awk '{ print length($0); }')
#g2=$(tail -n 4 $f | grep -v "^[@+]" | tail -1 | awk '{ print length($0); }')
#if [ q1 -eq g2 ]; then
#echo $f >> failedPairing3; fi
#done


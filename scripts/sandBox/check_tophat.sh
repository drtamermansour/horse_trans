sample_list="$1"

## Check for successful trimming
mkdir -p job_reports/failed_reports
mv tophat.[eo]* job_reports/. 2>/dev/null
cd job_reports
for f in tophat.e*; do if [ -f $f ]; then
  x=$(grep -c "Run complete" $f) #line=$(tail -1 $f)
  if [ $x -eq 0 ]; then
    f2=$(echo $f | sed 's/.e/.o/')
    start=$(grep -n "submit_args" $f2 | cut -d":" -f1)
    end=$(($(grep -n "start_time" $f2 | cut -d":" -f1)-1))
    # grep the paragraph of submit_args
    sed ''"$start"','"$end"'!d' $f2 > temp1
    # remove the spaces at the start of the lines
    cat temp1 | sed 's/^[ \t]*//' > temp2
    # remove the new line characters
    sed ':a;N;$!ba;s/\n//g' temp2 > temp3
    # grap the file name
    cat temp3 | cut -d" " -f4 | cut -d"," -f5 | cut -d"=" -f2 >> $sample_list
    #sed ''"$start"','"$end"'!d' $f | sed 's/^[ \t]*//' | sed ':a;N;$!ba;s/\n//g' | cut -d" " -f4 | cut -d"=" -f2
    mv $f failed_reports/.
  fi; fi; done;
rm -f temp[1-3]